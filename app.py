import streamlit as st
import pandas as pd
import sqlite3
import gzip
import os

st.title("🧬 VCF Variant Annotation Tool")

# =============================
# LOAD DATABASE
# =============================
@st.cache_data
def load_database():
    try:
        df = pd.read_csv("annotated_table_updated.csv")
    except:
        st.error("❌ annotated_table_updated.csv missing")
        st.stop()

    df.columns = df.columns.astype(str)
    df.columns = df.columns.str.strip()
    df.columns = df.columns.str.replace(r"\s+", "_", regex=True)

    if "SNP_ID" in df.columns:
        df = df.rename(columns={"SNP_ID": "snp_id"})
    if "P-value" in df.columns:
        df = df.rename(columns={"P-value": "p_value"})
    if "Phenotype_Description" in df.columns:
        df = df.rename(columns={"Phenotype_Description": "phenotype"})

    required = ["snp_id", "Sub_Trait", "p_value"]
    for col in required:
        if col not in df.columns:
            st.error(f"❌ Missing column: {col}")
            st.write(df.columns)
            st.stop()

    df["snp_id"] = df["snp_id"].astype(str).str.strip()
    return df

db_df = load_database()

st.subheader("📊 Database Preview")
st.dataframe(db_df.head())

# =============================
# PARSE VCF
# =============================
def parse_vcf(file):
    data = []

    for line in file:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        if len(cols) < 10:
            continue

        try:
            snp_id = cols[2].strip()
            chromosome = cols[0]
            position = int(cols[1])
            ref = cols[3]
            alt = cols[4]

            gt_raw = cols[9].split(":")[0]

            if gt_raw == "0/0":
                genotype = "0|0"
            elif gt_raw in ["0/1", "1/0"]:
                genotype = "0|1"
            elif gt_raw == "1/1":
                genotype = "1|1"
            else:
                genotype = "NA"

            data.append([snp_id, chromosome, position, ref, alt, genotype])
        except:
            continue

    return pd.DataFrame(data, columns=[
        "snp_id", "chromosome", "position",
        "ref", "alt", "genotype"
    ])

# =============================
# LOAD INTERNAL REFERENCE (SAFE)
# =============================
def load_reference():
    try:
        if not os.path.exists("reference.vcf.gz"):
            return pd.DataFrame()

        with gzip.open("reference.vcf.gz", "rt") as f:
            return parse_vcf(f)

    except:
        return pd.DataFrame()

# =============================
# MATCHING (SAFE SQL)
# =============================
def annotate(vcf_df):
    try:
        conn = sqlite3.connect(":memory:")

        vcf_df.to_sql("vcf", conn, index=False, if_exists="replace")
        db_df.to_sql("db", conn, index=False, if_exists="replace")

        query = """
        SELECT 
            v.snp_id,
            v.chromosome,
            v.position,
            v.ref,
            v.alt,
            v.genotype,
            d.Gene,
            d.Sub_Trait,
            d.p_value,
            d.phenotype
        FROM vcf v
        JOIN db d ON v.snp_id = d.snp_id
        """

        result = pd.read_sql_query(query, conn)
        conn.close()

        return result

    except Exception as e:
        st.error(f"Annotation error: {e}")
        return pd.DataFrame()

# =============================
# SAS SCORE (SAFE)
# =============================
def sas_score(df):
    if "SAS_MAF" not in df.columns or len(df) == 0:
        return 0
    return df["SAS_MAF"].mean() * 100

# =============================
# IBS CALCULATION
# =============================
def calculate_ibs(df):
    score, count = 0, 0

    for _, row in df.iterrows():
        g1 = str(row["genotype_user"])
        g2 = str(row["genotype_ref"])

        if g1 == g2:
            score += 2
        elif ("0|1" in g1) or ("0|1" in g2):
            score += 1

        count += 1

    return (score / (2 * count)) * 100 if count else 0

# =============================
# TOP SNP
# =============================
def get_top_snps(group):
    return group.head(10) if len(group) >= 10 else group.head(5)

# =============================
# INPUT VCF (UPLOAD ONLY - STABLE)
# =============================
st.subheader("📂 Upload VCF")

uploaded_vcf = st.file_uploader("Upload VCF file", type=["vcf"])

# =============================
# RUN PIPELINE
# =============================
if uploaded_vcf is not None:

    st.write("🔄 Processing...")

    try:
        vcf_df = parse_vcf(uploaded_vcf)

        if len(vcf_df) == 0:
            st.error("❌ No valid SNPs found in VCF")
            st.stop()

        matched_df = annotate(vcf_df)

    except Exception as e:
        st.error(f"❌ Processing failed: {e}")
        st.stop()

    st.success("✅ Done")

    st.write("📊 Total SNPs:", len(vcf_df))
    st.write("✅ Matched SNPs:", len(matched_df))

    # =============================
    # SAS SCORE
    # =============================
    st.metric("🧬 SAS Score (%)", f"{sas_score(matched_df):.2f}")

    # =============================
    # IBS (AUTO INTERNAL)
    # =============================
    st.subheader("🧬 IBS Similarity")

    reference_df = load_reference()

    if not reference_df.empty:

        comp = pd.merge(
            vcf_df,
            reference_df,
            on="snp_id",
            suffixes=("_user", "_ref")
        )

        st.write("🔗 Common SNPs:", len(comp))

        if len(comp) > 0:
            st.metric("IBS (%)", f"{calculate_ibs(comp):.2f}")
        else:
            st.warning("No overlapping SNPs")

    else:
        st.warning("⚠️ reference.vcf.gz not found (IBS disabled)")

    # =============================
    # SHOW MATCHED
    # =============================
    st.subheader("🔍 Matched SNPs")
    st.dataframe(matched_df)

    # =============================
    # TOP SNPs
    # =============================
    st.subheader("🏆 Top SNPs")

    if len(matched_df) > 0:
        df_top = matched_df.copy()

        df_top["p_value"] = pd.to_numeric(df_top["p_value"], errors="coerce")
        df_top = df_top.dropna(subset=["Sub_Trait", "p_value"])
        df_top = df_top.sort_values("p_value")

        top = df_top.groupby("Sub_Trait", group_keys=False).apply(get_top_snps)

        st.dataframe(top)

    # =============================
    # DOWNLOAD
    # =============================
    st.download_button(
        "⬇ Download Results",
        matched_df.to_csv(index=False),
        "results.csv"
    )

# =============================
# SEARCH
# =============================
st.subheader("🔍 Search SNP")

search = st.text_input("Enter rsID")

if search:
    res = db_df[db_df["snp_id"].str.contains(search, case=False, na=False)]

    if len(res) > 0:
        st.dataframe(res)
    else:
        st.warning("No match found")
