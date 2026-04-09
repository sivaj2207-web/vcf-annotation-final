import streamlit as st
import pandas as pd
import sqlite3
import gzip

# =============================
# TITLE
# =============================
st.title("🧬 VCF Variant Annotation Tool")

# =============================
# LOAD DATABASE
# =============================
@st.cache_data
def load_database():
    df = pd.read_csv("annotated_table_updated.csv")

    df.columns = df.columns.astype(str)
    df.columns = df.columns.str.strip()
    df.columns = df.columns.str.replace(r"\s+", "_", regex=True)

    if "SNP_ID" in df.columns:
        df = df.rename(columns={"SNP_ID": "snp_id"})
    if "P-value" in df.columns:
        df = df.rename(columns={"P-value": "p_value"})
    if "Phenotype_Description" in df.columns:
        df = df.rename(columns={"Phenotype_Description": "phenotype"})

    for col in df.columns:
        if "SAS" in col.upper():
            df = df.rename(columns={col: "SAS_MAF"})

    df["snp_id"] = df["snp_id"].astype(str).str.strip()
    return df

db_df = load_database()

# =============================
# PREVIEW
# =============================
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
            if snp_id == ".":
                continue

            chromosome = cols[0]
            position = int(cols[1])
            ref = cols[3]
            alt = cols[4]

            gt = cols[9].split(":")[0]

            if gt == "0/0":
                genotype = "0|0"
            elif gt in ["0/1", "1/0"]:
                genotype = "0|1"
            elif gt == "1/1":
                genotype = "1|1"
            else:
                genotype = "NA"

            data.append([snp_id, chromosome, position, ref, alt, genotype])
        except:
            continue

    df = pd.DataFrame(data, columns=[
        "snp_id", "chromosome", "position",
        "ref", "alt", "genotype"
    ])

    df["snp_id"] = df["snp_id"].str.lower().str.strip()
    return df

# =============================
# INPUT
# =============================
st.subheader("📂 Input VCF")

uploaded_vcf = st.file_uploader("Upload VCF file", type=["vcf"])

vcf_df = None
if uploaded_vcf:
    vcf_df = parse_vcf(uploaded_vcf)

# =============================
# MATCHING
# =============================
def annotate(vcf_df):
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
        d.phenotype,
        d.SAS_MAF
    FROM vcf v
    JOIN db d
    ON v.snp_id = d.snp_id
    """

    result = pd.read_sql_query(query, conn)
    conn.close()
    return result

# =============================
# SAS
# =============================
def sas_score(df):
    if "SAS_MAF" in df.columns and len(df) > 0:
        df["SAS_MAF"] = pd.to_numeric(df["SAS_MAF"], errors="coerce")
        return df["SAS_MAF"].mean() * 100
    return 0

# =============================
# IBS
# =============================
def load_reference():
    try:
        with gzip.open("final.vcf.gz", "rt") as f:
            return parse_vcf(f).head(3000)  # limit reference
    except:
        return pd.DataFrame()

def calculate_ibs(df):
    if len(df) == 0:
        return 0

    match = (df["genotype_user"] == df["genotype_ref"]).sum()
    het = ((df["genotype_user"] == "0|1") | (df["genotype_ref"] == "0|1")).sum()

    return ((2 * match + het) / (2 * len(df))) * 100

# =============================
# TOP SNP
# =============================
def get_top_snps(group):
    return group.head(10) if len(group) >= 10 else group.head(5)

# =============================
# RUN
# =============================
if vcf_df is not None:

    st.write("🔄 Processing...")

    try:
        matched_df = annotate(vcf_df)   # ✅ FULL VCF used here
    except Exception as e:
        st.error(f"❌ {e}")
        st.stop()

    st.success("✅ Done")

    st.write("📊 Total SNPs:", len(vcf_df))
    st.write("✅ Matched SNPs:", len(matched_df))

    # SAS
    st.metric("🧬 SAS Score (%)", f"{sas_score(matched_df):.2f}")

    # =============================
    # IBS (ONLY SAMPLE USED)
    # =============================
    st.subheader("🧬 IBS Similarity")

    vcf_sample = vcf_df.head(3000)  # 🔥 only for IBS
    ref_df = load_reference()

    if not ref_df.empty:
        comp = pd.merge(
            vcf_sample,
            ref_df,
            on="snp_id",
            suffixes=("_user", "_ref")
        ).head(2000)

        st.write("🔗 Common SNPs:", len(comp))

        if len(comp) > 0:
            ibs = calculate_ibs(comp)
            st.metric("🧬 IBS (%)", f"{ibs:.2f}")
        else:
            st.warning("No overlapping SNPs")
    else:
        st.warning("Reference missing")

    # =============================
    # OUTPUT
    # =============================
    st.subheader("🔍 Matched SNPs")
    st.dataframe(matched_df)

    st.subheader("🏆 Top SNPs per Sub-Trait")

    if len(matched_df) > 0:
        df_top = matched_df.copy()
        df_top["p_value"] = pd.to_numeric(df_top["p_value"], errors="coerce")
        df_top = df_top.dropna(subset=["Sub_Trait", "p_value"])
        df_top = df_top.sort_values(by="p_value")

        top_snps = df_top.groupby("Sub_Trait", group_keys=False).apply(get_top_snps)
        st.dataframe(top_snps)

    st.download_button(
        "⬇ Download Results",
        matched_df.to_csv(index=False),
        "matched_snps.csv"
    )

# =============================
# SEARCH
# =============================
st.subheader("🔍 Search SNP by rsID")

search = st.text_input("Enter rsID")

if search:
    result = db_df[
        db_df["snp_id"].str.contains(search.strip(), case=False, na=False)
    ]

    if len(result) > 0:
        st.success(f"Found {len(result)} SNP(s)")
        st.dataframe(result)
    else:
        st.warning("No SNP found")
