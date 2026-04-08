import streamlit as st
import pandas as pd
import sqlite3

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

    # Clean column names
    df.columns = df.columns.astype(str)
    df.columns = df.columns.str.strip()
    df.columns = df.columns.str.replace(r"\s+", "_", regex=True)

    # Rename safely
    if "SNP_ID" in df.columns:
        df = df.rename(columns={"SNP_ID": "snp_id"})
    if "P-value" in df.columns:
        df = df.rename(columns={"P-value": "p_value"})
    if "Phenotype_Description" in df.columns:
        df = df.rename(columns={"Phenotype_Description": "phenotype"})

    # Safety check
    required = ["snp_id", "Sub_Trait", "p_value"]
    for col in required:
        if col not in df.columns:
            st.error(f"❌ Missing column: {col}")
            st.write(df.columns)
            st.stop()

    df["snp_id"] = df["snp_id"].astype(str).str.strip()

    return df


db_df = load_database()

# =============================
# PREVIEW
# =============================
st.subheader("📊 Database Preview")
st.dataframe(db_df.head())

# =============================
# VCF UPLOAD
# =============================
uploaded_vcf = st.file_uploader("Upload VCF file", type=["vcf"])

# =============================
# PARSE VCF
# =============================
def parse_vcf(file):
    data = []

    for line in file:
        line = line.decode("utf-8")

        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")

        if len(cols) < 10:
            continue

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

    return pd.DataFrame(data, columns=[
        "snp_id", "chromosome", "position",
        "ref", "alt", "genotype"
    ])

# =============================
# MATCHING (SAFE SQL)
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
        d.phenotype
    FROM vcf v
    JOIN db d
    ON v.snp_id = d.snp_id
    """

    result = pd.read_sql_query(query, conn)
    conn.close()

    return result

# =============================
# SAS SCORE (SAFE)
# =============================
def sas_score(df):
    if "SAS_MAF" in df.columns and len(df) > 0:
        return df["SAS_MAF"].mean() * 100
    return 0

# =============================
# TOP SNP LOGIC
# =============================
def get_top_snps(group):
    if len(group) >= 10:
        return group.head(10)
    else:
        return group.head(5)

# =============================
# RUN PIPELINE
# =============================
if uploaded_vcf is not None:

    st.write("🔄 Processing VCF...")

    try:
        vcf_df = parse_vcf(uploaded_vcf)
        matched_df = annotate(vcf_df)
    except Exception as e:
        st.error(f"❌ Error: {e}")
        st.stop()

    st.success("✅ Annotation Complete")

    st.write("📊 Total SNPs:", len(vcf_df))
    st.write("✅ Matched SNPs:", len(matched_df))

    # =============================
    # SAS SCORE
    # =============================
    sas = sas_score(matched_df)
    st.metric("🧬 SAS Score (%)", f"{sas:.2f}")

    # =============================
    # SHOW MATCHED
    # =============================
    st.subheader("🔍 Matched SNPs")
    st.dataframe(matched_df)

    # =============================
    # TOP SNPs
    # =============================
    st.subheader("🏆 Top SNPs per Sub-Trait")

    if len(matched_df) > 0:
        df_top = matched_df.copy()

        df_top["p_value"] = pd.to_numeric(df_top["p_value"], errors="coerce")
        df_top = df_top.dropna(subset=["Sub_Trait", "p_value"])
        df_top = df_top.sort_values(by="p_value")

        top_snps = df_top.groupby("Sub_Trait", group_keys=False).apply(get_top_snps)

        st.dataframe(top_snps)

    else:
        st.warning("No matched SNPs")

    # =============================
    # DOWNLOAD
    # =============================
    csv = matched_df.to_csv(index=False).encode("utf-8")

    st.download_button(
        "⬇ Download Results",
        csv,
        "matched_snps.csv",
        "text/csv"
    )

# =============================
# RSID SEARCH
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
