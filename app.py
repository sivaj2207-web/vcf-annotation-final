import streamlit as st
import pandas as pd
import sqlite3
import gzip
import os

st.title("🧬 VCF Variant Annotation Tool")

# =============================
# LOAD DATABASE (SAFE)
# =============================
@st.cache_data
def load_database():
    try:
        df = pd.read_csv("annotated_table_updated.csv")
    except:
        st.error("❌ Database file missing")
        st.stop()

    df.columns = df.columns.str.strip().str.replace(" ", "_")

    if "SNP_ID" in df.columns:
        df.rename(columns={"SNP_ID": "snp_id"}, inplace=True)
    if "P-value" in df.columns:
        df.rename(columns={"P-value": "p_value"}, inplace=True)

    return df

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
            snp_id = cols[2]
            genotype = cols[9].split(":")[0]

            if genotype == "0/0":
                genotype = "0|0"
            elif genotype in ["0/1", "1/0"]:
                genotype = "0|1"
            elif genotype == "1/1":
                genotype = "1|1"

            data.append([snp_id, genotype])
        except:
            continue

    return pd.DataFrame(data, columns=["snp_id", "genotype"])

# =============================
# IBS CALCULATION
# =============================
def calculate_ibs(df):
    score = 0

    for _, row in df.iterrows():
        g1 = str(row["genotype_user"])
        g2 = str(row["genotype_ref"])

        if g1 == g2:
            score += 2
        elif "0|1" in g1 or "0|1" in g2:
            score += 1

    return (score / (2 * len(df))) * 100 if len(df) else 0

# =============================
# UI INPUT
# =============================
uploaded_vcf = st.file_uploader("📂 Upload your VCF", type=["vcf"])

# =============================
# RUN ONLY AFTER UPLOAD
# =============================
if uploaded_vcf:

    st.write("Processing...")

    db_df = load_database()

    # parse user VCF
    vcf_df = parse_vcf(uploaded_vcf)

    # =============================
    # ANNOTATION
    # =============================
    conn = sqlite3.connect(":memory:")
    vcf_df.to_sql("vcf", conn, index=False, if_exists="replace")
    db_df.to_sql("db", conn, index=False, if_exists="replace")

    query = """
    SELECT v.*, d.*
    FROM vcf v
    JOIN db d ON v.snp_id = d.snp_id
    """

    matched_df = pd.read_sql_query(query, conn)

    st.success("✅ Done")

    st.write("Matched SNPs:", len(matched_df))
    st.dataframe(matched_df.head())

    # =============================
    # IBS (SAFE LOAD)
    # =============================
    st.subheader("🧬 IBS Similarity")

    if os.path.exists("reference.vcf.gz"):

        try:
            with gzip.open("reference.vcf.gz", "rt") as f:
                ref_df = parse_vcf(f)

            comp = pd.merge(
                vcf_df,
                ref_df,
                on="snp_id",
                suffixes=("_user", "_ref")
            )

            if len(comp) > 0:
                ibs = calculate_ibs(comp)
                st.metric("IBS (%)", f"{ibs:.2f}")
            else:
                st.warning("No overlapping SNPs")

        except Exception as e:
            st.error(f"IBS error: {e}")

    else:
        st.warning("⚠️ reference.vcf.gz not found")
