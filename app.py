import streamlit as st
import pandas as pd
import sqlite3
import gzip

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

    required = ["snp_id", "Sub_Trait", "p_value"]
    for col in required:
        if col not in df.columns:
            st.error(f"Missing column: {col}")
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
# INPUT VCF (USER)
# =============================
st.subheader("📂 Input VCF")

input_method = st.radio("Choose input method:", ["Upload VCF", "Use file path"])
vcf_df = None

if input_method == "Upload VCF":
    uploaded_vcf = st.file_uploader("Upload VCF", type=["vcf"])
    if uploaded_vcf:
        vcf_df = parse_vcf(uploaded_vcf)

else:
    path = st.text_input("Enter file path")
    if path:
        try:
            with open(path, "rb") as f:
                vcf_df = parse_vcf(f)
        except Exception as e:
            st.error(e)

# =============================
# REFERENCE INPUT (IBS)
# =============================
st.subheader("📂 Reference VCF (for IBS)")

ref_option = st.radio("Reference input:", ["Upload", "Path"])
reference_df = pd.DataFrame()

if ref_option == "Upload":
    ref_file = st.file_uploader("Upload reference VCF", type=["vcf", "gz"])

    if ref_file:
        try:
            if ref_file.name.endswith(".gz"):
                with gzip.open(ref_file, "rt") as f:
                    reference_df = parse_vcf(f)
            else:
                reference_df = parse_vcf(ref_file)
        except Exception as e:
            st.error(e)

else:
    ref_path = st.text_input("Reference file path")
    if ref_path:
        try:
            if ref_path.endswith(".gz"):
                with gzip.open(ref_path, "rt") as f:
                    reference_df = parse_vcf(f)
            else:
                with open(ref_path, "rb") as f:
                    reference_df = parse_vcf(f)
        except Exception as e:
            st.error(e)

# =============================
# MATCHING
# =============================
def annotate(vcf_df):
    conn = sqlite3.connect(":memory:")
    vcf_df.to_sql("vcf", conn, index=False, if_exists="replace")
    db_df.to_sql("db", conn, index=False, if_exists="replace")

    query = """
    SELECT 
        v.snp_id, v.chromosome, v.position, v.ref, v.alt, v.genotype,
        d.Gene, d.Sub_Trait, d.p_value, d.phenotype, d.SAS_MAF
    FROM vcf v
    JOIN db d ON v.snp_id = d.snp_id
    """

    result = pd.read_sql_query(query, conn)
    conn.close()
    return result

# =============================
# SAS SCORE
# =============================
def sas_score(df):
    if "SAS_MAF" not in df.columns or len(df) == 0:
        return 0
    return df["SAS_MAF"].mean() * 100

# =============================
# IBS
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
# RUN PIPELINE
# =============================
if vcf_df is not None:

    st.write("Processing...")

    matched_df = annotate(vcf_df)

    st.success("Done")

    st.write("Total SNPs:", len(vcf_df))
    st.write("Matched SNPs:", len(matched_df))

    # SAS
    st.metric("🧬 SAS Score (%)", f"{sas_score(matched_df):.2f}")

    # IBS
    st.subheader("🧬 IBS Similarity")

    if not reference_df.empty:
        comp = pd.merge(vcf_df, reference_df, on="snp_id", suffixes=("_user", "_ref"))

        st.write("Common SNPs:", len(comp))

        if len(comp) > 0:
            st.metric("IBS (%)", f"{calculate_ibs(comp):.2f}")
        else:
            st.warning("No overlap")
    else:
        st.warning("Provide reference VCF")

    # Matched SNPs
    st.subheader("🔍 Matched SNPs")
    st.dataframe(matched_df)

    # Top SNPs
    st.subheader("🏆 Top SNPs")

    df_top = matched_df.copy()
    df_top["p_value"] = pd.to_numeric(df_top["p_value"], errors="coerce")
    df_top = df_top.dropna(subset=["Sub_Trait", "p_value"])
    df_top = df_top.sort_values("p_value")

    top = df_top.groupby("Sub_Trait", group_keys=False).apply(get_top_snps)
    st.dataframe(top)

    # Download
    st.download_button(
        "Download CSV",
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
    st.dataframe(res if len(res) else pd.DataFrame({"Result": ["No match"]}))
