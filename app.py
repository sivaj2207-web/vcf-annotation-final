import streamlit as st
import pandas as pd
import sqlite3
import requests
import io

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

    df["snp_id"] = df["snp_id"].astype(str).str.strip().str.lower()
    return df

db_df = load_database()

# =============================
# PREVIEW
# =============================
st.subheader("📊 Database Preview")
st.dataframe(db_df.head())

# =============================
# PARSE VCF (SAFE)
# =============================
def parse_vcf(file, limit=10000):
    data = []
    count = 0

    for line in file:
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")
        if len(cols) < 10:
            continue

        try:
            snp_id = cols[2].strip().lower()
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
            count += 1

            if count >= limit:  # prevents crash
                break

        except:
            continue

    return pd.DataFrame(data, columns=[
        "snp_id", "chromosome", "position",
        "ref", "alt", "genotype"
    ])

# =============================
# INPUT OPTIONS
# =============================
st.subheader("📂 Input VCF")

input_method = st.radio(
    "Choose input method:",
    ["Upload VCF", "Use file path", "Paste URL"]
)

vcf_df = None

# -------- Upload --------
if input_method == "Upload VCF":
    uploaded_vcf = st.file_uploader("Upload VCF file", type=["vcf"])
    if uploaded_vcf:
        st.write(f"File size: {uploaded_vcf.size / (1024*1024):.2f} MB")
        if st.button("🚀 Process File"):
            with st.spinner("⏳ Processing VCF..."):
                vcf_df = parse_vcf(uploaded_vcf, limit=10000)

# -------- File Path --------
elif input_method == "Use file path":
    file_path = st.text_input("Enter file path")

    if file_path:
        if st.button("🚀 Process File"):
            try:
                with st.spinner("⏳ Reading file..."):
                    with open(file_path, "rb") as f:
                        vcf_df = parse_vcf(f, limit=10000)
            except Exception as e:
                st.error(f"❌ {e}")

# -------- URL --------
elif input_method == "Paste URL":
    url = st.text_input("Paste VCF file URL")

    if url:
        if st.button("🚀 Process File"):
            try:
                with st.spinner("⏳ Downloading + processing VCF..."):
                    response = requests.get(url, stream=True)

                    if response.status_code == 200:
                        file_like = io.BytesIO(response.content)
                        vcf_df = parse_vcf(file_like, limit=10000)
                    else:
                        st.error("❌ Failed to download file")

            except Exception as e:
                st.error(f"❌ {e}")

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
# SAS SCORE
# =============================
def sas_score(df):
    if "SAS_MAF" in df.columns and len(df) > 0:
        df["SAS_MAF"] = pd.to_numeric(df["SAS_MAF"], errors="coerce")
        return df["SAS_MAF"].mean() * 100
    return 0

# =============================
# TOP SNP
# =============================
def get_top_snps(group):
    return group.head(10) if len(group) >= 10 else group.head(5)

# =============================
# RUN PIPELINE
# =============================
if vcf_df is not None:

    st.success("✅ Processing Complete")

    matched_df = annotate(vcf_df)

    st.write("📊 Total SNPs:", len(vcf_df))
    st.write("✅ Matched SNPs:", len(matched_df))

    st.metric("🧬 SAS Score (%)", f"{sas_score(matched_df):.2f}")

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
    else:
        st.warning("No matched SNPs")

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
        db_df["snp_id"].str.contains(search.strip().lower(), case=False, na=False)
    ]

    if len(result) > 0:
        st.success(f"Found {len(result)} SNP(s)")
        st.dataframe(result)
    else:
        st.warning("No SNP found")
