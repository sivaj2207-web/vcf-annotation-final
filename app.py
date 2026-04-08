import streamlit as st
import pandas as pd

# =============================
# TITLE
# =============================
st.title("🧬 VCF Variant Annotation Tool")

# =============================
# LOAD DATABASE (CSV)
# =============================
@st.cache_data
def load_database():
    df = pd.read_csv("annotated_table.csv")

    df.columns = df.columns.str.strip()

    df = df.rename(columns={
        "SNP_ID": "snp_id",
        "Chromosome No.": "chromosome",
        "Chromosome Position (GRCh38)": "position",
        "Reference_Allele/MAJOR_ALLELE": "ref",
        "Alternate_Allele": "alt",
        "Phenotypic Description": "phenotype",
        "Cytogenic region": "cytogenic_region",
        "P-value": "p_value",
        "Beta/Odds Ratio (Effect size)": "beta_or",
        "Sample size": "sample_size"
    })

    return df


# =============================
# FILE UPLOAD
# =============================
uploaded_vcf = st.file_uploader("Upload VCF file", type=["vcf"])


# =============================
# PARSE VCF
# =============================
def parse_vcf(file):
    vcf_data = []

    for line in file:
        line = line.decode("utf-8")

        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")

        chr = cols[0]
        pos = int(cols[1])
        snp_id = cols[2]
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

        vcf_data.append([snp_id, chr, pos, ref, alt, genotype])

    return pd.DataFrame(vcf_data, columns=[
        "snp_id", "chromosome", "position", "ref", "alt", "genotype"
    ])


# =============================
# ANNOTATION (NO SQLITE)
# =============================
def annotate(vcf_df):
    db_df = load_database()

    result = pd.merge(
        vcf_df,
        db_df,
        on="snp_id",
        how="left"
    )

    return result


# =============================
# SAS SCORING
# =============================
def sas_score(df):
    matched = df[df["SAS_MAF"].notnull()]

    if len(matched) == 0:
        return 0

    return matched["SAS_MAF"].mean() * 100


# =============================
# RUN PIPELINE
# =============================
if uploaded_vcf is not None:

    st.write("🔄 Processing VCF...")

    vcf_df = parse_vcf(uploaded_vcf)
    result_df = annotate(vcf_df)

    score = sas_score(result_df)

    st.success("✅ Annotation Complete")

    st.metric("🧬 SAS Score (%)", f"{score:.2f}")

    st.dataframe(result_df.head(100))

    csv = result_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇ Download Full Results",
        csv,
        "annotated_output.csv",
        "text/csv"
    )


# =============================
# DATABASE VIEWER
# =============================
st.subheader("📊 Database Viewer")

db_df = load_database()

if st.checkbox("Show full database"):
    st.dataframe(db_df)

st.write("Total SNPs:", len(db_df))
