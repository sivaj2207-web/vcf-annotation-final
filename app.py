import streamlit as st
import pandas as pd
import sqlite3
import gzip

# =============================
# TITLE
# =============================
st.title("🧬 VCF Variant Annotation Tool")

# =============================
# LOAD DATABASE (SAS + MATCHING)
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
# PARSE VCF (SAFE FOR LARGE FILES)
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

            if count >= limit:  # 🔥 prevents crash
                break

        except:
            continue

    return pd.DataFrame(data, columns=[
        "snp_id", "chromosome", "position",
        "ref", "alt", "genotype"
    ])

# =============================
# USER VCF INPUT
# =============================
st.subheader("📂 Upload User VCF")

user_vcf_file = st.file_uploader("Upload User VCF", type=["vcf"])

vcf_df = None
if user_vcf_file:
    vcf_df = parse_vcf(user_vcf_file, limit=10000)

# =============================
# IBS REFERENCE INPUT (NEW)
# =============================
st.subheader("🧬 IBS Reference Input")

ibs_ref_file = st.file_uploader(
    "Upload Reference VCF for IBS (.vcf or .vcf.gz)",
    type=["vcf", "gz"]
)

def load_reference(file):
    try:
        if file is None:
            return pd.DataFrame()

        if file.name.endswith(".gz"):
            return parse_vcf(gzip.open(file, "rt"), limit=10000)
        else:
            return parse_vcf(file, limit=10000)

    except Exception as e:
        st.warning(f"Reference load failed: {e}")
        return pd.DataFrame()

# =============================
# MATCHING (DB ONLY)
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
# IBS CALCULATION
# =============================
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
# RUN PIPELINE
# =============================
if vcf_df is not None:

    st.write("🔄 Processing...")

    # MATCHING + SAS (FULL DATA)
    matched_df = annotate(vcf_df)

    st.success("✅ Done")

    st.write("📊 Total SNPs:", len(vcf_df))
    st.write("✅ Matched SNPs:", len(matched_df))

    st.metric("🧬 SAS Score (%)", f"{sas_score(matched_df):.2f}")

    # =============================
    # IBS (SEPARATE REFERENCE)
    # =============================
    st.subheader("🧬 IBS Similarity")

    ref_df = load_reference(ibs_ref_file)

    if ibs_ref_file is None:
        st.warning("Upload IBS reference file")
    elif not ref_df.empty:

        # 🔥 find overlap first
        common_ids = set(vcf_df["snp_id"]) & set(ref_df["snp_id"])

        vcf_common = vcf_df[vcf_df["snp_id"].isin(common_ids)]
        ref_common = ref_df[ref_df["snp_id"].isin(common_ids)]

        # 🔥 random sampling
        if len(vcf_common) > 3000:
            vcf_common = vcf_common.sample(3000, random_state=42)
        if len(ref_common) > 3000:
            ref_common = ref_common.sample(3000, random_state=42)

        comp = pd.merge(
            vcf_common,
            ref_common,
            on="snp_id",
            suffixes=("_user", "_ref")
        )

        st.write("🔗 Common SNPs:", len(comp))

        if len(comp) > 0:
            ibs = calculate_ibs(comp)
            st.metric("🧬 IBS (%)", f"{ibs:.2f}")
        else:
            st.warning("No overlapping SNPs")

    else:
        st.warning("Reference loading failed")

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
        db_df["snp_id"].str.contains(search.strip().lower(), case=False, na=False)
    ]

    if len(result) > 0:
        st.success(f"Found {len(result)} SNP(s)")
        st.dataframe(result)
    else:
        st.warning("No SNP found")
