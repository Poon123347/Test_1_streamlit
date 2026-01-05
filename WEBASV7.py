import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Align import PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from io import StringIO

# ===============================
# PAGE CONFIG
# ===============================
st.set_page_config(
    page_title="ND5 Phylogenetic Explorer",
    layout="centered"
)

# ===============================
# NCBI CONFIG
# ===============================
Entrez.email = "poonthakorn@gmail.com"

# ===============================
# SPECIES ACCESSIONS (PREDEFINED)
# ===============================
species_accessions = {
    "Human": "NC_012920.1",
    "Chimpanzee": "NC_001643.1",
    "Bonobo": "NC_001644.1",
    "Western Gorilla": "NC_001645.1",
    "Eastern Gorilla": "NC_011120.1",
    "Sumatran Orangutan": "NC_002083.1",
    "Bornean Orangutan": "NC_001646.1",
    "Common Gibbon": "NC_002082.1",
    "Rhesus Macaque": "NC_005943.1",
    "Crab-eating Macaque": "NC_012670.1",
    "Barbary Ape": "NC_002764.1",
    "HamadrayasBaboons": "NC_001992.1",
    "GreenMonkey": "NC_008066.1",
    "TantalusMonkey": "NC_009748.1",
    "DuskyLeafMonkey": "NC_006900.1",
    "ProboscisMonkey": "NC_008216.1",
    "RedShankedDouc": "NC_008220.1",
    "TonkinSnubNoseMonkey": "NC_015485.1",
    "ColobusMonkey": "NC_006901.1",
    "BlackColobus": "NC_008219.1",
    "Babirusa": "NC_012732.1",
    "Pig": "NC_000845.1",
    "Cow": "NC_001567.1",
    "Horse": "NC_001640.1",
    "Dog": "NC_002008.4",
    "Cat": "NC_001700.1",
    "Sheep": "NC_001941.1",
    "Goat": "NC_005044.2",
    "Camel": "NC_009628.1",
    "Pig-tailedMacaque": "NC_001715.1",
    "NorwayRat": "NC_001665.2",
    "Mouse": "NC_005089.1",
    "Rabbit": "NC_001913.1",
    "Elephant": "NC_003319.1",
    "WoollyMammoth": "NC_007596.1",
    "BrownBear": "NC_003427.1",
    "PolarBear": "NC_003428.1",
    "Squirrel": "NC_004742.1",
    "Beaver": "NC_011112.1",
    "Koala": "NC_026731.1",
    "Platypus": "NC_000891.1",
    "Chicken": "NC_001323.1",
    "Turkey": "NC_010195.1",
    "Duck": "NC_009684.1",
    "Ostrich": "NC_002785.1",
    "ZebraFinch": "NC_007897.1",
    "Lizard": "NC_001323.1",  # same reference but can replace with reptile of interest
    "Salmon": "NC_002980.1",
    "Zebrafish": "NC_002333.2"
}

# ===============================
# LANGUAGE SUPPORT
# ===============================
LANG = {
    "English": {
        "title": "🧬 ND5 Phylogenetic Explorer",
        "species_select": "Select species (min 2)",
        "fetch_nd5": "🔍 Fetch ND5",
        "identity": "📊 Identity Matrix",
        "build_tree": "🌳 Build Phylogenetic Tree",
        "no_data": "No ND5 data yet"
    },
    "ภาษาไทย": {
        "title": "🧬 ตัวสำรวจสายวิวัฒนาการ ND5",
        "species_select": "เลือกสายพันธุ์ (อย่างน้อย 2)",
        "fetch_nd5": "🔍 ดึง ND5",
        "identity": "📊 ตาราง % ความเหมือน",
        "build_tree": "🌳 สร้างแผนภาพต้นไม้",
        "no_data": "ยังไม่มีข้อมูล ND5"
    }
}

# ===============================
# SIDEBAR
# ===============================
language = st.sidebar.selectbox("🌐 Language / ภาษา", ["English", "ภาษาไทย"])
T = LANG[language]

st.title(T["title"])

# ===============================
# SPECIES SELECTOR
# ===============================
selected_species = st.multiselect(
    T["species_select"],
    options=list(species_accessions.keys())
)

selected_accessions = {
    name: species_accessions[name]
    for name in selected_species
}

# ===============================
# ND5 EXTRACTION
# ===============================
def extract_nd5(accession_id):
    with Entrez.efetch(
        db="nucleotide",
        id=accession_id,
        rettype="gb",
        retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "genbank")

    for feat in record.features:
        if feat.type == "CDS" and "gene" in feat.qualifiers:
            if feat.qualifiers["gene"][0].lower() == "nd5":
                return feat.extract(record.seq)

    return None

# ===============================
# PERCENT IDENTITY
# ===============================
def percent_identity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    score = aligner.score(seq1, seq2)
    return score / max(len(seq1), len(seq2)) * 100

# ===============================
# FETCH ND5
# ===============================
if st.button(T["fetch_nd5"]):
    if len(selected_accessions) < 2:
        st.warning(T["species_select"])
    else:
        nd5_seqs = {}
        with st.spinner("Fetching ND5..."):
            for name, acc in selected_accessions.items():
                seq = extract_nd5(acc)
                if seq:
                    nd5_seqs[name] = seq

        st.session_state["nd5_seqs"] = nd5_seqs

# ===============================
# SHOW ND5 LENGTH
# ===============================
if "nd5_seqs" in st.session_state:
    st.subheader("ND5 length")
    for name, seq in st.session_state["nd5_seqs"].items():
        st.text(f"{name} → {len(seq)} bp")

# ===============================
# IDENTITY MATRIX
# ===============================
if st.button(T["identity"]):
    seqs = st.session_state.get("nd5_seqs", {})
    if not seqs:
        st.warning(T["no_data"])
    else:
        names = list(seqs.keys())
        matrix = []

        for i, a in enumerate(names):
            row = []
            for j, b in enumerate(names):
                if i == j:
                    row.append(100.0)
                else:
                    row.append(round(percent_identity(seqs[a], seqs[b]), 2))
            matrix.append(row)

        df = pd.DataFrame(matrix, index=names, columns=names)
        st.dataframe(df)

# ===============================
# BUILD PHYLOGENETIC TREE
# ===============================
if st.button(T["build_tree"]):
    seqs = st.session_state.get("nd5_seqs", {})
    if not seqs:
        st.warning(T["no_data"])
    else:
        fasta_io = StringIO()
        for name, seq in seqs.items():
            fasta_io.write(f">{name}\n{seq}\n")
        fasta_io.seek(0)

        alignment = AlignIO.read(fasta_io, "fasta")
        dm = DistanceCalculator("identity").get_distance(alignment)
        tree = DistanceTreeConstructor().nj(dm)
        tree.ladderize()

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        Phylo.draw(tree, axes=ax, do_show=False)
        st.pyplot(fig)
