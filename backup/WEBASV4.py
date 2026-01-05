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
# COMMON NAME DICTIONARY
# ===============================
COMMON_NAMES = {
    "Homo sapiens":        {"English": "Human",     "ภาษาไทย": "มนุษย์"},
    "Pan troglodytes":     {"English": "Chimpanzee","ภาษาไทย": "ชิมแพนซี"},
    "Gorilla gorilla":     {"English": "Gorilla",   "ภาษาไทย": "กอริลลา"},
    "Pan paniscus":        {"English": "Bonobo",    "ภาษาไทย": "โบนาโบ"},
    "Elephas maximus":     {"English": "Asian elephant", "ภาษาไทย": "ช้างเอเชีย"},
    "Loxodonta africana":  {"English": "African elephant", "ภาษาไทย": "ช้างแอฟริกา"},
    "Trichechus manatus":  {"English": "Manatee",   "ภาษาไทย": "มานาตี"},
    "Dugong dugon":        {"English": "Dugong",    "ภาษาไทย": "ดูงอง"},
    "Procavia capensis":   {"English": "Rock hyrax","ภาษาไทย": "ไฮแร็กซ์หิน"},
    
    # …คุณสามารถเพิ่มรายชื่ออื่นได้เอง
}

# ===============================
# LANGUAGE SUPPORT
# ===============================
LANG = {
    "English": {
        "title": "🧬 ND5 Phylogenetic Explorer",
        "fetch_list": "📥 Fetch Species (20)",
        "species_select": "Select species (min 2)",
        "fetch_nd5": "🔍 Fetch ND5",
        "identity": "📊 Identity Matrix",
        "build_tree": "🌳 Build Phylogenetic Tree",
        "theme": "Theme",
        "light": "Light",
        "dark": "Dark",
        "no_data": "No ND5 data yet",
        "fetching": "Fetching species info..."
    },
    "ภาษาไทย": {
        "title": "🧬 ตัวสำรวจสายวิวัฒนาการ ND5",
        "fetch_list": "📥 ดึงสายพันธุ์ (20 ตัว)",
        "species_select": "เลือกสายพันธุ์ (อย่างน้อย 2)",
        "fetch_nd5": "🔍 ดึง ND5",
        "identity": "📊 ตาราง % ความเหมือน",
        "build_tree": "🌳 สร้างแผนภาพต้นไม้",
        "theme": "ธีม",
        "light": "สว่าง",
        "dark": "มืด",
        "no_data": "ยังไม่มีข้อมูล ND5",
        "fetching": "กำลังดึงข้อมูลสายพันธุ์..."
    }
}

# ===============================
# SIDEBAR CONTROLS
# ===============================
language = st.sidebar.selectbox("🌐 Language / ภาษา", ["English", "ภาษาไทย"])
T = LANG[language]

theme = st.sidebar.selectbox(T["theme"], [T["light"], T["dark"]], index=0)
if theme == T["dark"]:
    st.markdown("""
    <style>
    .stApp { background: #0e1117; color: white; }
    div[data-testid="stDataFrame"] { background: #1e1e1e; }
    </style>
    """, unsafe_allow_html=True)

st.title(T["title"])

# ===============================
# NCBI CONFIG
# ===============================
Entrez.email = "poonthakorn@gmail.com"

# ===============================
# FETCH 20 SPECIES
# ===============================
species_map = st.session_state.get("species_map", {})

if st.button(T["fetch_list"]):
    query = '"mitochondrion"[Filter] AND refseq[Filter] AND complete genome'
    search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=20)
    search_record = Entrez.read(search_handle)
    id_list = search_record["IdList"]

    new_map = {}
    with st.spinner(T["fetching"]):
        for uid in id_list:
            with Entrez.efetch(db="nucleotide", id=uid, rettype="gb", retmode="text") as fetch_handle:
                rec = SeqIO.read(fetch_handle, "genbank")
                sci = rec.annotations["organism"]  # scientific
                acc = rec.annotations["accessions"][0]
                new_map[f"{sci} ({acc})"] = acc

    st.session_state["species_map"] = new_map
    species_map = new_map

# ===============================
# SHOW SEARCHABLE SPECIES LIST
# ===============================
selected = []
if species_map:
    display_list = []
    display_to_key = {}
    for sci_acc in species_map.keys():
        sci_name = sci_acc.split(" (")[0]
        if sci_name in COMMON_NAMES:
            display_label = COMMON_NAMES[sci_name][language]
        else:
            display_label = sci_name  # fallback

        if display_label not in display_to_key:
            display_to_key[display_label] = sci_acc
            display_list.append(display_label)

    selected_display = st.multiselect(
        T["species_select"],
        options=display_list
    )

    selected = [display_to_key[d] for d in selected_display]

# ===============================
# ND5 EXTRACTION
# ===============================
def extract_nd5(accession_id):
    with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text") as handle:
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
    nd5_seqs = {}
    for key in selected:
        acc = species_map[key]
        seq = extract_nd5(acc)
        if seq:
            nd5_seqs[key] = seq
    st.session_state["nd5_seqs"] = nd5_seqs

# ===============================
# SHOW ND5 LENGTH INFO
# ===============================
if "nd5_seqs" in st.session_state:
    for k, seq in st.session_state["nd5_seqs"].items():
        sci_name = k.split(" (")[0]
        display = COMMON_NAMES.get(sci_name, {}).get(language, sci_name)
        st.text(f"{display} → {len(seq)} bp")

# ===============================
# IDENTITY TABLE
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
                    row.append(100)
                else:
                    row.append(round(percent_identity(seqs[a], seqs[b]), 2))
            matrix.append(row)
        df = pd.DataFrame(matrix, index=[
            COMMON_NAMES.get(n.split(" (")[0], {}).get(language, n.split(" (")[0])
            for n in names
        ], columns=[
            COMMON_NAMES.get(n.split(" (")[0], {}).get(language, n.split(" (")[0])
            for n in names
        ])
        st.dataframe(df)

# ===============================
# BUILD TREE
# ===============================
if st.button(T["build_tree"]):
    seqs = st.session_state.get("nd5_seqs", {})
    if not seqs:
        st.warning(T["no_data"])
    else:
        fasta_io = StringIO()
        for k, seq in seqs.items():
            fasta_io.write(f">{k}\n{seq}\n")
        fasta_io.seek(0)

        alignment = AlignIO.read(fasta_io, "fasta")
        dm = DistanceCalculator("identity").get_distance(alignment)
        tree = DistanceTreeConstructor().nj(dm)
        tree.ladderize()

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        Phylo.draw(tree, axes=ax, do_show=False)
        st.pyplot(fig)
