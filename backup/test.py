import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO, pairwise2
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np

# ===============================
# Page config
# ===============================
st.set_page_config(page_title="DNA & ND5 Explorer", layout="wide")

# ===============================
# Language dictionary
# ===============================
LANG = {
    "ภาษาไทย": {
        "title": "📊 ตัวสำรวจ DNA และยีน ND5",
        "load_file": "📥 แหล่งข้อมูล FASTA",
        "paste_label": "วางข้อความ FASTA ที่นี่ (หรืออัปโหลดด้านล่าง)",
        "upload_label": "อัปโหลดไฟล์ FASTA (ตัวเลือก)",
        "fetch_list": "📥 ดึงตัวอย่างสายพันธุ์จาก NCBI (20)",
        "species_select": "เลือกสายพันธุ์",
        "fetch_nd5": "🔍 ดึงลำดับ ND5",
        "identity": "📊 ตารางความเหมือน (%)",
        "theme": "ธีม",
        "light": "สว่าง",
        "dark": "มืด",
        "no_data": "ยังไม่มีข้อมูล ND5",
        "fetching_species": "กำลังดึงรายชื่อสายพันธุ์...",
        "fetching_nd5": "กำลังดึงลำดับ ND5...",
        "select_prompt": "กรุณาเลือกสายพันธุ์จากด้านซ้าย",
        "need_two": "ต้องมีอย่างน้อย 2 ลำดับ",
        "failed_fetch": "ไม่สามารถดึงข้อมูลได้",
        "fetched": "ดึง ND5 ได้ {n}",
        "fasta_loaded": "โหลด FASTA แล้ว: {n} รายการ",
        "no_local_fasta": "ยังไม่มีไฟล์ FASTA ในอินพุต — วางข้อความหรืออัปโหลดได้",
        "fetched_species": "ได้รายการ {n} สายพันธุ์จาก NCBI",
        "download_csv": "ดาวน์โหลด CSV"
    },
    "English": {
        "title": "📊 DNA & ND5 Explorer",
        "load_file": "📥 FASTA input",
        "paste_label": "Paste FASTA text here (or upload below)",
        "upload_label": "Upload FASTA file (optional)",
        "fetch_list": "📥 Fetch example species from NCBI (20)",
        "species_select": "Select species",
        "fetch_nd5": "🔍 Fetch ND5",
        "identity": "📊 Identity matrix (%)",
        "theme": "Theme",
        "light": "Light",
        "dark": "Dark",
        "no_data": "No ND5 data",
        "fetching_species": "Fetching species list...",
        "fetching_nd5": "Fetching ND5 sequences...",
        "select_prompt": "Please select species from sidebar",
        "need_two": "At least 2 sequences required",
        "failed_fetch": "Failed to fetch data",
        "fetched": "Fetched ND5: {n}",
        "fasta_loaded": "FASTA loaded: {n} records",
        "no_local_fasta": "No FASTA input — paste or upload",
        "fetched_species": "Fetched {n} species from NCBI",
        "download_csv": "Download CSV"
    }
}

# ===============================
# Sidebar: language + theme (UI)
# ===============================
with st.sidebar:
    language = st.selectbox("🌐 Language / ภาษา", ["ภาษาไทย", "English"], index=0)
    # theme selectbox uses localized labels
    theme = st.selectbox(LANG[language]["theme"], [LANG[language]["light"], LANG[language]["dark"]], index=0)

T = LANG[language]

# ===============================
# Improved CSS (covers selectbox/multiselect/listbox/alerts/buttons/cursor)
# ===============================
THEME_CSS = """
<style>
/* Smooth transitions */
html, body, .stApp, .block-container, section, div {
  transition: background-color 0.2s ease, color 0.2s ease, border-color 0.2s ease;
}

/* Base app colors */
.stApp, .reportview-container, .block-container, main {
  background-color: VAR_BG !important;
  color: VAR_TEXT !important;
}

/* Sidebar */
section[data-testid="stSidebar"] {
  background-color: VAR_SIDEBAR_BG !important;
  color: VAR_TEXT !important;
  border-right: 1px solid VAR_BORDER !important;
  padding-top: 0.6rem !important;
}

/* Panels, cards, expanders */
div[data-testid="stVerticalBlock"], div[data-testid="stMarkdownContainer"],
div[data-testid="stCodeBlock"], .stDataFrame, .stExpander, section > div[role="region"] {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
  border-radius: 10px !important;
  padding: 0.5rem 0.8rem !important;
}

/* Code blocks */
pre, code, div[data-testid="stCodeBlock"] pre, div[data-testid="stCodeBlock"] code {
  background-color: VAR_CODE_BG !important;
  color: VAR_CODE_TEXT !important;
  border: 1px solid VAR_BORDER !important;
  border-radius: 8px !important;
}

/* Inputs: text, textarea */
input, textarea, .stTextInput, .stTextArea {
  background-color: VAR_INPUT_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
  border-radius: 8px !important;
}

/* Selectbox / Multiselect / Radio / Checkbox containers */
div[data-testid="stSelectbox"], div[data-testid="stMultiSelect"], div[data-testid="stRadio"], div[data-testid="stCheckbox"] {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
  border-radius: 8px !important;
  padding: 0.15rem 0.3rem !important;
}

/* Make selectbox/multiselect button area use pointer cursor */
div[data-testid="stSelectbox"] button, div[data-testid="stMultiSelect"] button,
div[data-testid="stRadio"] button, div[data-testid="stCheckbox"] button {
  cursor: pointer !important;
  color: VAR_TEXT !important;
  background: transparent !important;
  border: none !important;
}

/* Listbox / options (dropdown) appearance */
ul[role="listbox"] li, div[role="option"] {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border-radius: 6px !important;
}

/* Hover for options */
ul[role="listbox"] li:hover, div[role="option"]:hover {
  background-color: rgba(255,255,255,0.03) !important;
}

/* Selected chips & tags */
.css-1q8dd3e, .css-1v3fvcr, .css-1d391kg {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
}

/* Buttons */
.stButton button {
  background-color: VAR_BTN_BG !important;
  color: VAR_BTN_TEXT !important;
  border: 1px solid VAR_BTN_BORDER !important;
  border-radius: 8px !important;
  cursor: pointer !important;
}

/* Alerts (info/success/warning/error) */
div[data-testid="stAlert"] {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
  border-radius: 10px !important;
  padding: 0.5rem 0.8rem !important;
  box-shadow: none !important;
}
div[data-testid="stAlert"] svg { color: VAR_TEXT !important; opacity: 0.9; }

/* Uploader */
div[data-testid="uploadDropzone"] {
  background-color: VAR_INPUT_BG !important;
  border: 1px dashed VAR_BORDER !important;
  border-radius: 8px !important;
  padding: 0.6rem !important;
}

/* DataFrame styling */
.stDataFrame table { border-collapse: separate !important; border-spacing: 0 6px !important; }
.stDataFrame td, .stDataFrame th {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
  padding: 0.35rem 0.6rem !important;
  border-radius: 6px !important;
}

/* Hide default top header for cleaner look */
header[data-testid="stHeader"] { display: none !important; }

/* Scrollbar */
::-webkit-scrollbar { height: 10px; width: 8px; }
::-webkit-scrollbar-thumb { background: rgba(0,0,0,0.15); border-radius: 8px; }

/* Focus outline (accessible) */
:focus { outline: 2px dashed rgba(255,255,255,0.06) !important; outline-offset: 2px !important; }

</style>
"""

# Color palettes
# ===============================
# COLOR PALETTES (FINAL)
# ===============================

DARK = {
    # Base
    "BG": "#0d1117",            # พื้นหลังหลัก (ไม่ดำสนิท)
    "TEXT": "#e6edf3",          # ตัวอักษรขาวนวล อ่านสบาย

    # Layout
    "SIDEBAR_BG": "#161b22",    # Sidebar โทนเข้มแยกชัด
    "PANEL_BG": "#11161d",      # Panel / Card (สว่างกว่า BG นิดเดียว)
    "BORDER": "#30363d",        # เส้นขอบสุภาพ

    # Code / Data
    "CODE_BG": "#0b0f14",
    "CODE_TEXT": "#d1d7e0",

    # Inputs
    "INPUT_BG": "#0d1117",

    # Buttons
    "BTN_BG": "#238636",        # เขียวสุภาพ (ไม่ดำ ไม่ฉูด)
    "BTN_TEXT": "#ffffff",
    "BTN_BORDER": "#2ea043"
}


LIGHT = {
    # Base
    "BG": "#f9fafb",            # ขาวอมเทา ไม่แสบตา
    "TEXT": "#111827",          # ดำอมฟ้า อ่านชัด

    # Layout
    "SIDEBAR_BG": "#f3f4f6",    # Sidebar อ่อน
    "PANEL_BG": "#ffffff",     # Panel ขาวสะอาด
    "BORDER": "#000000",        # ขอบเทานุ่ม

    # Code / Data
    "CODE_BG": "#f1f5f9",
    "CODE_TEXT": "#111827",

    # Inputs
    "INPUT_BG": "#ffffff",

    # Buttons
    "BTN_BG": "#ffffff",        # น้ำเงินสุภาพ
    "BTN_TEXT": "#ffffff",
    "BTN_BORDER": "#000000"
}


_palette = DARK if theme == T["dark"] else LIGHT

css = THEME_CSS.replace("VAR_BG", _palette["BG"]) \
                .replace("VAR_TEXT", _palette["TEXT"]) \
                .replace("VAR_SIDEBAR_BG", _palette["SIDEBAR_BG"]) \
                .replace("VAR_PANEL_BG", _palette["PANEL_BG"]) \
                .replace("VAR_CODE_BG", _palette["CODE_BG"]) \
                .replace("VAR_CODE_TEXT", _palette["CODE_TEXT"]) \
                .replace("VAR_INPUT_BG", _palette["INPUT_BG"]) \
                .replace("VAR_BORDER", _palette["BORDER"]) \
                .replace("VAR_BTN_BG", _palette["BTN_BG"]) \
                .replace("VAR_BTN_TEXT", _palette["BTN_TEXT"]) \
                .replace("VAR_BTN_BORDER", _palette["BTN_BORDER"])

st.markdown(css, unsafe_allow_html=True)

# ===============================
# App title
# ===============================
st.title(T["title"])

# ===============================
# NCBI & session state
# ===============================
Entrez.email = "poonthakorn@gmail.com"
st.session_state.setdefault("species_map", {})
st.session_state.setdefault("nd5_seqs", {})

# ===============================
# FASTA helpers (allow N)
# ===============================
def clean_seq_text(s: str):
    allowed = set("ACGTN")
    return "".join(ch for ch in s.upper() if ch in allowed)

def read_fasta_text(text: str):
    records = {}
    hdr = None
    seq_parts = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if hdr:
                records[hdr] = clean_seq_text("".join(seq_parts))
            hdr = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if hdr:
        records[hdr] = clean_seq_text("".join(seq_parts))
    return records

# ===============================
# Cached NCBI helpers
# ===============================
@st.cache_data(show_spinner=False)
def fetch_species_from_ncbi(retmax=20):
    try:
        handle = Entrez.esearch(db="nucleotide", term="mitochondrion[Filter] AND complete genome", retmax=retmax)
        rec = Entrez.read(handle)
    except Exception:
        return {}
    species_map = {}
    for uid in rec.get("IdList", []):
        try:
            with Entrez.efetch(db="nucleotide", id=uid, rettype="gb", retmode="text") as h:
                gb = SeqIO.read(h, "genbank")
            name = gb.annotations.get("organism", "Unknown")
            acc = gb.annotations.get("accessions", [""])[0]
            species_map[f"{name} ({acc})"] = acc
        except Exception:
            continue
    return species_map

@st.cache_data(show_spinner=False)
def extract_nd5_from_accession(accession):
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as h:
            record = SeqIO.read(h, "genbank")
        for f in record.features:
            if f.type == "CDS" and "gene" in f.qualifiers:
                if f.qualifiers["gene"][0].lower() == "nd5":
                    return str(f.extract(record.seq))
    except Exception:
        return None
    return None

# ===============================
# Sidebar inputs & fetch list
# ===============================
with st.sidebar:
    st.header(T["load_file"])
    fasta_text = st.text_area(T["paste_label"], height=160)
    uploaded = st.file_uploader(T["upload_label"], type=["fasta", "fa", "txt"])
    if uploaded is not None:
        try:
            raw = uploaded.getvalue().decode("utf-8", "ignore")
        except Exception:
            raw = uploaded.getvalue().decode("latin-1", "ignore")
        fasta_text = raw
        st.success(T["fasta_loaded"].format(n=len(read_fasta_text(raw))))

    if st.button(T["fetch_list"]):
        with st.spinner(T["fetching_species"]):
            species_map = fetch_species_from_ncbi(retmax=20)
            st.session_state["species_map"] = species_map
            st.success(T["fetched_species"].format(n=len(species_map)))

# ===============================
# Read local FASTA
# ===============================
local_records = {}
if 'fasta_text' in locals() and fasta_text and fasta_text.strip():
    local_records = read_fasta_text(fasta_text)
    st.sidebar.success(T["fasta_loaded"].format(n=len(local_records)))

# ===============================
# Selection UI
# ===============================
selected = []
col1, col2 = st.columns(2)

with col1:
    if local_records:
        st.markdown(f"**{T['load_file']} (จากข้อความ/ไฟล์)**")
        chosen_local = st.multiselect(T["species_select"] + " (FASTA)", list(local_records.keys()))
        selected.extend(chosen_local)

with col2:
    if st.session_state.get("species_map"):
        st.markdown(f"**{T['fetch_list']}**")
        chosen_ncbi = st.multiselect(T["species_select"] + " (NCBI)", list(st.session_state["species_map"].keys()))
        selected.extend(chosen_ncbi)

if not selected:
    st.info(T["select_prompt"])

# ===============================
# Fetch ND5
# ===============================
if st.button(T["fetch_nd5"]):
    with st.spinner(T["fetching_nd5"]):
        for key in selected:
            if key in local_records:
                st.session_state["nd5_seqs"][key] = local_records[key]
            else:
                acc = st.session_state.get("species_map", {}).get(key)
                if acc:
                    seq = extract_nd5_from_accession(acc)
                    if seq:
                        st.session_state["nd5_seqs"][key] = seq
                    else:
                        st.warning(f"{T['failed_fetch']}: {key}")
                else:
                    st.warning(f"{T['failed_fetch']}: {key}")
    st.success(T["fetched"].format(n=len(st.session_state.get("nd5_seqs", {}))))

# ===============================
# Identity calculation (pairwise2 globalxx, count only positions where both non-gap)
# ===============================
def identity_percent(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    alns = pairwise2.align.globalxx(a, b, one_alignment_only=True)
    if not alns:
        return 0.0
    aligned_a, aligned_b = alns[0][0], alns[0][1]
    matches = 0
    positions = 0
    for x, y in zip(aligned_a, aligned_b):
        if x != "-" and y != "-":
            positions += 1
            if x == y:
                matches += 1
    if positions == 0:
        return 0.0
    return matches / positions * 100.0

# ===============================
# Show identity matrix + heatmap + download
# ===============================
if st.button(T["identity"]):
    seqs = st.session_state.get("nd5_seqs", {})
    if len(seqs) < 2:
        st.warning(T["need_two"])
    else:
        names = list(seqs.keys())
        matrix = np.zeros((len(names), len(names)))
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                if i == j:
                    matrix[i, j] = 100.0
                else:
                    matrix[i, j] = round(identity_percent(seqs[a], seqs[b]), 2)

        # display names (localized common names for Thai)
        COMMON_NAMES = {
            "Homo sapiens": "มนุษย์",
            "Pan troglodytes": "ชิมแพนซี",
            "Gorilla gorilla": "กอริลลา",
            "Pan paniscus": "โบนาโบ",
            "Pongo abelii": "อุรังอุตัง",
        }
        display_names = []
        for n in names:
            sci = n.split(" (")[0]
            display_names.append(COMMON_NAMES.get(sci, sci) if language == "ภาษาไทย" else sci)

        df = pd.DataFrame(matrix, index=display_names, columns=display_names)
        st.dataframe(df)

        csv_bytes = df.to_csv().encode("utf-8")
        st.download_button(T["download_csv"], data=csv_bytes, file_name="nd5_identity_matrix.csv", mime="text/csv")

        # heatmap
        fig, ax = plt.subplots(figsize=(6, 4))
        im = ax.imshow(matrix, aspect='auto')
        ax.set_xticks(np.arange(len(display_names)))
        ax.set_yticks(np.arange(len(display_names)))
        ax.set_xticklabels(display_names, rotation=45, ha='right')
        ax.set_yticklabels(display_names)
        ax.set_title(T["identity"])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        st.pyplot(fig)

# ===============================
# Theory / educational content
# ===============================
THEORY = {
    "ภาษาไทย": [
        ("🔬 กลไกระดับโมเลกุลของ DNA (ระดับลึก)",
         """
**DNA ถูกบรรจุเป็นโครมาติน (Chromatin)**  
โดยมี **นิวคลีโอโซม (Nucleosome)** เป็นหน่วยพื้นฐาน  
ประกอบด้วย DNA ประมาณ **147 เบส** พันรอบฮิสโตน **8 ตัว**

### สภาพโครมาติน
- **ยูโครมาติน (Euchromatin)** — จัดเรียงหลวม → เข้าถึงง่าย → ถอดรหัสได้  
- **เฮเทอโรโครมาติน (Heterochromatin)** — จัดเรียงแน่น → ปิดการถอดรหัส

### การดัดแปลงฮิสโตน
- **Acetylation** (เช่น H3K9ac, H3K27ac) → เปิดโครมาติน → ส่งเสริมการถอดรหัส  
- **Methylation** (เช่น H3K27me3) → โครมาตินแน่น → ยับยั้งการถอดรหัส

### DNA methylation
- เกิดที่ **CpG sites** → recruit โปรตีนยับยั้ง → ปิด promoter → ลดการแสดงออกของยีน
"""),
        ("🧠 กลไกระดับโมเลกุลของ DNA (ภาษาง่าย)",
         """
ลองนึกว่า **DNA = หนังสือคู่มือชีวิต**  
- DNA พันรอบฮิสโตน เหมือนด้ายพันแกน  
- พันแน่น → อ่านไม่ออก (เฮเทอโรโครมาติน)  
- พันหลวม → อ่านง่าย (ยูโครมาติน)

**เครื่องจักรอ่าน**: RNA polymerase II ทำงานร่วมกับ transcription factors ที่ promoter/enhancer เพื่อเริ่ม transcription → mRNA → ตัดต่อ → แปลเป็นโปรตีน

**RNA ที่ไม่เข้ารหัส** (microRNA, lncRNA) ทำหน้าที่ปรับหรือยับยั้งการแปลหลังการถอดรหัส
"""),
    ],
    "English": [
        ("🔬 Molecular mechanisms of DNA (detailed)",
         """
**DNA is packaged as chromatin.**  
The basic unit is the **nucleosome**: ~147 bp of DNA wrapped around an octamer of histones.

### Chromatin states
- **Euchromatin** — open, accessible, transcriptionally active.  
- **Heterochromatin** — compact, transcriptionally repressed.

### Histone modifications
- **Acetylation** (e.g. H3K9ac, H3K27ac) → opens chromatin and promotes transcription.  
- **Methylation** (e.g. H3K27me3) → can compact chromatin and repress transcription.

### DNA methylation
- Occurs at **CpG** sites, recruits repressive proteins, and often silences promoters.
"""),
        ("🧠 Molecular mechanisms of DNA (plain language)",
         """
Imagine **DNA as a huge instruction manual**.  
DNA is wrapped around histone 'spools' — if wrapped tightly, pages are hard to read (heterochromatin); if loosely, they're easy to read (euchromatin).

**RNA polymerase II** and transcription factors bind promoters/enhancers to start transcription → pre-mRNA → splicing → mRNA → translation.

Non-coding RNAs (microRNA, lncRNA) regulate gene expression post-transcriptionally.
"""),
    ]
}

st.markdown("## 🧬 ความรู้ชีววิทยาระดับโมเลกุล & วิวัฒนาการ" if language == "ภาษาไทย" else "## 🧬 Molecular biology & evolution")
for title, content in THEORY[language]:
    with st.expander(title, expanded=False):
        st.markdown(content)
