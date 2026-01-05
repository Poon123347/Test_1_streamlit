import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from io import StringIO

# ===============================
# ตั้งค่าเพจ
# ===============================
st.set_page_config(page_title="DNA & ND5 Explorer", layout="wide")

# ===============================
# พจนานุกรมภาษา (ไทย + อังกฤษ)
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
        "fetched": "ดึง ND5 ได้",
        "fasta_loaded": "โหลด FASTA แล้ว",
        "no_local_fasta": "ยังไม่มีไฟล์ FASTA ในอินพุต — วางข้อความหรืออัปโหลดได้",
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
        "fetched": "Fetched ND5",
        "fasta_loaded": "FASTA loaded",
        "no_local_fasta": "No FASTA input — paste or upload",
    }
}

# ===============================
# Sidebar: language + theme
# ===============================
with st.sidebar:
    language = st.selectbox("🌐 Language / ภาษา", ["ภาษาไทย", "English"], index=0)
    theme = st.selectbox(LANG[language]["theme"], [LANG[language]["light"], LANG[language]["dark"]], index=0)

T = LANG[language]

# ===============================
# Theme CSS template (will replace placeholders) & hide top header
# ===============================
THEME_CSS = """
<style>
html, body, .stApp, .block-container, section, div {
  transition: background-color 0.45s ease, color 0.45s ease, border-color 0.45s ease;
}
.stApp, .reportview-container, .block-container, main {
  background-color: VAR_BG !important;
  color: VAR_TEXT !important;
}
section[data-testid="stSidebar"] {
  background-color: VAR_SIDEBAR_BG !important;
  color: VAR_TEXT !important;
  border-right: 1px solid VAR_BORDER !important;
}
div[data-testid="stVerticalBlock"], div[data-testid="stMarkdownContainer"],
div[data-testid="stCodeBlock"], .stDataFrame, .css-1dp5vir, .stExpander {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
  border-radius: 10px !important;
}
pre, code, div[data-testid="stCodeBlock"] pre, div[data-testid="stCodeBlock"] code {
  background-color: VAR_CODE_BG !important;
  color: VAR_CODE_TEXT !important;
  border: 1px solid VAR_BORDER !important;
}
input, textarea, select, .stTextInput, .stTextArea, .stSelectbox, .stMultiSelect {
  background-color: VAR_INPUT_BG !important;
  color: VAR_TEXT !important;
  border: 1px solid VAR_BORDER !important;
}
.stButton button {
  background-color: VAR_BTN_BG !important;
  color: VAR_BTN_TEXT !important;
  border: 1px solid VAR_BTN_BORDER !important;
  border-radius: 8px !important;
}
.css-1q8dd3e, .css-1d391kg, .css-1v3fvcr {
  background-color: VAR_PANEL_BG !important;
  color: VAR_TEXT !important;
  border-color: VAR_BORDER !important;
}
/* hide Streamlit top header (toolbar) for a clean look */
header[data-testid="stHeader"] {
  display: none !important;
}
.stApp {
  padding-top: 0 !important;
}
::-webkit-scrollbar { height: 10px; width: 8px; }
::-webkit-scrollbar-thumb { background: rgba(0,0,0,0.15); border-radius: 8px; }
</style>
"""

# Dark / Light palettes
DARK = {
  "BG": "#0e1117","TEXT": "#e6e6e6","SIDEBAR_BG": "#161b22","PANEL_BG": "#0f1419",
  "CODE_BG": "#0b0f13","CODE_TEXT": "#e6e6e6","INPUT_BG": "#0d1117","BORDER": "#30363d",
  "BTN_BG": "#11632d","BTN_TEXT": "#ffffff","BTN_BORDER": "#238636"
}
LIGHT = {
  "BG": "#ffffff","TEXT": "#111111","SIDEBAR_BG": "#f6f6f6","PANEL_BG": "#fbfbfb",
  "CODE_BG": "#f5f5f5","CODE_TEXT": "#111111","INPUT_BG": "#ffffff","BORDER": "#e6e6e6",
  "BTN_BG": "#f0f0f0","BTN_TEXT": "#111111","BTN_BORDER": "#d0d0d0"
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
# Title
# ===============================
st.title(T["title"])

# ===============================
# NCBI setup & session state
# ===============================
Entrez.email = "poonthakorn@gmail.com"
st.session_state.setdefault("species_map", {})   # map: "Sci (ACC)" -> ACC
st.session_state.setdefault("nd5_seqs", {})     # map: label -> seq (string)

# ===============================
# PairwiseAligner
# ===============================
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1.0
aligner.mismatch_score = 0.0
aligner.open_gap_score = 0.0
aligner.extend_gap_score = 0.0

# ===============================
# FASTA helpers
# ===============================
def clean_seq_text(s: str):
    return "".join(ch for ch in s.upper() if ch in "ACGT")

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
# Sidebar input & fetch list
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
        st.success(f"{T['fasta_loaded']}: {len(read_fasta_text(raw))}")

    if st.button(T["fetch_list"]):
        with st.spinner(T["fetching_species"]):
            try:
                query = "mitochondrion[Filter] AND complete genome"
                handle = Entrez.esearch(db="nucleotide", term=query, retmax=20)
                rec = Entrez.read(handle)
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
                st.session_state["species_map"] = species_map
                st.success(f"ได้รายการ {len(species_map)} สายพันธุ์จาก NCBI")
            except Exception as e:
                st.warning(f"{T['failed_fetch']}: {e}")

# ===============================
# Read FASTA from paste/upload
# ===============================
local_records = {}
if 'fasta_text' in locals() and fasta_text and fasta_text.strip():
    local_records = read_fasta_text(fasta_text)
    st.sidebar.success(f"{T['fasta_loaded']}: {len(local_records)}")

# ===============================
# Species selection UI
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
# ND5 fetch helpers
# ===============================
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

def guess_species_from_header(header: str):
    parts = header.split()
    if len(parts) >= 2:
        return parts[0] + " " + parts[1]
    return parts[0] if parts else header

# ===============================
# Fetch ND5 button behavior
# ===============================
if st.button(T["fetch_nd5"]):
    nd5_seqs = st.session_state.get("nd5_seqs", {})
    with st.spinner(T["fetching_nd5"]):
        for key in selected:
            if key in local_records:
                species_name = guess_species_from_header(key)
                seq = None
                try:
                    handle = Entrez.esearch(db="nucleotide", term=f"{species_name}[Organism] AND mitochondrion[Filter] AND complete genome", retmax=1)
                    rec = Entrez.read(handle)
                    if rec.get("IdList"):
                        uid = rec["IdList"][0]
                        with Entrez.efetch(db="nucleotide", id=uid, rettype="gb", retmode="text") as h2:
                            record = SeqIO.read(h2, "genbank")
                        for f in record.features:
                            if f.type == "CDS" and "gene" in f.qualifiers and f.qualifiers["gene"][0].lower() == "nd5":
                                seq = str(f.extract(record.seq))
                except Exception:
                    seq = None
                if seq:
                    st.session_state["nd5_seqs"][key] = seq
                else:
                    st.warning(f"{T['failed_fetch']}: {key}")
            else:
                acc = st.session_state["species_map"].get(key)
                if acc:
                    seq = extract_nd5_from_accession(acc)
                    if seq:
                        st.session_state["nd5_seqs"][key] = seq
                    else:
                        st.warning(f"{T['failed_fetch']}: {key}")
    st.success(f"{T['fetched']}: {len(st.session_state.get('nd5_seqs', {}))}")

# ===============================
# Identity calculation using PairwiseAligner
# ===============================
def global_identity(s1: str, s2: str) -> float:
    if not s1 or not s2:
        return 0.0
    try:
        aln = aligner.align(s1, s2)[0]
    except Exception:
        return 0.0
    aligned_blocks_1, aligned_blocks_2 = aln.aligned
    matches = 0
    total_aligned = 0
    for (a_start, a_end), (b_start, b_end) in zip(aligned_blocks_1, aligned_blocks_2):
        length = min(a_end - a_start, b_end - b_start)
        for i in range(length):
            if s1[a_start + i] == s2[b_start + i]:
                matches += 1
            total_aligned += 1
    denom = max(len(s1), len(s2))
    return (matches / denom * 100.0) if denom > 0 else 0.0

# ===============================
# Show identity matrix
# ===============================
if st.button(T["identity"]):
    seqs = st.session_state.get("nd5_seqs", {})
    if len(seqs) < 2:
        st.warning(T["need_two"])
    else:
        names = list(seqs.keys())
        matrix = []
        for i, a in enumerate(names):
            row = []
            for j, b in enumerate(names):
                if i == j:
                    row.append(100.0)
                else:
                    pid = global_identity(seqs[a], seqs[b])
                    row.append(round(pid, 2))
            matrix.append(row)
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

# ===============================
# Theory content (localized)
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
        ("🧬 ความแตกต่างของ DNA มนุษย์กับสัตว์อื่น (ระดับลึก)",
         """
จีโนมของชนิดต่าง ๆ ต่างกันทั้งลำดับและการจัดวาง แม้หลายยีนจะอนุรักษ์ เช่น HOX genes  
- มนุษย์–ชิมแปนซี ~98–99% แต่ความต่างเล็ก ๆ สามารถทำให้เกิดความแตกต่างทางสรีรวิทยาได้มาก  
- หนูมี synteny กับมนุษย์สูง แต่แยกสายกว่า (75–90 ล้านปี)  
- ไก่มีความเหมือนยีนพื้นฐานกับมนุษย์ ~60%
"""),
        ("🐒 ความแตกต่างของ DNA (ภาษาง่าย)",
         """
DNA เปรียบเหมือนพิมพ์เขียว: แม้บางสายพันธุ์จะเหมือนกันมาก แต่รายละเอียดเล็ก ๆ (ฟุตโน้ต, enhancer, non-coding RNA) ทำให้ผลลัพธ์แตกต่าง
"""),
        ("🌳 ความสัมพันธ์เชิงวิวัฒนาการ (ระดับลึก)",
         """
### Molecular clock & phylogeny
- มนุษย์–ชิมแปนซี: แยก ~6–7 ล้านปี  
- มนุษย์–หนู: แยก ~75–90 ล้านปี  
- มนุษย์–ไก่: แยก >300 ล้านปี

mtDNA (เช่น ND5) เปลี่ยนช้าจึงใช้สืบสายวิวัฒนาการได้ดี
"""),
        ("📚 แหล่งอ้างอิง",
         """
- NCBI Bookshelf – Chromatin & Epigenetics  
  https://www.ncbi.nlm.nih.gov/books/NBK532999/

- Enhancer / Silencer  
  https://www.ncbi.nlm.nih.gov/books/NBK459456/

- Human vs Chimp DNA (livescience)  
  https://www.livescience.com/archaeology/human-evolution/do-humans-and-chimps-really-share-nearly-99-percent-of-their-dna

- Mouse–Human Genomics (PMC)  
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6413734/

- Chicken–Human Genome (Genome.gov)  
  https://www.genome.gov/12514316
""")
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
        ("🧬 Differences between human and other animals (detailed)",
         """
Genomes differ in sequence and organization. Many core genes (e.g., HOX) are conserved, but regulatory differences and small sequence changes can yield big phenotypic differences.

- Human vs Chimp: ~98–99% similarity; split ~6–7 million years ago.  
- Mouse: conserved synteny but diverged ~75–90 million years.  
- Chicken: ~60% gene similarity to human in many regions.
"""),
        ("🐒 Differences (plain language)",
         """
DNA is like a blueprint. Small edits (a few letters or regulatory notes) can create large differences between species, even when overall sequences look very similar.
"""),
        ("🌳 Evolutionary relationships (detailed)",
         """
### Molecular clock & phylogeny
- Human–Chimp: ~6–7 million years  
- Human–Mouse: ~75–90 million years  
- Human–Chicken: >300 million years

mtDNA genes (e.g., ND5) evolve relatively slowly and are useful for phylogenetic inference.
"""),
        ("📚 References",
         """
- NCBI Bookshelf – Chromatin & Epigenetics  
  https://www.ncbi.nlm.nih.gov/books/NBK532999/

- Enhancer / Silencer  
  https://www.ncbi.nlm.nih.gov/books/NBK459456/

- Human vs Chimp DNA (livescience)  
  https://www.livescience.com/archaeology/human-evolution/do-humans-and-chimps-really-share-nearly-99-percent-of-their-dna

- Mouse–Human Genomics (PMC)  
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6413734/

- Chicken–Human Genome (Genome.gov)  
  https://www.genome.gov/12514316
""")
    ]
}

st.markdown("## 🧬 ความรู้ชีววิทยาระดับโมเลกุล & วิวัฒนาการ" if language == "ภาษาไทย" else "## 🧬 Molecular biology & evolution")

for title, content in THEORY[language]:
    with st.expander(title, expanded=False):
        st.markdown(content)

# End of file
