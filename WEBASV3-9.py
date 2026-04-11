import matplotlib
matplotlib.use("Agg")  # kill my prob
import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from io import StringIO
from io import BytesIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from copy import deepcopy
from matplotlib.collections import LineCollection
plt.rcParams["toolbar"] = "None"
from matplotlib import font_manager, rcParams
import re
import os


FONT_PATH = os.path.abspath("fonts/NotoSansThai-Regular.ttf")

if not os.path.exists(FONT_PATH):
    raise FileNotFoundError(f"Font not found: {FONT_PATH}")

font_manager.fontManager.addfont(FONT_PATH)
font_prop = font_manager.FontProperties(fname=FONT_PATH)

rcParams["font.family"] = font_prop.get_name()
rcParams["axes.unicode_minus"] = False

# ========== PAGE CONFIG ==========
st.set_page_config(
    page_title="ND5 Phylogenetic Explorer",
    page_icon="🧬",
    layout="centered",
    initial_sidebar_state="collapsed"  # 👈 IMPORTANT
)

language = st.selectbox("Language", ["English", "ภาษาไทย"])

DARK = {
    "BG": "#0d1117",
    "TEXT": "#e6edf3",
    "TEXT_MUTED": "#9ca3af",
    "SIDEBAR_BG": "#161b22",
    "PANEL_BG": "#11161d",
    "BORDER": "#30363d",
    "CODE_BG": "#0b0f14",
    "CODE_TEXT": "#d1d7e0",
    "BTN_BG": "#238636",
    "BTN_TEXT": "#ffffff",
    "TREE_BORDER_BG" : "#9ca3af",
    "TABLE_TEXT" : "#535353"
}

THEME = DARK 
st.markdown('<div style="height:48px"></div>', unsafe_allow_html=True)

gdg_accent  = "#238636"
gdg_sel     = "rgba(35, 134, 54, 0.30)"
gdg_region  = "rgba(35, 134, 54, 0.18)"
gdg_overlay = "rgba(35, 134, 54, 0.12)"
gdg_hover   = "rgba(35, 134, 54, 0.08)"

st.markdown(f"""
<style>

/* ================== THEME VARIABLES ================== */
:root {{
    --bg: {THEME["BG"]};
    --text: {THEME["TEXT"]};
    --text-muted: {THEME["TEXT_MUTED"]};
    --sidebar-bg: {THEME["SIDEBAR_BG"]};
    --panel-bg: {THEME["PANEL_BG"]};
    --border: {THEME["BORDER"]};
    --btn-bg: {THEME["BTN_BG"]};
    --btn-text: {THEME["BTN_TEXT"]};
    --code-bg: {THEME["CODE_BG"]};
    --code-text: {THEME["CODE_TEXT"]};
}}

/* ================== GLOBAL ================== */
.stApp {{
    background-color: var(--bg) !important;
    color: var(--text) !important;
    padding-top: 0 !important;
}}

html, body, div, span, p, label,
h1, h2, h3, h4, h5, h6 {{
    color: var(--text) !important;
}}

header[data-testid="stHeader"] {{
    display: block !important;
    background: var(--bg) !important;
}}

/* ================== BUTTONS ================== */
.stButton button,
.stDownloadButton button {{
    background-color: var(--btn-bg) !important;
    color: var(--btn-text) !important;
    border-radius: 12px !important;
    padding: 0.45rem 1.2rem !important;
    border: none !important;
}}

div[data-testid="stButton"] > button {{
    width: 300% !important;
    min-width: 160px !important;
    box-sizing: border-box !important;
}}

.stDownloadButton svg {{
    fill: var(--btn-text) !important;
}}

/* ================== EXPANDER ================== */
div[data-testid="stExpander"] > details > summary {{
    background-color: var(--panel-bg) !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
    padding: 12px !important;
}}

div[data-testid="stExpander"] > details > div[role="region"] {{
    background-color: var(--panel-bg) !important;
    border: 1px solid var(--border) !important;
    border-top: none !important;
    border-radius: 0 0 12px 12px !important;
    padding: 16px !important;
}}

/* ================== SELECT ================== */
div[data-baseweb="select"] > div {{
    background-color: var(--panel-bg) !important;
    border: 1px solid var(--border) !important;
    min-height: 48px !important;
    padding: 4px 8px !important;
}}

div[data-baseweb="select"] span,
div[data-baseweb="select"] input {{
    color: var(--text) !important;
}}

li[role="option"] {{
    background-color: var(--panel-bg) !important;
    color: var(--text) !important;
}}

/* ================== DATAFRAME ================== */
div[data-testid="stDataFrame"] {{
    width: 300% !important;
    overflow-x: auto !important;
}}

/* ================== METRICS ================== */
div[data-testid="stMetricLabel"],
div[data-testid="stMetricValue"] {{
    color: var(--text) !important;
}}

/* ================== CODE ================== */
pre, code {{
    background-color: var(--code-bg) !important;
    color: var(--code-text) !important;
    border-radius: 10px !important;
}}

/* ================== SVG ================== */
svg text,
svg tspan {{
    fill: var(--text) !important;
}}

/* ================== DATAFRAME TEXT ================== */
div[data-testid="stDataFrame"] * {{
    color: var(--text) !important;
}}
div[data-testid="stDataFrame"] div[role="columnheader"] {{
    color: var(--text) !important;
    font-weight: 600 !important;
}}
div[data-testid="stDataFrame"] div[role="gridcell"] {{
    background-color: var(--panel-bg) !important;
}}

/* ================== ARROWS / ICONS ================== */
div[data-testid="stExpander"] summary svg {{
    fill: var(--text) !important;
}}
div[data-baseweb="select"] svg {{
    fill: var(--text) !important;
}}

/* ================== GLIDE DATAGRID SELECTION (THEME-AWARE) ================== */
div[data-testid="stDataFrame"] {{
    --gdg-accent-color:           {gdg_accent}  !important;
    --gdg-selection-color:        {gdg_sel}     !important;
    --gdg-selection-region-color: {gdg_region}  !important;
    --gdg-selection-overlay-color:{gdg_overlay} !important;
}}

div[data-testid="stDataFrame"] div[role="gridcell"]:hover {{
    background-color: {gdg_hover} !important;
}}

div[data-testid="stDataFrame"] .gdg-selection,
div[data-testid="stDataFrame"] .selection,
div[data-testid="stDataFrame"] [class*="selection"] {{
    background: {gdg_sel} !important;
    outline-color: {gdg_accent} !important;
}}

</style>
""", unsafe_allow_html=True)

# ========== NCBI CONFIG ==========
Entrez.email = "poonthakorn@gmail.com"

species_accessions = {
    # === Primates ===
    "Human": "NC_012920.1",
    "Neanderthal": "NC_011137.1",       # ancient human
    "Chimpanzee": "NC_001643.1",
    "Bonobo": "NC_001644.1",
    "Eastern Gorilla": "KF914213.1",
    "Western Gorilla": "NC_011120.1",
    "Bornean Orangutan": "NC_001646.1",
    "Sumatran Orangutan": "NC_002083.1",
    "Common Gibbon": "NC_002082.1",
    "Rhesus Macaque": "NC_005943.1",

    # === Old World monkeys ===
    "Crab-eating Macaque": "NC_012670.1",
    "Barbary Ape": "NC_002764.1",
    "Green Monkey": "NC_008066.1",
    "Tantalus Monkey": "NC_009748.1",
    "Dusky Leaf Monkey": "NC_006900.1",
    "Proboscis Monkey": "NC_008216.1",
    "Red-shanked Douc": "NC_008220.1",
    "Tonkin Snub-nose Monkey": "NC_015485.1",
    "Colobus Monkey": "NC_006901.1",

    # === Other Mammals ===
    "Pig": "NC_000845.1",
    "Cow": "NC_001567.1",
    "Horse": "NC_001640.1",
    "Dog": "NC_002008.4",
    "Cat": "NC_001700.1",
    "Sheep": "NC_001941.1",
    "Goat": "NC_005044.2",
    "Mouse": "NC_005089.1",
    "Rat": "NC_001665.2",
    "Rabbit": "NC_001913.1",

    # === Additional Mammals ===
    "Elephant": "NC_000934.1",
    "Whale": "NC_006853.1",
    "Lion": "NC_028302.1",
    "Tiger": "NC_010642.1",
    "Brown Bear": "NC_003427.1",
    "Dolphin": "NC_012057.1",
    "Seal": "NC_008416.1",

    # === Birds ===
    "Chicken": "NC_001323.1",
    "Duck": "NC_009684.1",
    "Goose": "NC_003049.1",
    "Turkey": "NC_010195.1",

    # === Fish ===
    "Zebrafish": "NC_002333.2",
    "Salmon": "NC_001960.1",
    "Tilapia": "NC_013723.1",
    "Catfish": "NC_015644.1",
    "Bluefin Tuna": "NC_010383.1",

    # === Reptiles / Amphibians ===
    "Frog": "NC_002805.1",
    "Lizard": "NC_018045.1",

    # === Invertebrates ===
    "Fruit Fly": "NC_024511.2",
    "Bee": "NC_001566.1",
    "Octopus": "NC_006353.1",
    "Roundworms": "NC_001328.1"
}
SPECIES_TH = {
    "Human": "มนุษย์",
    "Neanderthal": "นีแอนเดอร์ธัล",
    "Chimpanzee": "ชิมแปนซี",
    "Bonobo": "โบโนโบ",
    "Eastern Gorilla": "กอริลลาตะวันออก",
    "Western Gorilla": "กอริลลาตะวันตก",
    "Bornean Orangutan": "อุรังอุตังบอร์เนียว",
    "Sumatran Orangutan": "อุรังอุตังสุมาตรา",
    "Common Gibbon": "ชะนี",
    "Rhesus Macaque": "ลิงแสม",
    "Crab-eating Macaque": "ลิงแสมหางยาว",
    "Barbary Ape": "ลิงบาบารี",
    "Green Monkey": "ลิงเขียว",
    "Tantalus Monkey": "ลิงแทนทาลัส",
    "Dusky Leaf Monkey": "ค่างแว่น",
    "Proboscis Monkey": "ลิงจมูกยาว",
    "Red-shanked Douc": "ค่างห้าสี",
    "Tonkin Snub-nose Monkey": "ลิงจมูกเชิดตองกิ้น",
    "Colobus Monkey": "ลิงโคลอบัส",
    "Pig": "หมู",
    "Cow": "วัว",
    "Horse": "ม้า",
    "Dog": "สุนัข",
    "Cat": "แมว",
    "Sheep": "แกะ",
    "Goat": "แพะ",
    "Mouse": "หนู",
    "Rat": "หนูท่อ",
    "Rabbit": "กระต่าย",
    "Elephant": "ช้าง",
    "Whale": "วาฬ",
    "Lion": "สิงโต",
    "Tiger": "เสือโคร่ง",
    "Brown Bear": "หมีสีน้ำตาล",
    "Dolphin": "โลมา",
    "Seal": "แมวน้ำ",
    "Chicken": "ไก่",
    "Duck": "เป็ด",
    "Goose": "ห่าน",
    "Turkey": "ไก่งวง",
    "Zebrafish": "ปลาม้าลาย",
    "Salmon": "ปลาแซลมอน",
    "Tilapia": "ปลาตะเพียน",
    "Catfish": "ปลาดุก",
    "Bluefin Tuna": "ปลาทูน่าครีบใหญ่",
    "Frog": "กบ",
    "Lizard": "กิ้งก่า",
    "Fruit Fly": "แมลงวันผลไม้",
    "Bee": "ผึ้ง",
    "Octopus": "หมึกสาย",
    "Roundworms": "พยาธิตัวกลม"
}


# ========== LANGUAGE SUPPORT ==========
LANGS = {
    "ภาษาไทย": {
        # ===== Header =====
        "title": "ตัวสำรวจสายวิวัฒนาการ ND5",
        "subtitle": "เปรียบเทียบ ND5 ระหว่างสายพันธุ์ – แผนภูมิ + ตารางความเหมือน",

        # ===== Controls =====
        "select_species": "กรุณาเลือกสายพันธุ์ตัวอย่าง (อย่างน้อย 2 ชนิด)",
        "choose_options": "เลือกสายพันธุ์",
        "available": "สายพันธุ์ทั้งหมด",
        "selected": "สายพันธุ์ที่เลือก",
        "fetch": "ดึงข้อมูล",

        # ===== Info =====
        "length": "ความยาว ND5 (bp)",
        "preview": "ดูตัวอย่างลำดับเบส",

        # ===== Tabs =====
        "results": "ผลลัพธ์",
        "tree": "ต้นไม้วิวัฒนาการ",
        "methods": "วิธีการและพื้นฐานทางทฤษฎี",

        # ===== Outputs =====
        "identity": "ตารางความเหมือน (%)",
        "phylo_tree": "ต้นไม้วิวัฒนาการ (ND5)",

        # ===== States =====
        "spinner": "กำลังดึงข้อมูล...",
        "no_results": "กรุณาดึงข้อมูล ND5 เพื่อดูผลลัพธ์",
        "no_tree": "กรุณาดึงข้อมูล ND5 ก่อน",

        # ===== Messages =====
        "missing_nd5": "⚠️ สายพันธุ์ต่อไปนี้ไม่มีข้อมูล ND5 ในฐานข้อมูล NCBI:",

        # ===== Expanders =====
        "methods_expander": "🔬 กลไกระดับโมเลกุล ขั้นตอน และวัตถุประสงค์",

        "species" : "สายพันธุ์"
    },

    "English": {
        # ===== Header =====
        "title": "ND5 Phylogenetic Explorer",
        "subtitle": "Compare ND5 across species – tree + identity matrix",

        # ===== Controls =====
        "select_species": "Please select species (min 2)",
        "choose_options": "Choose species",
        "available": "Available species",
        "selected": "Selected",
        "fetch": "Fetch",

        # ===== Info =====
        "length": "ND5 length (bp)",
        "preview": "Preview sequences",

        # ===== Tabs =====
        "results": "Results",
        "tree": "Phylogenetic Tree",
        "methods": "Methods & Background",

        # ===== Outputs =====
        "identity": "Identity Matrix (%)",
        "phylo_tree": "Phylogenetic Tree (ND5)",

        # ===== States =====
        "spinner": "Fetching data...",
        "no_results": "Fetch ND5 to see results.",
        "no_tree": "Fetch ND5 to generate tree.",

        # ===== Messages =====
        "missing_nd5": "⚠️ The following species do not have ND5 data in NCBI:",

        # ===== Expanders =====
        "methods_expander": "🔬 Molecular mechanisms, pipeline and objectives",

        "species" : "species"
    }
}

# ===============================
# METHODS TAB
# ===============================
METHODS_TEXT = {
    "English": """
Molecular Mechanisms of DNA (In-depth Overview)

DNA Packaging and Chromatin Structure  
Eukaryotic DNA is packaged into chromatin, whose fundamental structural unit is the nucleosome.
Each nucleosome consists of approximately 146–147 base pairs of DNA wrapped around an octamer of histone proteins.
Chromatin exists in two major states: euchromatin, which is loosely packed and transcriptionally active,
and heterochromatin, which is tightly packed and transcriptionally repressive.

Histone Modifications  
Histone protein tails can undergo chemical modifications such as acetylation and methylation,
which affect chromatin structure and regulate gene expression.
For example, histone acetylation at H3K27ac is commonly associated with actively transcribed regions,
whereas histone methylation can either activate or repress transcription depending on the modified residue and context.

DNA Methylation (CpG Islands)  
DNA methylation typically occurs at cytosine residues within CpG dinucleotides,
particularly in CpG islands near gene promoter regions.
Methylation in these regions can inhibit transcription factor binding
and lead to gene silencing, serving as an important epigenetic regulatory mechanism.

Chromatin Remodeling  
Chromatin remodeling complexes use ATP to reposition or remove nucleosomes,
thereby regulating DNA accessibility to transcriptional machinery.

Transcription and Cis-Regulatory Elements  
Transcription begins at promoter regions through the binding of RNA polymerase II
and transcription factors.
Enhancers increase transcriptional activity, whereas silencers suppress gene expression.

Epigenetics  
Epigenetics refers to heritable changes in gene expression that do not involve
alterations in the DNA sequence.
Major epigenetic mechanisms include DNA methylation and histone modifications,
which define cell-type-specific gene expression patterns.

Non-coding DNA and RNA  
Although most of the human genome does not encode proteins,
non-coding elements such as microRNAs and long non-coding RNAs
play essential roles in post-transcriptional gene regulation.

---

Human DNA Analysis Using Bioinformatics

Definition of Bioinformatics  
Bioinformatics is the application of computational tools and algorithms
to analyze biological sequence data, enabling the study of genetic variation,
gene function, evolution, and disease-associated mutations.

---

Objectives  
1 To study human DNA in order to understand human lineage and genetic relationships  
2 To construct phylogenetic trees to examine the relationships between human DNA and animal genetic variation and to study the functional roles of DNA

Expected Benefits  
1 To enhance understanding of human lineage and relationships through DNA analysis  
2 To gain insight into the origin and evolutionary development of humans and provide a foundation for further study in medical research

---

Bioinformatics Analysis Pipeline

Data Source  
DNA sequences used in this study are obtained from the public database NCBI GenBank
using accession numbers for each species.
All data are publicly available and contain no personal or sensitive information.

Data Preparation and Quality Control  
1 Extraction of the ND5 (MT-ND5) gene from mitochondrial genomes  
2 ND5 gene length varies slightly among species such as 1812 bp 1815 bp or 1821 bp  
3 These differences arise from evolutionary insertions and deletions differences in start or stop codon annotation and species specific mutations accumulated over time  
4 Such variation is biologically normal and reflects evolutionary divergence rather than sequencing error  
5 Sequences with excessive ambiguous bases or internal stop codons are excluded

Multiple Sequence Alignment  
Multiple sequence alignment is performed to ensure homologous nucleotide positions
are accurately compared across species.

Evolutionary Distance and Phylogenetic Tree Construction  
1 Pairwise sequence similarity and evolutionary distances are calculated  
2 A phylogenetic tree is constructed to visualize evolutionary relationships based on ND5 sequence variation

---

Biological Function of the ND5 Gene  
The mitochondrial MT-ND5 gene encodes a subunit of Complex I
in the mitochondrial electron transport chain,
which plays a critical role in ATP production.

Medical Significance  
Mutations in the ND5 gene are associated with mitochondrial disorders
such as MELAS and Leigh syndrome,
highlighting its importance in both evolutionary biology and medical genetics.

Reproducibility and Future Applications  
The analytical workflow can be reproduced and extended to other datasets,
providing a foundation for future studies in evolutionary biology,
population genetics, and biomedical research.
""",

    "ภาษาไทย": """
กลไกระดับโมเลกุลของ DNA (ภาพรวมเชิงลึก)

การจัดบรรจุ DNA และโครมาติน  
DNA ของยูคาริโอตถูกบรรจุอยู่ในรูปของโครมาติน โดยมีหน่วยพื้นฐานคือ นิวคลีโอโซม
ซึ่งประกอบด้วย DNA ประมาณ 146–147 คู่เบสพันรอบโปรตีนฮิสโตนจำนวน 8 ตัว
โครมาตินแบ่งออกเป็นสองสภาพหลัก คือ ยูโครมาตินซึ่งมีโครงสร้างหลวมและเอื้อต่อการถอดรหัส
และเฮเทอโรโครมาตินซึ่งมีโครงสร้างแน่นและยับยั้งการแสดงออกของยีน

การปรับแต่งฮิสโตน  
หางของโปรตีนฮิสโตนสามารถถูกปรับแต่งทางเคมี เช่น การอะซิทิลและการเมทิล
ซึ่งส่งผลต่อระดับการเข้าถึง DNA และการแสดงออกของยีน
ตัวอย่างเช่น การอะซิทิลของฮิสโตนตำแหน่ง H3K27ac มักเกี่ยวข้องกับบริเวณที่มีการถอดรหัสสูง
ในขณะที่การเมทิลของฮิสโตนบางตำแหน่งอาจกระตุ้นหรือยับยั้งการถอดรหัส ขึ้นอยู่กับบริบทของตำแหน่งนั้น

การเมทิลของ DNA  
การเมทิลของ DNA มักเกิดขึ้นที่ไซโตซีนในตำแหน่ง CpG
โดยเฉพาะบริเวณ CpG islands ใกล้โปรโมเตอร์ของยีน
การเมทิลสามารถยับยั้งการจับของ transcription factors
และนำไปสู่การปิดการแสดงออกของยีน

การปรับโครงสร้างโครมาติน  
โปรตีนกลุ่ม chromatin remodeling complexes ใช้พลังงานจาก ATP
ในการเคลื่อนย้ายหรือปรับตำแหน่งของนิวคลีโอโซม
ทำให้ DNA เปิดหรือปิดต่อกระบวนการถอดรหัส

เอพิเจเนติกส์  
เอพิเจเนติกส์คือการเปลี่ยนแปลงการแสดงออกของยีน
โดยไม่เปลี่ยนลำดับเบสของ DNA
กลไกหลักได้แก่ DNA methylation และ histone modifications

DNA และ RNA ที่ไม่เข้ารหัส  
แม้ DNA ส่วนใหญ่ไม่เข้ารหัสโปรตีน
แต่ส่วนที่ไม่เข้ารหัส เช่น microRNA และ long non coding RNA
มีบทบาทสำคัญในการควบคุมการแสดงออกของยีน

---

การวิเคราะห์ดีเอ็นเอของมนุษย์ด้วยชีวสารสนเทศ

ความหมายของชีวสารสนเทศ  
ชีวสารสนเทศคือการประยุกต์ใช้คอมพิวเตอร์และอัลกอริทึม
เพื่อวิเคราะห์ข้อมูลลำดับ DNA ศึกษาความแตกต่างทางพันธุกรรม
วิวัฒนาการ และการกลายพันธุ์ที่เกี่ยวข้องกับโรค

---

วัตถุประสงค์ของการศึกษา  
1 ศึกษาดีเอ็นเอของมนุษย์เพื่อเข้าใจสายพันธุ์และความสัมพันธ์  
2 สร้างแผนภูมิต้นไม้วิวัฒนาการเพื่อศึกษาความสัมพันธ์ของดีเอ็นเอและหน้าที่ของยีน

ประโยชน์ที่คาดว่าจะได้รับ  
1 เข้าใจสายพันธุ์และความสัมพันธ์ของมนุษย์ผ่านดีเอ็นเอ  
2 เข้าใจต้นกำเนิดและพัฒนาการของมนุษย์และนำไปศึกษาต่อทางการแพทย์

---

กระบวนการวิเคราะห์ด้วยชีวสารสนเทศ

แหล่งที่มาของข้อมูล  
ข้อมูลลำดับดีเอ็นเอได้มาจากฐานข้อมูลสาธารณะ NCBI GenBank
โดยใช้หมายเลข accession ของแต่ละสายพันธุ์

การเตรียมข้อมูลและการควบคุมคุณภาพ  
1 ดึงลำดับยีน ND5 จากจีโนมไมโตคอนเดรีย  
2 ความยาวของยีน ND5 แตกต่างกันเล็กน้อยในแต่ละสายพันธุ์  
3 ความแตกต่างเกิดจากกระบวนการวิวัฒนาการและการกลายพันธุ์  
4 ความแตกต่างดังกล่าวเป็นเรื่องปกติทางชีววิทยา  
5 คัดกรองลำดับที่มีเบสไม่ทราบชนิดจำนวนมากหรือมีการหยุดการแปลผิดปกติ

การจัดเรียงหลายลำดับ  
ทำการจัดเรียงลำดับ ND5 เพื่อให้ตำแหน่งที่เปรียบเทียบกันสอดคล้องกัน

การคำนวณระยะทางและการสร้างต้นไม้วิวัฒนาการ  
1 คำนวณค่าความเหมือนและระยะทางเชิงวิวัฒนาการ  
2 สร้างแผนภูมิต้นไม้วิวัฒนาการจากข้อมูลยีน ND5

---

หน้าที่ทางชีววิทยาของยีน ND5  
ยีน MT ND5 เข้ารหัสหน่วยย่อยของ Complex I
ในระบบขนส่งอิเล็กตรอนของไมโตคอนเดรีย
ซึ่งมีบทบาทสำคัญในการสร้างพลังงานของเซลล์

ความสำคัญทางการแพทย์  
การกลายพันธุ์ในยีน ND5 เกี่ยวข้องกับโรคไมโตคอนเดรีย
เช่น MELAS และ Leigh syndrome

ความสามารถในการทำซ้ำและการศึกษาต่อยอด  
กระบวนการวิเคราะห์สามารถนำไปทำซ้ำ
และประยุกต์ใช้กับข้อมูลหรือยีนอื่น
เพื่อการศึกษาวิจัยในอนาคต
"""
}



# ========== SIDEBAR ==========
# ===============================

T = LANGS[language]

if language == "ภาษาไทย":
    species_display = {
        eng: SPECIES_TH.get(eng, eng)
        for eng in species_accessions
    }
else:
    species_display = {
        eng: eng
        for eng in species_accessions
    }


# DEBUG: show missing keys early

REQUIRED_KEYS = [
    "select_species", "fetch", "identity",
    "phylo_tree", "preview",
    "available", "selected"
]
for k in REQUIRED_KEYS:
    assert k in T, f"Missing LANG key: {k}"



# optional: small logo (uncomment and provide local path or URL)
# st.sidebar.image("logo.png", width=120)

st.title(T["title"])
st.caption(T["subtitle"])

# ========== UI LAYOUT: control panel ==========
with st.container():

    # --- two-column row ---
    controls_col, info_col = st.columns([2.2, 1])

    # ===== LEFT: species selector + fetch button =====
    with controls_col:
        selected_labels = st.multiselect(
            T["select_species"],
            options=list(species_display.values()),
            placeholder=T["choose_options"],
            key="species_multiselect"
        )

        # label → English key
        label_to_eng = {v: k for k, v in species_display.items()}
        selected_species = [label_to_eng[label] for label in selected_labels]

        # ✅ FETCH BUTTON (CORRECT LOCATION)
        fetch_clicked = st.button(
            T["fetch"],
        )

    # ===== RIGHT: info panel =====
    with info_col:
        st.metric(T["available"], len(species_accessions))
        st.write("")
        st.write(f"{T['selected']} {len(selected_species)}")


# ========== CACHING FOR FETCH ==========
@st.cache_data(show_spinner=False)
def extract_nd5_cached(accession):
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()

        for feature in record.features:
            if feature.type == "CDS":
                gene    = feature.qualifiers.get("gene",    [""])[0].lower().strip()
                product = feature.qualifiers.get("product", [""])[0].lower().strip()
                note    = feature.qualifiers.get("note",    [""])[0].lower().strip()

                is_nd5 = any([
                    gene in ["nd5", "nad5", "mt-nd5", "mtnd5"],
                    "nad5"  in product,
                    "nd5"   in product,
                    "ndh5"  in product,
                    "subunit 5" in product and "nadh" in product,
                    "nad5"  in note,
                    "nd5"   in note,
                ])

                if is_nd5:
                    return feature.extract(record.seq)

        # ── fallback: try gene features too (some records skip CDS) ──
        for feature in record.features:
            if feature.type == "gene":
                gene = feature.qualifiers.get("gene", [""])[0].lower().strip()
                if gene in ["nd5", "nad5", "mt-nd5", "mtnd5"]:
                    return feature.extract(record.seq)

        return None

    except Exception:
        return None



# ========== HELPERS ==========
def percent_identity(seq1, seq2):
    matches, valid = 0, 0
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        valid += 1
        if a == b:
            matches += 1
    return round(matches / valid * 100, 2) if valid > 0 else 0.0

def build_fasta_text(seq_dict):
    s = ""
    for name, seq in seq_dict.items():
        s += f">{name}\n{seq}\n"
    return s

def clean_label(x):
    if x.name:
        return re.sub(r"\s*\([^)]*\)", "", x.name)
    return ""

def translate_species(name, language):
    if language == "ภาษาไทย":
        return SPECIES_TH.get(name, name)
    return name


# ========== PROCESS: FETCH ND5 ==========
if fetch_clicked:
    if len(selected_species) < 2:
        st.warning(T["select_species"])
    else:
        nd5_seqs = {}
        missing_species = []   # 👈 เพิ่ม

        for eng_name in selected_species:
            acc = species_accessions[eng_name]
            seq = extract_nd5_cached(acc)

            if seq is None:
                missing_species.append(species_display.get(eng_name, eng_name))

                continue

            display_name = species_display.get(eng_name, eng_name)
            nd5_seqs[display_name] = seq

        nd5_lengths = {k: len(v) for k, v in nd5_seqs.items()}
                    # ===== FORCE SAME ND5 LENGTH =====
        
        if nd5_seqs:
            max_len = max(len(seq) for seq in nd5_seqs.values())
            for k in nd5_seqs:
                nd5_seqs[k] = nd5_seqs[k] + "-" * (max_len - len(nd5_seqs[k]))

        st.session_state["nd5_seqs"] = nd5_seqs
        st.session_state["nd5_lengths"] = nd5_lengths
        st.session_state["valid_species"] = list(nd5_seqs.keys())

        st.session_state["nd5_seqs"] = nd5_seqs
        if missing_species:
            st.warning(
                f"{T['missing_nd5']}\n\n"
                + ", ".join(missing_species)
            )
        # AUTO-COMPUTE IDENTITY MATRIX
        names = list(nd5_seqs.keys())
        matrix = []
        for i, a in enumerate(names):
            row = []
            for j, b in enumerate(names):
                row.append(100.0 if i == j else percent_identity(nd5_seqs[a], nd5_seqs[b]))
            matrix.append(row)

        st.session_state["identity_df"] = pd.DataFrame(
            matrix, index=names, columns=names
        )

        # AUTO-BUILD TREE DATA
        fasta_io = StringIO()
        for name, seq in nd5_seqs.items():
            fasta_io.write(f">{name}\n{seq}\n")
        fasta_io.seek(0)

        records = [
            SeqRecord(Seq(str(seq)), id=name, description="")
            for name, seq in nd5_seqs.items()
        ]

        alignment = MultipleSeqAlignment(records)

        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        tree.ladderize()

        for clade in tree.find_clades():
            if clade.name and clade.name.startswith("Inner"):
                clade.name = None

        st.session_state["tree"] = tree


# ========== SHOW ND5 LENGTHS (compact) ==========
if "nd5_seqs" in st.session_state:
    seqs = st.session_state["nd5_seqs"]

    st.subheader(T["length"])

    cols = st.columns(4)
    for i, (name, seq) in enumerate(seqs.items()):
        real_len = st.session_state["nd5_lengths"][name]
        cols[i % 4].metric(name, f"{real_len} bp")

    # PREVIEW DIRECTLY UNDER
    with st.expander(T["preview"]):
        for name, seq in seqs.items():
            st.markdown(f"**{name}**")
            st.code(seq[:120] + ("..." if len(seq) > 120 else ""))

# ========== TABS FOR RESULTS / TREE / METHODS ==========

tab_results, tab_tree, tab_methods = st.tabs([
    T["results"],
    T["tree"],
    T["methods"]
])

# ===============================
# RESULTS TAB
# ===============================
with tab_results:
    st.subheader(T["results"])

    df = st.session_state.get("identity_df")

    if df is not None:

        # ===== translate row + column labels =====
        df_display = df.copy()

        df_display.index = [
            translate_species(i, language) for i in df_display.index
        ]
        df_display.columns = [
            translate_species(c, language) for c in df_display.columns
        ]

        # ===== insert translated Species header =====
        df_display.insert(0, T["species"], df_display.index)

        # ===== FIX: numeric text color by theme =====
        number_text_color = THEME["TABLE_TEXT"]

        styled_df = (
            df_display.style

            # ----- number format -----
            .format("{:.2f}", subset=df_display.columns[1:])

            # ----- gradient (numbers only) -----
            .background_gradient(
                cmap="Blues",
                subset=df_display.columns[1:],
                text_color_threshold=0   # 🔥 important fix
            )

            # ----- Species column -----
            .set_properties(
                subset=[T["species"]],
                **{
                    "text-align": "left",
                    "font-weight": "600",
                    "background-color": THEME["SIDEBAR_BG"],
                    "color": THEME["TEXT"],
                    "white-space": "nowrap",
                    "min-width": "220px",
                }
            )

            # ----- headers & borders -----
            .set_table_styles([
                {
                    "selector": "th",
                    "props": [
                        ("background-color", THEME["PANEL_BG"]),
                        ("color", THEME["TEXT"]),
                        ("border", f"1px solid {THEME['BORDER']}"),
                        ("font-weight", "700"),
                        ("text-align", "center"),
                        ("white-space", "nowrap"),
                        ("padding", "6px 10px"),
                    ],
                },
                {
                    "selector": "td",
                    "props": [
                        ("border", f"1px solid {THEME['BORDER']}"),
                        ("padding", "6px 10px"),
                    ],
                },
                {
                    "selector": "table",
                    "props": [
                        ("border-collapse", "collapse"),
                        ("border", f"1px solid {THEME['BORDER']}"),
                    ],
                },
            ])
        )

        st.dataframe(
            styled_df,
            width="stretch",
            hide_index=True
        )
        
        # ===== PAIRWISE DETAIL VIEW =====
        if df is not None and "nd5_seqs" in st.session_state:
            st.markdown("---")

            # ===== INSTRUCTION =====
            st.caption(
                "Select two species below to compare their ND5 identity score."
                if language == "English"
                else "เลือกสปีชีส์สองชนิดเพื่อเปรียบเทียบค่าความเหมือนของ ND5"
            )

            all_names = list(st.session_state["nd5_seqs"].keys())

            col_a, col_b = st.columns(2)
            with col_a:
                species_a = st.selectbox(
                    "Species A" if language == "English" else "สายพันธุ์ที่ 1",
                    options=all_names,
                    key="pair_a"
                )
            with col_b:
                species_b = st.selectbox(
                    "Species B" if language == "English" else "สายพันธุ์ที่ 2",
                    options=all_names,
                    index=min(1, len(all_names) - 1),
                    key="pair_b"
                )

            with st.expander(
                "🔍 Pairwise Detail" if language == "English" else "🔍 รายละเอียดคู่สปีชีส์",
                expanded=True
            ):
                if species_a == species_b:
                    st.info("Select two different species." if language == "English" else "กรุณาเลือกสปีชีส์ที่แตกต่างกัน")
                else:
                    # pull identity score from existing df
                    # df index/columns are in display names already
                    score = st.session_state["identity_df"].loc[species_a, species_b]

                    if language == "English":
                        delta_label = (
                            "identical" if score == 100 else
                            "very close" if score >= 95 else
                            "moderately similar" if score >= 80 else
                            "distantly related"
                        )
                    else:
                        delta_label = (
                            "เหมือนกันทุกประการ" if score == 100 else
                            "ใกล้ชิดมาก" if score >= 95 else
                            "คล้ายคลึงปานกลาง" if score >= 80 else
                            "ห่างไกลทางวิวัฒนาการ"
                        )

                    st.metric(
                        label=f"{species_a}  ↔  {species_b}",
                        value=f"{score:.2f}%",
                        delta=delta_label
                    )
    else:
        st.info(T["no_results"])


# ===============================
# METHODS TAB
# ===============================
with tab_methods:
    st.subheader(T["methods"])

    with st.expander(T["methods_expander"], expanded=False):
        st.markdown(METHODS_TEXT[language])

# ===============================
# TREE TAB
# ===============================
with tab_tree:
    st.subheader(T["phylo_tree"])

    tree = deepcopy(st.session_state.get("tree"))

    if tree is None:
        st.info(T["no_tree"])
    else:
        # ===== NORMALIZE BRANCH LENGTH =====
        lengths = [
            clade.branch_length
            for clade in tree.find_clades()
            if clade.branch_length
        ]

        if lengths:
            max_len = max(lengths)
            for clade in tree.find_clades():
                if clade.branch_length:
                    clade.branch_length /= max_len

        # ===== FIGURE (RESPONSIVE SAFE) =====
        fig = plt.figure(
            figsize=(10, 5),
            facecolor=THEME["TREE_BORDER_BG"]
        )
        ax = fig.add_subplot(111)
        ax.set_facecolor(THEME["PANEL_BG"])

        Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            show_confidence=False,
            label_func=lambda x: x.name if x.name else ""
        )

        TEXT_COLOR = "#9CA3AF"
        BRANCH_COLOR = "#CBD5E1"

        for text in ax.texts:
            text.set_fontproperties(font_prop)
            text.set_color(TEXT_COLOR)

        n = len(tree.get_terminals())
        font_size = max(8, 14 - n)

        thai_font = None
        if os.path.exists(FONT_PATH):
            thai_font = font_manager.FontProperties(fname=FONT_PATH)

        for text in ax.texts:
            text.set_fontsize(font_size)
            text.set_color(TEXT_COLOR)
            if thai_font:
                text.set_fontproperties(thai_font)
                text.set_fontfamily("sans-serif")

        for collection in ax.collections:
            if isinstance(collection, LineCollection):
                collection.set_color(BRANCH_COLOR)
                collection.set_linewidth(2.2)
                collection.set_alpha(1.0)

        ax.margins(x=0.06, y=0.22)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        for spine in ax.spines.values():
            spine.set_visible(False)

        fig.tight_layout()
        st.pyplot(fig, width = "stretch")

        buf = BytesIO()
        fig.savefig(
            buf,
            format="png",
            dpi=300,
            bbox_inches="tight",
            facecolor=fig.get_facecolor()
        )

        _, c, _ = st.columns([1, 2, 1])
        with c:
            st.download_button(
                "Download phylogenetic tree (PNG)",
                data=buf.getvalue(),
                file_name="nd5_phylogenetic_tree.png",
                mime="image/png"
            )