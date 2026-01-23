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

font_prop = font_manager.FontProperties(fname=FONT_PATH)

# ตั้งค่าทั้งระบบ
rcParams["font.family"] = font_prop.get_name()
rcParams["axes.unicode_minus"] = False

# ========== PAGE CONFIG ==========
st.set_page_config(
    page_title="ND5 Phylogenetic Explorer",
    page_icon="🧬",
    layout="centered",
    initial_sidebar_state="expanded"
)
# ===============================
# THEME STATE
# ===============================
if "theme_mode" not in st.session_state:
    st.session_state.theme_mode = "Dark"

with st.container():
    st.markdown("### ⚙️ Settings")
    st.radio("🎨 Theme", ["Dark", "Light"], key="theme_mode")
    language = st.selectbox("🌐 Language", ["English", "ภาษาไทย"])


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
    "TREE_BORDER_BG" : "#9ca3af"
}

LIGHT = {
    "BG": "#ffffff",
    "TEXT": "#111827",
    "TEXT_MUTED": "#6b7280",
    "SIDEBAR_BG": "#f3f4f6",
    "PANEL_BG": "#ffffff",
    "BORDER": "#d1d5db",
    "CODE_BG": "#f9fafb",
    "CODE_TEXT": "#111827",
    "BTN_BG": "#2563eb",
    "BTN_TEXT": "#ffffff",
    "TREE_BORDER_BG" : "#d8e6fe"
}
THEME = DARK if st.session_state.theme_mode == "Dark" else LIGHT
# ---------- THEME STATE (SAFE DEFAULT) ----------
st.markdown(
f"""
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
}}

html, body, div, span, p, label,
h1, h2, h3, h4, h5, h6 {{
    color: var(--text) !important;
}}

/* ================== SIDEBAR ================== */
section[data-testid="stSidebar"] {{
    background-color: var(--sidebar-bg) !important;
    border-right: 1px solid var(--border) !important;
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

.stDownloadButton svg {{
    fill: var(--btn-text) !important;
}}

/* ================== EXPANDERS ================== */
div[data-testid="stExpander"] summary {{
    background-color: var(--panel-bg) !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
    padding: 12px !important;
}}

div[data-testid="stExpander"] div[role="region"] {{
    background-color: var(--panel-bg) !important;
    border: 1px solid var(--border) !important;
    border-top: none !important;
    border-radius: 0 0 12px 12px !important;
    padding: 16px !important;
}}

/* ================== PANELS ================== */
div[data-testid="stMetric"],
div[data-testid="stDataFrame"] {{
    background-color: var(--panel-bg) !important;
    border: 1px solid var(--border) !important;
    border-radius: 12px !important;
}}

/* ================== DATAFRAME FIX (IMPORTANT) ================== */

/* full width */
[data-testid="stDataFrame"] {{
    width: 100% !important;
}}

/* header wrap */
[data-testid="stDataFrame"] th {{
    white-space: normal !important;
    word-break: break-word !important;
    text-align: center !important;
    font-size: 14px !important;
    line-height: 1.3 !important;
}}

/* cell wrap */
[data-testid="stDataFrame"] td {{
    white-space: normal !important;
    word-break: break-word !important;
    text-align: center !important;
    font-size: 14px !important;
}}

/* index column */
[data-testid="stDataFrame"] th:first-child {{
    min-width: 140px !important;
}}

/* ================== CODE BLOCKS ================== */
pre, code {{
    background-color: var(--code-bg) !important;
    color: var(--code-text) !important;
    border-radius: 10px !important;
}}

/* ================== SELECT ================== */
div[data-baseweb="select"] > div {{
    background-color: var(--panel-bg) !important;
    border: 1px solid var(--border) !important;
}}

div[data-baseweb="select"] span,
div[data-baseweb="select"] input {{
    color: var(--text) !important;
}}

li[role="option"] {{
    background-color: var(--panel-bg) !important;
    color: var(--text) !important;
}}

/* ================== METRICS ================== */
div[data-testid="stMetricLabel"],
div[data-testid="stMetricValue"] {{
    color: var(--text) !important;
}}

/* ================== DIVIDER ================== */
hr {{
    border-top: 1px solid var(--border) !important;
}}

/* ================== SVG (PHYLO TREE) ================== */
svg text,
svg tspan {{
    fill: var(--text) !important;
}}

/* ================== DESKTOP ONLY ================== */
@media (min-width: 768px) {{

    /* ซ่อน toolbar อย่างเดียว */
    div[data-testid="stToolbar"],
    div[data-testid="stElementToolbar"] {{
        display: none !important;
    }}

    /* ❌ อย่าซ่อน header */
    header {{
        visibility: visible !important;
    }}
/* DESKTOP ONLY */
@media (min-width: 768px) {{

    div[data-testid="stToolbar"],
    div[data-testid="stElementToolbar"] {{
        display: none !important;
    }}

    header {{
        visibility: visible !important;
    }}
}}



</style>
""",
unsafe_allow_html=True
)




# ========== NCBI CONFIG ==========
Entrez.email = "poonthakorn@gmail.com"

# ========== SPECIES ACCESSIONS ==========
species_accessions = {
    # === Hominidae ===
    "Human": "NC_012920.1",
    "Chimpanzee": "NC_001643.1",
    "Bonobo": "NC_001644.1",
    "Eastern Gorilla": "NC_011120.1",
    "Bornean Orangutan": "NC_001646.1",

    # === Lesser apes ===
    "Common Gibbon": "NC_002082.1",

    # === Old World monkeys ===
    "Rhesus Macaque": "NC_005943.1",
    "Crab-eating Macaque": "NC_012670.1",
    "Barbary Ape": "NC_002764.1",
    "Hamadrayas Baboon": "NC_001992.1",
    "Green Monkey": "NC_008066.1",
    "Tantalus Monkey": "NC_009748.1",
    "Dusky Leaf Monkey": "NC_006900.1",
    "Proboscis Monkey": "NC_008216.1",
    "Red-shanked Douc": "NC_008220.1",
    "Tonkin Snub-nose Monkey": "NC_015485.1",
    "Colobus Monkey": "NC_006901.1",

    # === Other mammals ===
    "Pig": "NC_000845.1",
    "Cow": "NC_001567.1",
    "Horse": "NC_001640.1",
    "Dog": "NC_002008.4",
    "Cat": "NC_001700.1",
    "Sheep": "NC_001941.1",
    "Goat": "NC_005044.2",

    # === Rodents ===
    "Mouse": "NC_005089.1",
    "Rat": "NC_001665.2",

    # === Added to reach 30 ===
    "Rabbit": "NC_001913.1",
    "Guinea Pig": "NC_000884.1",
    "Elephant": "NC_000934.1",
    "Whale": "NC_006853.1"
}


SPECIES_TH = {
    "Human": "มนุษย์",
    "Chimpanzee": "ชิมแปนซี",
    "Bonobo": "โบโนโบ",
    "Eastern Gorilla": "กอริลลาตะวันออก",
    "Bornean Orangutan": "อุรังอุตังบอร์เนียว",
    "Common Gibbon": "ชะนี",
    "Rhesus Macaque": "ลิงแสม",
    "Crab-eating Macaque": "ลิงแสมหางยาว",
    "Barbary Ape": "ลิงบาบารี",
    "Hamadrayas Baboon": "ลิงบาบูนฮามาดรายัส",
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
    "Guinea Pig": "หนูตะเภา",
    "Elephant": "ช้าง",
    "Whale": "วาฬ"
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
        "selected": "ที่เลือก",
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
        "methods_expander": "🔬 กลไกระดับโมเลกุล ขั้นตอน และวัตถุประสงค์"
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
        "methods_expander": "🔬 Molecular mechanisms, pipeline and objectives"
    }
}

# ===============================
# METHODS TAB
# ===============================
METHODS_TEXT = {
    "English": """
🔬 Molecular Mechanisms of DNA (Comprehensive Overview)

DNA Packaging and Chromatin Structure  
Eukaryotic DNA is packaged into chromatin, whose fundamental unit is the **nucleosome**.
Each nucleosome consists of approximately 146–147 base pairs of DNA wrapped around an octamer of histone proteins.
Chromatin exists in two major forms: **euchromatin**, which is loosely packed and transcriptionally active,
and **heterochromatin**, which is densely packed and transcriptionally repressive.

Histone Modifications  
Histone tails undergo post-translational modifications such as **acetylation** and **methylation**,
which influence chromatin accessibility and gene expression.
For example, histone acetylation (e.g., **H3K27ac**) is generally associated with transcriptionally active regions,
while histone methylation can either activate or repress transcription depending on the residue and context
(e.g., H3K4, H3K9, H3K27).

DNA Methylation (CpG Islands)  
DNA methylation commonly occurs at cytosine residues in CpG dinucleotides.
Hypermethylation of CpG islands in promoter regions can inhibit transcription factor binding
and lead to gene silencing, playing an important role in epigenetic regulation.

Chromatin Remodeling  
ATP-dependent chromatin remodeling complexes reposition or remove nucleosomes,
thereby regulating the accessibility of DNA to transcriptional machinery.

Transcription and Cis-Regulatory Elements  
Transcription initiates at promoter regions with the binding of RNA polymerase II
and transcription factors. **Enhancers** increase transcriptional activity,
whereas **silencers** suppress gene expression.

Epigenetics  
Epigenetics refers to heritable changes in gene expression that do not involve alterations
in the DNA sequence. DNA methylation and histone modifications are major epigenetic mechanisms
that define cell-type-specific gene expression patterns.

Non-coding DNA and RNA  
Although most of the human genome does not encode proteins,
non-coding elements such as **microRNAs (miRNAs)** and **long non-coding RNAs (lncRNAs)**
play crucial roles in post-transcriptional gene regulation by inhibiting translation
or promoting mRNA degradation.

---

Human DNA Analysis Using Bioinformatics

### Definition of Bioinformatics  
**Bioinformatics** is the application of computational tools and algorithms
to analyze biological sequence data, enabling the study of genetic variation,
gene function, evolution, and disease-associated mutations.

---

Objectives of the Study
1. To study human DNA in order to understand human genetic variation and species relationships.  
2. To construct **phylogenetic trees** to investigate evolutionary relationships
between humans and other species based on DNA sequence differences,
and to explore functional aspects of DNA in animals.  
3. To analyze sequence variation in the mitochondrial **ND5 gene**
and assess its evolutionary and biological significance.

Expected Benefits
1. To enhance understanding of human genetic diversity and evolutionary relationships
through DNA sequence analysis.  
2. To provide insights into the origin and evolutionary divergence of humans
and to establish a foundation for further studies in medical and genetic research,
particularly in mitochondrial-related diseases.

---

Bioinformatics Analysis Pipeline

Data Source  
DNA sequences used in this study are retrieved from the public database
**NCBI GenBank** using accession numbers for each species.
All sequences are publicly available and do not contain personal or sensitive information.

Data Preparation and Quality Control  
- Extraction of the **ND5 (MT-ND5)** coding sequence from mitochondrial genomes.  
- All sequences are standardized to a common length after extraction, 
    based on the longest ND5 sequence, to ensure accurate alignment and comparison.  
- Sequences with insufficient length, excessive ambiguous bases (N),
or internal stop codons are excluded from analysis.
- The length of the ND5 gene varies slightly among species (e.g., 1,812 bp, 1,815 bp, or 1,821 bp).
- This variation is caused by evolutionary insertions and deletions (indels), differences in start or stop codon annotation
and species-specific mutations accumulated over time.

Such length variation is biologically normal and reflects evolutionary divergence rather than sequencing error.

Multiple Sequence Alignment  
Multiple sequence alignment (MSA) is performed to align ND5 sequences
from different species, ensuring homologous positions are compared.

Evolutionary Distance and Phylogenetic Tree Construction  
- Pairwise sequence similarity (% identity) and evolutionary distances are calculated.  
- A **phylogenetic tree** is constructed to visualize evolutionary relationships
among species based on ND5 sequence variation.

---

Biological Function of the ND5 Gene
The mitochondrial **MT-ND5** gene encodes a subunit of **Complex I**
(NADH dehydrogenase) in the mitochondrial electron transport chain.
This complex plays a critical role in oxidative phosphorylation
and ATP production.

Why the ND5 Gene Is Used in Evolutionary Studies
The mitochondrial ND5 gene is commonly used in phylogenetic and evolutionary studies for several reasons.
-First, mitochondrial DNA is maternally inherited and does not undergo recombination, allowing evolutionary 
changes to accumulate in a relatively predictable manner.

-Second, ND5 evolves at a moderate rate, making it suitable for resolving evolutionary relationships 
among closely related species such as mammals and primates.

-Third, ND5 is a protein-coding gene involved in cellular energy production, meaning its sequence is evolutionarily conserved while still allowing informative mutations to occur.

These characteristics make ND5 particularly useful for studying species divergence, 
evolutionary distance, and phylogenetic relationships.

---

Medical Significance
Mutations in the ND5 gene are associated with mitochondrial disorders
such as **MELAS** and **Leigh syndrome**, which impair cellular energy production.
Therefore, ND5 is an important gene for both evolutionary studies
and medical genetics.

---

Reproducibility and Future Applications
The analytical workflow and computational steps can be reproduced
and extended to larger datasets. This approach provides a foundation
for future research in evolutionary biology, population genetics,
and biomedical studies.
""",           
    "ภาษาไทย": """
## 🔬 กลไกระดับโมเลกุลของ DNA (ภาพรวมเชิงลึก)

การจัดบรรจุ DNA และโครมาติน  
DNA ของยูคาริโอตถูกบรรจุอยู่ในรูปของโครมาติน โดยมีหน่วยพื้นฐานคือ **นิวคลีโอโซม**
ซึ่งประกอบด้วย DNA ประมาณ 146–147 คู่เบสพันรอบโปรตีนฮิสโตนจำนวน 8 ตัว
โครมาตินสามารถแบ่งออกเป็น 2 สภาพหลัก ได้แก่  
**ยูโครมาติน (euchromatin)** ซึ่งมีโครงสร้างหลวมและเอื้อต่อการถอดรหัส
และ **เฮเทอโรโครมาติน (heterochromatin)** ซึ่งมีโครงสร้างแน่นและยับยั้งการแสดงออกของยีน

การปรับแต่งฮิสโตน (Histone Modifications)  
หางของโปรตีนฮิสโตนสามารถถูกปรับแต่งทางเคมี เช่น **การอะซิทิล (acetylation)**
และ **การเมทิล (methylation)** ซึ่งส่งผลต่อระดับการเข้าถึง DNA และการแสดงออกของยีน
ตัวอย่างเช่น การอะซิทิลของฮิสโตนตำแหน่ง **H3K27ac** มักเกี่ยวข้องกับบริเวณที่มีการถอดรหัสสูง
ในขณะที่การเมทิลของฮิสโตนบางตำแหน่งอาจกระตุ้นหรือยับยั้งการถอดรหัส ขึ้นอยู่กับบริบทของตำแหน่งนั้น

การเมทิลของ DNA (CpG Islands)  
การเมทิลของ DNA มักเกิดขึ้นที่ไซโตซีนในตำแหน่ง CpG dinucleotides
โดยเฉพาะบริเวณ CpG islands ที่อยู่ใกล้โปรโมเตอร์ของยีน
การเมทิลในบริเวณดังกล่าวสามารถยับยั้งการจับของ transcription factors
และนำไปสู่การปิดการแสดงออกของยีน ซึ่งเป็นกลไกสำคัญของการควบคุมยีนในระดับเอพิเจเนติกส์

การปรับโครงสร้างโครมาติน (Chromatin Remodeling)  
โปรตีนกลุ่ม chromatin remodeling complexes ใช้พลังงานจาก ATP
ในการเคลื่อนย้ายหรือปรับตำแหน่งของนิวคลีโอโซม
ทำให้ DNA เปิดหรือปิดต่อเครื่องจักรถอดรหัสได้ตามความเหมาะสม

การถอดรหัสและองค์ประกอบควบคุมยีน (Cis-Regulatory Elements)  
กระบวนการถอดรหัสเริ่มต้นที่บริเวณโปรโมเตอร์ โดย RNA polymerase II
และ transcription factors  
**Enhancer** มีบทบาทในการเพิ่มระดับการถอดรหัส
ขณะที่ **silencer** ทำหน้าที่ยับยั้งการแสดงออกของยีน

เอพิเจเนติกส์ (Epigenetics)  
เอพิเจเนติกส์หมายถึงการเปลี่ยนแปลงการแสดงออกของยีน
โดยไม่เปลี่ยนแปลงลำดับเบสของ DNA
กลไกหลักได้แก่ DNA methylation และ histone modifications
ซึ่งช่วยกำหนดรูปแบบการแสดงออกของยีนที่แตกต่างกันในแต่ละชนิดของเซลล์

DNA และ RNA ที่ไม่เข้ารหัส  
แม้ว่า DNA ส่วนใหญ่ของจีโนมมนุษย์จะไม่เข้ารหัสโปรตีน
แต่ส่วนที่ไม่เข้ารหัสเหล่านี้ เช่น **microRNA (miRNA)** และ **long non-coding RNA (lncRNA)**
มีบทบาทสำคัญในการควบคุมการแสดงออกของยีน
โดยการยับยั้งการแปลโปรตีนหรือกระตุ้นการสลายของ mRNA

---

🧬 การวิเคราะห์ดีเอ็นเอของมนุษย์ด้วยชีวสารสนเทศ

ความหมายของชีวสารสนเทศ  
**ชีวสารสนเทศ (Bioinformatics)** คือการประยุกต์ใช้คอมพิวเตอร์
และอัลกอริทึมในการวิเคราะห์ข้อมูลชีวภาพ โดยเฉพาะข้อมูลลำดับ DNA
เพื่อศึกษาโครงสร้าง ความแตกต่างทางพันธุกรรม วิวัฒนาการ
และการกลายพันธุ์ที่เกี่ยวข้องกับโรค

---

🎯 วัตถุประสงค์ของการศึกษา
1. เพื่อศึกษาดีเอ็นเอของมนุษย์และทำความเข้าใจความหลากหลายทางพันธุกรรม
รวมถึงความสัมพันธ์ระหว่างมนุษย์กับสายพันธุ์อื่น  
2. เพื่อสร้าง **แผนภูมิวิวัฒนาการ (phylogenetic tree)**
สำหรับอธิบายความสัมพันธ์เชิงวิวัฒนาการจากความแตกต่างของลำดับดีเอ็นเอ
และศึกษาหน้าที่ของยีนในสัตว์ชนิดต่าง ๆ  
3. เพื่อวิเคราะห์ความแปรผันของลำดับยีนไมโตคอนเดรีย **ND5**
และประเมินความสำคัญทางชีววิทยาและวิวัฒนาการของยีนดังกล่าว

---

✅ ประโยชน์ที่คาดว่าจะได้รับ
1. ทำให้เข้าใจความสัมพันธ์และสายวิวัฒนาการของมนุษย์
ผ่านการวิเคราะห์ข้อมูลลำดับดีเอ็นเอ  
2. ทำให้ทราบถึงต้นกำเนิดและการพัฒนาการของมนุษย์ในเชิงวิวัฒนาการ
และสามารถนำความรู้ไปต่อยอดในการศึกษาทางการแพทย์
โดยเฉพาะโรคที่เกี่ยวข้องกับไมโตคอนเดรีย

---

🛠️ กระบวนการวิเคราะห์ด้วยชีวสารสนเทศ

แหล่งที่มาของข้อมูล  
ลำดับดีเอ็นเอที่ใช้ในการศึกษานี้ถูกดึงมาจากฐานข้อมูลสาธารณะ
**NCBI GenBank** โดยใช้หมายเลข accession ของแต่ละสายพันธุ์
ข้อมูลทั้งหมดเป็นข้อมูลที่เปิดเผยเพื่อการวิจัย

การเตรียมข้อมูลและการควบคุมคุณภาพ  
- ดึงลำดับยีน **ND5 (MT-ND5)** จากจีโนมไมโตคอนเดรีย  
- ความยาวของยีน ND5 แตกต่างกันเล็กน้อยในแต่ละสปีชีส์ (เช่น 1,812 bp, 1,815 bp หรือ 1,821 bp)
ความแตกต่างนี้เกิดจากการแทรกและการลบในระหว่างวิวัฒนาการ ความแตกต่างในการกำหนดรหัสเริ่มต้นหรือรหัสหยุด 
และการกลายพันธุ์เฉพาะสปีชีส์ที่สะสมมาตามกาลเวลา
ความแตกต่างของความยาวดังกล่าวเป็นเรื่องปกติทางชีววิทยาและสะท้อนถึงความแตกต่างทางวิวัฒนาการมากกว่าข้อผิดพลาดในการจัดลำดับดีเอ็นเอ
- คัดกรองลำดับที่มีเบสไม่ทราบชนิด (N) จำนวนมาก
หรือมีสัญญาณการหยุดการแปลภายในกรอบการอ่าน

การจัดเรียงหลายลำดับ (Multiple Sequence Alignment)  
ทำการจัดเรียงลำดับ ND5 ของแต่ละสายพันธุ์
เพื่อให้ตำแหน่งที่เปรียบเทียบกันเป็นตำแหน่งที่สอดคล้องกันทางชีววิทยา

การคำนวณระยะทางและการสร้างต้นไม้วิวัฒนาการ  
- คำนวณค่าความเหมือนของลำดับเบสและระยะทางทางวิวัฒนาการ  
- สร้าง **แผนภูมิต้นไม้วิวัฒนาการ**
เพื่อแสดงความสัมพันธ์ของสายพันธุ์ต่าง ๆ จากข้อมูลยีน ND5



---

หน้าที่ทางชีววิทยาของยีน ND5
ยีน **MT-ND5** ทำหน้าที่เข้ารหัสหน่วยย่อยของ **Complex I**
ในระบบขนส่งอิเล็กตรอนของไมโตคอนเดรีย
ซึ่งมีบทบาทสำคัญในกระบวนการสร้างพลังงาน (ATP)
ของเซลล์

เหตุผลที่เลือกใช้ยีน ND5 ในการศึกษาวิวัฒนาการ

ยีนไมโตคอนเดรีย ND5 ถูกนำมาใช้อย่างแพร่หลายในการศึกษาทางวิวัฒนาการและสายสัมพันธ์ของสิ่งมีชีวิต
เนื่องจากมีคุณสมบัติสำคัญดังต่อไปนี้
    ประการแรก ดีเอ็นเอของไมโตคอนเดรียถ่ายทอดทางมารดาและไม่มีการเกิด recombination
ทำให้การเปลี่ยนแปลงของลำดับเบสสะสมไปตามกาลเวลาอย่างเป็นลำดับ
เหมาะสำหรับการศึกษาความสัมพันธ์เชิงวิวัฒนาการ
    ประการที่สอง ยีน ND5 มีอัตราการกลายพันธุ์ในระดับปานกลาง
สามารถแยกความแตกต่างระหว่างสายพันธุ์ที่มีความใกล้ชิดกัน
เช่น กลุ่มสัตว์เลี้ยงลูกด้วยนมและไพรเมต ได้อย่างมีประสิทธิภาพ
    ประการที่สาม ND5 เป็นยีนที่เข้ารหัสโปรตีนซึ่งมีบทบาทสำคัญในการสร้างพลังงานของเซลล์
ทำให้ลำดับยีนมีความคงตัวในระดับหนึ่ง
แต่ยังคงมีตำแหน่งที่เปลี่ยนแปลงได้เพียงพอสำหรับใช้เป็นข้อมูลเชิงวิวัฒนาการ

ด้วยเหตุผลเหล่านี้ ยีน ND5 จึงเหมาะสมอย่างยิ่งสำหรับการศึกษา
ความแตกต่างระหว่างสายพันธุ์ ระยะทางเชิงวิวัฒนาการ
และการสร้างแผนภูมิต้นไม้วิวัฒนาการ (phylogenetic tree)

---

🚨 ความสำคัญทางการแพทย์
การกลายพันธุ์ในยีน ND5 มีความเกี่ยวข้องกับโรคไมโตคอนเดรีย
เช่น **MELAS** และ **Leigh syndrome**
ซึ่งส่งผลต่อการสร้างพลังงานของเซลล์
ดังนั้นยีน ND5 จึงมีความสำคัญทั้งในด้านชีววิทยาเชิงวิวัฒนาการ
และด้านพันธุศาสตร์ทางการแพทย์

---

🔁 ความสามารถในการทำซ้ำและการศึกษาต่อยอด
กระบวนการวิเคราะห์ที่ใช้ในการศึกษานี้สามารถนำไปทำซ้ำ
และประยุกต์ใช้กับข้อมูลจากสายพันธุ์อื่นหรือยีนอื่นได้
ซึ่งเป็นพื้นฐานสำหรับการวิจัยในอนาคต
ด้านชีววิทยาวิวัฒนาการ พันธุศาสตร์ประชากร
และชีวการแพทย์
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
            use_container_width=True
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

        found = False   # 👈 เพิ่มบรรทัดนี้

        for feature in record.features:
            if feature.type == "gene":
                gene = feature.qualifiers.get("gene", [""])[0].lower()
                if gene in ["nd5", "nad5"]:
                    found = True
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    return record.seq[start:end]

        # ❌ ไม่มี ND5 จริง
        if not found:
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
                    # ===== FORCE SAME ND5 LENGTH =====
        if nd5_seqs:
            max_len = max(len(seq) for seq in nd5_seqs.values())
            for k in nd5_seqs:
                nd5_seqs[k] = nd5_seqs[k] + "-" * (max_len - len(nd5_seqs[k]))

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
        cols[i % 4].metric(name, f"{len(seq)} bp")

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
            figsize=(10, 5),  # ขนาดกลาง ใช้ได้ทั้ง mobile + desktop
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
        # ===== COLORS =====
        
        TEXT_COLOR = "#9CA3AF" if st.session_state.theme_mode == "Dark" else "#000000"
        BRANCH_COLOR = "#CBD5E1" if st.session_state.theme_mode == "Dark" else "#000000"
        for text in ax.texts:
            text.set_fontproperties(font_prop)
            text.set_color(TEXT_COLOR)
        # ===== AUTO SCALE TEXT =====
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

        # ===== BRANCH STYLE =====
        for collection in ax.collections:
            if isinstance(collection, LineCollection):
                collection.set_color(BRANCH_COLOR)
                collection.set_linewidth(2.2)
                collection.set_alpha(1.0)

        # ===== CLEAN AXIS =====
        ax.margins(x=0.06, y=0.22)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        for spine in ax.spines.values():
            spine.set_visible(False)

        fig.tight_layout()

        # ✅ สำคัญ: ให้ Streamlit จัดขนาดเอง (มือถือไม่พัง)
        st.pyplot(fig, use_container_width=True)

        # ===== DOWNLOAD =====
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

# ===============================
# RESULTS TAB 
# ===============================
with tab_results:
    st.markdown(f"## {T['identity']}")   
    st.markdown("---")

    df = st.session_state.get("identity_df")
    if df is not None:
        st.dataframe(
            df.style
            .format("{:.2f}")
            .background_gradient(cmap="Blues"),
        width="stretch")

    else:
        st.info(T["no_results"])