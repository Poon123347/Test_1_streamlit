import matplotlib
matplotlib.use("Agg")
import re
import time
import streamlit as st
import pandas as pd
from Bio import Entrez, SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from io import BytesIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from copy import deepcopy
from matplotlib.collections import LineCollection
plt.rcParams["toolbar"] = "None"
from matplotlib import font_manager, rcParams
import os
import numpy as np
import matplotlib.colors as mcolors
from urllib.error import HTTPError, URLError

FONT_PATH = os.path.abspath("fonts/NotoSansThai-Regular.ttf")

if not os.path.exists(FONT_PATH):
    raise FileNotFoundError(f"Font not found: {FONT_PATH}")

font_manager.fontManager.addfont(FONT_PATH)
font_prop = font_manager.FontProperties(fname=FONT_PATH)

rcParams["font.family"] = font_prop.get_name()
rcParams["axes.unicode_minus"] = False

# ========== PAGE CONFIG ==========
st.set_page_config(
    page_title="ตัวสำรวจสายวิวัฒนาการ ND5",
    page_icon="🧬",
    layout="centered",
    initial_sidebar_state="expanded"  # 👈 IMPORTANT
)

language = st.sidebar.selectbox("Language/ภาษา", ["English", "ภาษาไทย"], index=1)

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
}

THEME = DARK 
cmap = plt.get_cmap("Greens")

st.markdown(f"""
<style>

/* ================== THEME VARIABLES ================== */
:root {{
    --bg: {THEME["BG"]};
    --text: #e5e7eb;
    --text-muted: {THEME["TEXT_MUTED"]};
    --sidebar-bg: {THEME["SIDEBAR_BG"]};
    --panel-bg: #1e293b;
    --border: {THEME["BORDER"]};
    --btn-bg: {THEME["BTN_BG"]};
    --btn-text: {THEME["BTN_TEXT"]};
    --code-bg: {THEME["CODE_BG"]};
    --code-text: {THEME["CODE_TEXT"]};
}}

/* ================== GLOBAL ================== */
.stApp {{
    background-color: #0f172a !important;
    color: var(--text) !important;
    padding-top: 0 !important;
}}

html, body, div, span, p, label,
h1, h2, h3, h4, h5, h6 {{
    color: var(--text) !important;
}}

header[data-testid="stHeader"] {{
    display: block !important;
    background: #0f172a !important; 
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

/* ================== ARROWS / ICONS ================== */
div[data-testid="stExpander"] summary svg {{
    fill: var(--text) !important;
}}
div[data-baseweb="select"] svg {{
    fill: var(--text) !important;
}}

/* ================== TYPOGRAPHY / SPACING ================== */
.stApp > div {{
    padding-left: 1rem !important;
    padding-right: 1rem !important;
}}
div[data-testid="stMetricValue"] {{
    font-size: 20px !important;
    font-weight: 500 !important;
    letter-spacing: -0.3px !important;
}}
div[data-testid="stMetricLabel"] {{
    font-size: 11px !important;
    text-transform: uppercase !important;
    letter-spacing: 0.05em !important;
    color: {THEME["TEXT_MUTED"]} !important;
}}
div[data-testid="stTabs"] button {{
    font-size: 13px !important;
    font-weight: 500 !important;
    letter-spacing: 0.01em !important;
}}
h1, h2, h3 {{
    letter-spacing: -0.4px !important;
    font-weight: 500 !important;
}}

""", unsafe_allow_html=True)

# ========== NCBI CONFIG ==========
Entrez.email = "poonthakorn@gmail.com"

species_accessions = {
    # === Hominins ===
    "Human (Homo sapiens)": "NC_012920.1",
    "Neanderthal (Homo neanderthalensis)": "NC_011137.1",
    "Denisovan (Homo sp. Denisova)": "NC_013993.1",

    # === Great Apes ===
    "Chimpanzee (Pan troglodytes)": "NC_001643.1",
    "Bonobo (Pan paniscus)": "NC_001644.1",
    "Western Gorilla (Gorilla gorilla)": "NC_011120.1",
    "Eastern Gorilla (Gorilla beringei)": "NC_037853.1", 
    "Bornean Orangutan (Pongo pygmaeus)": "NC_001646.1",
    "Sumatran Orangutan (Pongo abelii)": "NC_002083.1",

    # === Lesser Apes ===
    "Common Gibbon (Hylobates lar)": "NC_002082.1",

    # === Old World Monkeys ===
    "Rhesus Macaque (Macaca mulatta)": "NC_005943.1",
    "Barbary Macaque (Macaca sylvanus)": "NC_002764.1",
    "Black-and-White Colobus Monkeys (Colobus guereza)": "NC_006901.1",
    "Green Monkey (Chlorocebus sabaeus)": "NC_008066.1",
    "Proboscis Monkey (Nasalis larvatus)": "NC_008216.1",
    
}

SPECIES_TH = {
    "Human": "มนุษย์",
    "Neanderthal": "มนุษย์นีแอนเดอร์ธัล",
    "Denisovan": "มนุษย์เดนิโซวาน",

    "Chimpanzee": "ลิงชิมแปนซี",
    "Bonobo": "ลิงโบโนโบ",

    "Eastern Gorilla": "ลิงกอริลลาตะวันออก",
    "Western Gorilla": "ลิงกอริลลาตะวันตก",

    "Bornean Orangutan": "ลิงอุรังอุตังบอร์เนียว",
    "Sumatran Orangutan": "ลิงอุรังอุตังสุมาตรา",

    "Common Gibbon": "ชะนี",

    "Rhesus Macaque": "ลิงวอก",
    "Barbary Macaque": "ลิงบาร์บารี",

    "Green Monkey": "ลิงเขียว",
    "Proboscis Monkey": "ลิงจมูกยาว",
    "Black-and-White Colobus Monkeys": "ลิงโคโลบัสขาวดำ",
}


# ========== LANGUAGE SUPPORT ==========
LANGS = {
    "ภาษาไทย": {
        # ===== Header =====
        "title": "ตัวสำรวจสายวิวัฒนาการ ND5",
        "subtitle": "เปรียบเทียบ ND5 ระหว่างสายพันธุ์ : แผนภูมิวิวัฒนาการ และตารางความคล้ายคลึงของลำดับดีเอ็นเอ",

        # ===== Controls =====
        "select_species": "กรุณาเลือกสายพันธุ์ตัวอย่าง (อย่างน้อย 2 ชนิด)",
        "choose_options": "เลือกสายพันธุ์",
        "available": "สายพันธุ์ทั้งหมด",
        "selected": "สายพันธุ์ที่เลือก",
        "fetch": "ดึงข้อมูล",

        # ===== Info =====
        "length": "ความยาว ND5 (คู่เบส)",
        "preview": "ดูตัวอย่างลำดับเบส",
        "workflow_title": "ขั้นตอนการใช้งาน",
        "workflow_step1": "เลือกอย่างน้อย 2 สายพันธุ์ ที่อยากนำมาเปรียบเทียบจากแถบด้านล่าง",
        "workflow_step2": "กดปุ่มดึงข้อมูล",
        "workflow_step3": "ดูตารางความคล้ายคลึงของลำดับดีเอ็นเอ และแท็บวิวัฒนาการด้านล่าง",
        "fetch_help": "ค้นหาสปีชีส์โดยชื่อ แล้วเลือกอย่างน้อยสองชนิด",
        "control_info": "เลือกสายพันธุ์ที่ต้องการแล้วกดดึงข้อมูล เพื่อดึงข้อมูลจาก NCBI",

        # ===== Tabs =====
        "results": "ผลลัพธ์",
        "tree": "วิวัฒนาการ",
        # ===== Outputs =====
        "identity": "ตารางความคล้ายคลึงของลำดับดีเอ็นเอ (%)",
        "phylo_tree": "แผนภูมิวิวัฒนาการ (ND5)",

        # ===== States =====
        "spinner": "กำลังดึงข้อมูล...",
        "no_results": "กรุณาดึงข้อมูล ND5 ก่อนเพื่อดูผลลัพธ์",
        "no_tree": "กรุณาดึงข้อมูล ND5 ก่อนเพื่อดูผลลัพธ์",

        # ===== Messages =====
        "missing_nd5": "⚠️ สายพันธุ์ต่อไปนี้ไม่มีข้อมูล ND5 ในฐานข้อมูล NCBI:",

        # ===== Expanders =====

        "tree_description": "แผนภาพนี้แสดงความสัมพันธ์ทางวิวัฒนาระหว่างสายพันธุ์ตามลำดับ DNA ND5 สายพันธุ์ที่เชื่อมต่อกันดีโดยกิ่งที่สั้นกว่านั้นมีความสัมพันธุ์ที่ใกล้เคียงกัน ซึ่งหมายว่า DNA มีความคล้ายคลึงมาก ระยะห่างที่ยาวขึ้นระหว่างกิ่งน้อยต่อบ่งชี้ถึงความแตกต่างทางพันธุ์มากขึ้น ตัวเลขบนเส้นแสดงภาวะห่างของระยะห่างทางพันธุ์สาภียน (สีเขียว = ใกล้เคียง, สีแดง = ห่างกันมาก)",
        "branch_length_note": "ความยาวกิ่งสะท้อนความแตกต่างของลำดับ ND5 โดยถ้าความยาวกิ่งระหว่างสายพันธุ์ไม่ต้างจากกันมาก หมายถึงสายพันธุ์นั้นๆมึความใกล้เคียงกันทางวิวัฒนาการ และความยาวกิ่งระหว่างสายพันธุ์กอยู่ห่างจากกัน หมายถึงสายพันธุ์นั้นๆมึระยะห่างทางวิวัฒนาการ",
        "branch_color_note": "",
        "branch_number_note": "",

        "taxa_axis": "สายพันธุ์",
        "branch_axis": "ความยาวกิ่ง",

        "download_tree": "ดาวน์โหลดแผนภูมิวิวัฒนาการ (PNG)",
        "download_nd5": "ดาวน์โหลดลำดับ ND5 (FASTA)",
        "species" : "สายพันธุ์"
    },

    "English": {
        # ===== Header =====
        "title": "ND5 Phylogenetic Explorer",
        "subtitle": "Compare ND5 across species – Phylogenetic tree + DNA Sequence Similarity Matrix",

        # ===== Controls =====
        "select_species": "Please select species (min 2)",
        "choose_options": "Choose species",
        "available": "Available species",
        "selected": "Selected",
        "fetch": "Fetch",

        # ===== Info =====
        "length": "ND5 length (Base pairs)",
        "preview": "Preview sequences",
        "workflow_title": "Quick workflow",
        "workflow_step1": "Select at least two species.",
        "workflow_step2": "Click Fetch to load ND5 data.",
        "workflow_step3": "View the DNA Sequence Similarity Matrix and Phylogenetic tab.",
        "fetch_help": "Search species by name and pick at least two for comparison.",
        "control_info": "Use the list to choose species, then fetch ND5 sequences from NCBI.",

        # ===== Tabs =====
        "results": "Results",
        "tree": "Phylogenetic",
        # ===== Outputs =====
        "identity": "DNA Sequence Similarity Matrix (%)",
        "phylo_tree": "Phylogenetic Tree (ND5)",

        # ===== States =====
        "spinner": "Fetching data...",
        "no_results": "Fetch ND5 to see results.",
        "no_tree": "Fetch ND5 to generate tree.",

        # ===== Messages =====
        "missing_nd5": "The following species do not have ND5 data in NCBI:",

        # ===== Expanders =====
        "tree_description": "This diagram shows the evolutionary relationships between species based on their ND5 DNA sequences. Species that are connected by shorter branches are more closely related, meaning their DNA is more similar. Longer distances between branches indicate greater genetic differences. The numbers on the lines represent the relative genetic distance between species (green branches = shorter, red branches = longer).",
        "branch_length_note": "Branch length reflects relative ND5 sequence difference — shorter means more similar, longer means more different.",
        "branch_color_note": "",
        "branch_number_note": "",

        "taxa_axis": "taxa",
        "branch_axis": "branch length",

        "download_tree": "Download phylogenetic tree (PNG)",
        "download_nd5": "Download ND5 sequences (FASTA)",

        "species" : "species"
    }
}

# ========== SIDEBAR ==========
# ===============================

T = LANGS[language]

# map full accession names to display labels and keep a short-name lookup
species_short_to_full = {}

for full_name in species_accessions:
    short_name = re.sub(r"\s*\(.*\)$", "", full_name)
    species_short_to_full[short_name] = full_name

if language == "ภาษาไทย":
    species_display = {
        full_name: SPECIES_TH.get(re.sub(r"\s*\(.*\)$", "", full_name), short_name)
        for full_name, short_name in ((full_name, re.sub(r"\s*\(.*\)$", "", full_name)) for full_name in species_accessions)
    }
else:
    species_display = {
        full_name: re.sub(r"\s*\(.*\)$", "", full_name)
        for full_name in species_accessions
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

st.markdown(f"""
<div style="padding: 1.5rem 0 1rem; border-bottom: 0.5px solid {THEME['BORDER']}; margin-bottom: 1.5rem;">
  <div style="display:flex; align-items:center; gap:10px; margin-bottom:6px;">
    <div style="width:32px;height:32px;border-radius:8px;background:#238636;
      display:flex;align-items:center;justify-content:center;
      color:#fff;font-size:13px;font-weight:600;letter-spacing:-0.5px;">N5</div>
    <span style="font-size:20px;font-weight:500;color:{THEME['TEXT']};letter-spacing:-0.3px;">{T['title']}</span>
  </div>
  <p style="font-size:13px;color:{THEME['TEXT_MUTED']};margin:0;padding-left:42px;">{T['subtitle']}</p>
</div>
""", unsafe_allow_html=True)

st.markdown(
    f"""
    <div style="padding: 0.75rem 1rem; border: 1px solid {THEME['BORDER']}; background: {THEME['PANEL_BG']}; border-radius: 16px; margin-bottom: 1.5rem;">
      <strong style="font-size:13px;color:{THEME['TEXT']};">{T['workflow_title']}</strong>
      <div style="font-size:13px;color:{THEME['TEXT_MUTED']}; margin-top:0.5rem; line-height:1.5;">
        • {T['workflow_step1']}<br>
        • {T['workflow_step2']}
      </div>
    </div>
    """,
    unsafe_allow_html=True,
)

# ========== UI LAYOUT: control panel ==========
with st.container():

    # --- two-column row ---
    controls_col, info_col = st.columns([2.2, 1])

    # ===== LEFT: species selector + fetch button =====
    with controls_col:

        CATEGORIES = {
            "Hominins" if language == "English" else "กลุ่มสายพันธุ์มนุษย์": [
                "Human", "Neanderthal", "Denisovan",
            ],

            "Closest relatives to human" if language == "English" else "ญาติใกล้ชิดมนุษย์": [
                "Chimpanzee", "Bonobo",
            ],

            "Great apes" if language == "English" else "ลิงใหญ่": [
                "Eastern Gorilla", "Western Gorilla",
            ],

            "Orangutans" if language == "English" else "อุรังอุตัง": [
                "Bornean Orangutan", "Sumatran Orangutan",
            ],

            "Lesser apes" if language == "English" else "ลิงเล็ก": [
                "Common Gibbon",
            ],

            "Old World monkeys" if language == "English" else "ลิงโลกเก่า": [
                "Rhesus Macaque",
                "Barbary Macaque",
                "Green Monkey",
                "Proboscis Monkey",
                "Black-and-White Colobus Monkeys",
            ],
        }

        st.markdown(f"**{T['select_species']}**")

        selected_short_names = []

        with st.form(key="species_form"):

            for category, members in CATEGORIES.items():
                with st.expander(category, expanded=False):

                    for short_name in members:
                        full_name = species_short_to_full.get(short_name)
                        if not full_name:
                            continue

                        label = species_display[full_name]

                        # unique key per checkbox
                        key = f"chk_{short_name}"

                        if st.checkbox(label, key=key):
                            selected_short_names.append(short_name)

            fetch_clicked = st.form_submit_button(T["fetch"])

        selected_species = [
            species_short_to_full[short_name]
            for short_name in selected_short_names
            if short_name in species_short_to_full
        ]

    # ===== RIGHT: info panel =====
    with info_col:
        st.markdown(
            f"""
            <div style="padding:1rem; border:1px solid {THEME['BORDER']}; border-radius:16px; background:{THEME['PANEL_BG']};">
              <p style="margin:0 0 0.75rem; font-size:13px; color:{THEME['TEXT_MUTED']};">{T['control_info']}</p>
              <div style="display:grid; gap:0.75rem;">
                <div style="padding:0.8rem; background:#111827; border-radius:12px;">
                  <strong style="font-size:12px; color:{THEME['TEXT']};">{T['available']}</strong>
                  <div style="font-size:20px; color:{THEME['TEXT']};">{len(species_accessions)}</div>
                </div>
                <div style="padding:0.8rem; background:#111827; border-radius:12px;">
                  <strong style="font-size:12px; color:{THEME['TEXT']};">{T['selected']}</strong>
                  <div style="font-size:20px; color:{THEME['TEXT']};">{len(selected_species)}</div>
                </div>
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        # ========== CACHING FOR FETCH ==========
@st.cache_data(show_spinner=False)
def extract_nd5_cached(accession):
    def feature_is_nd5(feature):
        gene = feature.qualifiers.get("gene", [""])[0].lower().strip()
        qualifiers = " ".join(
            " ".join(values) for values in feature.qualifiers.values()
        ).lower()

        if gene in ["nd5", "nad5", "mt-nd5", "mtnd5", "nadh5"]:
            return True

        return any(
            term in qualifiers for term in [
                "nad5",
                "nd5",
                "nadh5",
                "nadh dehydrogenase subunit 5",
                "subunit 5"
            ]
        )

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
            if feature.type == "CDS" and feature_is_nd5(feature):
                return feature.extract(record.seq)

        # ── fallback: try gene features too (some records skip CDS) ──
        for feature in record.features:
            if feature.type == "gene" and feature_is_nd5(feature):
                return feature.extract(record.seq)

        # ── fallback: search broader feature types for ND5 annotation ──
        for feature in record.features:
            if feature.type in {"misc_feature", "region", "exon", "mat_peptide", "rRNA", "tRNA", "source"} and feature_is_nd5(feature):
                return feature.extract(record.seq)

        return None

    except (ValueError, IOError, HTTPError, URLError) as error:
        print(error)  # or st.error(error) for debugging
        return None



@st.cache_data(show_spinner=False)
def build_nd5_fasta_download(accessions, species_names):
    records = []
    for accession, name in zip(accessions, species_names):
        try:
            nd5_seq = extract_nd5_cached(accession)
            if nd5_seq is None:
                continue
            records.append(
                SeqRecord(
                    nd5_seq,
                    id=name,
                    description=f"ND5 from {accession}"
                )
            )
            time.sleep(0.4)
        except Exception as e:
            print(f"Error building FASTA for {accession}: {e}")

    output = BytesIO()
    SeqIO.write(records, output, "fasta")
    return output.getvalue()


# ========== HELPERS ==========
def percent_identity(seq1, seq2):
    matches, valid = 0, 0
    for base1, base2 in zip(seq1, seq2):
        if base1 == "-" or base2 == "-":
            continue
        valid += 1
        if base1 == base2:
            matches += 1
    return round(matches / valid * 100, 2) if valid > 0 else 0.0

def translate_species(species_name, selected_language):
    if selected_language == "ภาษาไทย":
        return SPECIES_TH.get(species_name, species_name)
    return species_name


def cell_color(cell_value):
    if not isinstance(cell_value, float):
        return THEME["SIDEBAR_BG"], THEME["TEXT"]
    if cell_value == 100.0:
        return "#165423", "#ffffff"
    branch_color_norm = max(0.15, min(0.85, (cell_value - 20) / 80))
    red, green, blue, *_ = cmap(branch_color_norm)
    background_color = mcolors.to_hex((red, green, blue))
    brightness = 0.299 * red + 0.587 * green + 0.114 * blue
    text_color = "#0d1117" if brightness > 0.55 else "#e6edf3"
    return background_color, text_color


def calc_gc(sequence):
    sequence = str(sequence).upper().replace("-", "")
    if len(sequence) == 0:
        return 0.0
    return round((sequence.count("G") + sequence.count("C")) / len(sequence) * 100, 2)


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
        raw_nd5_seqs = nd5_seqs.copy()
                    # ===== FORCE SAME ND5 LENGTH =====
        
        if nd5_seqs:
            max_len = max(len(seq) for seq in nd5_seqs.values())
            for k in nd5_seqs:
                nd5_seqs[k] = nd5_seqs[k] + "-" * (max_len - len(nd5_seqs[k]))

        st.session_state["nd5_lengths"] = nd5_lengths
        st.session_state["raw_nd5_seqs"] = raw_nd5_seqs
        st.session_state["valid_species"] = list(nd5_seqs.keys())
        st.session_state["nd5_seqs"] = nd5_seqs
        
        if missing_species:
            st.warning(
                f"{T['missing_nd5']}\n\n"
                + ", ".join(missing_species)
            )
        # AUTO-COMPUTE SIMILARITY MATRIX
        names = list(nd5_seqs.keys())
        matrix = []
        for i, name_a in enumerate(names):
            row = []
            for j, name_b in enumerate(names):
                row.append(100.0 if i == j else percent_identity(nd5_seqs[name_a], nd5_seqs[name_b]))
            matrix.append(row)

        st.session_state["identity_df"] = pd.DataFrame(
            matrix, index=names, columns=names
        )

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
        mask = ~np.eye(len(nd5_seqs), dtype=bool)
        avg_identity = st.session_state["identity_df"].values[mask].mean()
        st.session_state["avg_identity"] = round(avg_identity, 1)

# ========== SHOW ND5 LENGTHS (compact) ==========
if "nd5_seqs" in st.session_state:
    seqs = st.session_state["nd5_seqs"]

    # ── divider label ──
    st.markdown(f"""
    <p style="font-size:11px;font-weight:500;color:{THEME['TEXT_MUTED']};
      text-transform:uppercase;letter-spacing:0.06em;margin:1.5rem 0 0.5rem;">
      {'Overview' if language == 'English' else 'ภาพรวม'}
    </p>
    """, unsafe_allow_html=True)

    # ── 2-column stats row ──
    s1, s2 = st.columns(2)

    # Average identity
    avg_identity = st.session_state.get("avg_identity")
    avg_value = f"{avg_identity}%" if avg_identity is not None else "—"

    s1.metric(
        "Avg similarity" if language == "English" else "ความคล้ายคลึงเฉลี่ย",
        avg_value
    )

    # ND5 length (range)
    lengths = st.session_state.get("nd5_lengths", {})

    if lengths:
        min_len = min(lengths.values())
        max_len = max(lengths.values())
        value = f"{min_len}-{max_len} bp" if min_len != max_len else f"{min_len} bp"
    else:
        value = "—"

    s2.metric(
        "ND5 length" if language == "English" else "ความยาว ND5",
        value
    )

    st.markdown("<div style='margin-bottom:1rem;'></div>", unsafe_allow_html=True)

    # ── per-species lengths ──
    st.markdown(f"""
    <p style="font-size:11px;font-weight:500;color:{THEME['TEXT_MUTED']};
      text-transform:uppercase;letter-spacing:0.06em;margin:0 0 0.5rem;">
      {T['length']}
    </p>
    """, unsafe_allow_html=True)

    cols = st.columns(4)
    for i, (name, _) in enumerate(seqs.items()):
        real_len = st.session_state["nd5_lengths"][name]
        cols[i % 4].metric(name, f"{real_len} bp")

    with st.expander(T["preview"]):
        for name, seq in seqs.items():
            st.markdown(f"**{name}**")
            st.code(seq[:120] + ("..." if len(seq) > 120 else ""))

        # Translated download button at the bottom of the preview
        try:
            raw_seqs = st.session_state.get("raw_nd5_seqs", seqs)
            fasta_buf = BytesIO()
            for name, seq in raw_seqs.items():
                fasta_buf.write(f">{name}\n".encode("utf-8"))
                s = str(seq).replace("-", "")
                fasta_buf.write(s.encode("utf-8"))

            st.download_button(
                T["download_nd5"],
                data=fasta_buf.getvalue(),
                file_name="nd5_sequences.fasta",
                mime="text/x-fasta",
                key="download_nd5_preview",
            )
        except Exception:
            pass

# ========== TABS FOR RESULTS / TREE ==========
tab_results, tab_tree = st.tabs([
    T["results"],
    T["tree"]
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

        bounds = [0, 70, 85, 95, 98, 100]
        colors = [
            "#d73027",
            "#fc8d59",
            "#fee08b",
            "#91cf60",
            "#1a9850",
        ]

        cmap = mcolors.ListedColormap(colors)
        norm = mcolors.BoundaryNorm(bounds, cmap.N)

        fig, ax = plt.subplots(
            figsize=(max(12, len(df_display) * 0.9), max(8, len(df_display) * 0.6)),
            dpi=120,
            facecolor=THEME["PANEL_BG"],
        )
        fig.patch.set_facecolor(THEME["PANEL_BG"])
        ax.set_facecolor(THEME["PANEL_BG"])

        img = ax.imshow(df_display.values, cmap=cmap, norm=norm, aspect='auto')

        ax.set_xticks(np.arange(df_display.shape[1]))
        ax.set_yticks(np.arange(df_display.shape[0]))
        ax.set_xticklabels(df_display.columns, rotation=45, ha='right', fontsize=12)
        ax.set_yticklabels(df_display.index, fontsize=12)
        ax.set_title(T["identity"], pad=12, fontsize=18, color="#9CA3AF")

        for i in range(df_display.shape[0]):
            for j in range(df_display.shape[1]):
                val = df_display.iat[i, j]
                rgba = cmap(norm(val))
                brightness = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
                text_color = "black" if brightness > 0.55 else "white"
                ax.text(
                    j,
                    i,
                    f"{val:.2f}%",
                    ha='center',
                    va='center',
                    color=text_color,
                    fontsize=10,
                )

        # Remove the side colorbar and adjust margins to better center the heatmap.
        fig.subplots_adjust(left=0.04, right=0.99, top=0.92, bottom=0.30)
        ax.tick_params(axis='x', labelsize=12, colors="#9CA3AF")
        ax.tick_params(axis='y', labelsize=12, colors="#9CA3AF")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xlabel("", color="#9CA3AF")
        ax.set_ylabel("", color="#9CA3AF")

        ax.set_xlim(-0.5, df_display.shape[1] - 0.5)
        ax.set_ylim(df_display.shape[0] - 0.5, -0.5)
        ax.grid(False)
        fig.tight_layout()

        st.pyplot(fig)
        plt.close(fig)

        # ===== PAIRWISE DETAIL VIEW =====
        if df is not None and "nd5_seqs" in st.session_state:
            st.markdown("---")

            # ===== INSTRUCTION =====
            st.caption(
                "Select two species below to compare their ND5 similarity score."
                if language == "English"
                else "เลือกสปีชีส์สองชนิดเพื่อเปรียบเทียบค่าความคล้ายคลึงของ ND5"
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
                        if score == 100:
                            delta_label = "Identical"
                            delta_desc = "No sequence differences"
                        elif score >= 95:
                            delta_label = "Highly similar"
                            delta_desc = "Very closely related"
                        elif score >= 80:
                            delta_label = "Moderately similar"
                            delta_desc = "Some evolutionary divergence"
                        else:
                            delta_label = "Divergent"
                            delta_desc = "Distant evolutionary relationship"
                    else:
                        if score == 100:
                            delta_label = "เหมือนกันทุกประการ"
                            delta_desc = "ไม่มีความแตกต่างในลำดับ"
                        elif score >= 95:
                            delta_label = "คล้ายคลึงมาก"
                            delta_desc = "มีความสัมพันธ์ใกล้ชิดมาก"
                        elif score >= 80:
                            delta_label = "คล้ายคลึงปานกลาง"
                            delta_desc = "มีความแตกต่างทางวิวัฒนาการบางส่วน"
                        else:
                            delta_label = "แตกต่าง"
                            delta_desc = "มีความสัมพันธ์ทางวิวัฒนาการที่ห่างไกล"

                    # Determine color based on similarity
                    # 100% = Dark Green, 95-99% = Green, 80-94% = Orange, <80% = Red
                    if score == 100:
                        color_hex = "#165423"  # Dark Green
                    elif score >= 95:
                        color_hex = "#238636"  # Green
                    elif score >= 80:
                        color_hex = "#f97316"  # Orange
                    else:
                        color_hex = "#dc2626"  # Red

                    # Custom metric with exact colors
                    st.markdown(f"""
                    <div style="
                        background: {THEME['PANEL_BG']};
                        border: 1px solid {THEME['BORDER']};
                        border-radius: 12px;
                        padding: 16px 20px;
                        margin: 8px 0;
                    ">
                        <div style="
                            font-size: 13px;
                            color: {THEME['TEXT_MUTED']};
                            text-transform: uppercase;
                            letter-spacing: 0.05em;
                            margin-bottom: 4px;
                        ">
                            {species_a} ↔ {species_b}
                        </div>
                        <div style="
                            font-size: 32px;
                            font-weight: 600;
                            color: {THEME['TEXT']};
                            margin: 4px 0;
                        ">
                            {score:.2f}%
                        </div>
                        <div style="
                            display: inline-flex;
                            align-items: center;
                            gap: 6px;
                            background: {color_hex}20;
                            border: 1px solid {color_hex};
                            border-radius: 20px;
                            padding: 4px 12px;
                            font-size: 13px;
                            font-weight: 500;
                            color: {color_hex};
                        ">
                            {delta_label}
                        </div>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    # Show interpretation below
                    st.caption(f"**{delta_label}**: {delta_desc}")
    else:
        st.info(T["no_results"])

# ===============================
# TREE TAB
# ===============================
with tab_tree:
    st.subheader(T["phylo_tree"])

    tree = deepcopy(st.session_state.get("tree"))

    if tree is None:
        st.info(T["no_tree"])
    else:
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

        n = len(tree.get_terminals())
        fig_height = max(5, min(20, n * 0.32))
        font_size = max(6, min(12, int(16 - n * 0.16)))

        fig = plt.figure(
            figsize=(10, fig_height),
            facecolor=THEME["PANEL_BG"]
        )
        fig.patch.set_facecolor(THEME["PANEL_BG"])
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
                collection.set_linewidth(2.4)
                collection.set_alpha(1.0)

        ax.margins(x=0.0, y=0.22)
        x0, x1 = ax.get_xlim()
        ax.set_xlim(left=max(0, x0 * 0.98), right=x1 * 1.02)
        ax.set_xlabel(
            T["branch_axis"],
            color=TEXT_COLOR,
            fontproperties=font_prop,
            labelpad=10,
            fontsize=9
        )
        ax.set_ylabel(
            T["taxa_axis"],
            color=TEXT_COLOR,
            fontproperties=font_prop,
            labelpad=10,
            fontsize=9
        )
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

        terminals = tree.get_terminals()
        if terminals:
            root_y = (len(terminals) - 1) / 2
            curr_xlim = ax.get_xlim()
            ax.set_xlim(left=-0.04, right=curr_xlim[1] + 0.02)

        for spine in ax.spines.values():
            spine.set_visible(False)

        fig.tight_layout()
        fig.subplots_adjust(left=0.04, right=0.98, top=0.95, bottom=0.08)
        st.pyplot(fig, width = "stretch")
        st.caption(T["branch_length_note"])

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
                T["download_tree"],
                data=buf.getvalue(),
                file_name="nd5_phylogenetic_tree.png",
                mime="image/png"
            )
