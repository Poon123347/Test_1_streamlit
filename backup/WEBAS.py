# streamlit_simple_dna.py
# ------------------------------------------------------------
# Simple DNA Compare (Human vs 4 Apes) — Auto-load local FASTA
# - วางไฟล์ FASTA (เช่น "coi_sequences.fasta") ไว้ในโฟลเดอร์เดียวกับสคริปต์นี้
# - ไม่ต้อง input/อัปโหลดไฟล์ แอปโหลดเองอัตโนมัติ
# - รันด้วย:  py -m streamlit run streamlit_simple_dna.py
# ------------------------------------------------------------

import streamlit as st
from pathlib import Path
from collections import defaultdict
import hashlib
try:
    from Bio import pairwise2
    _HAS_BIO = True
except Exception:
    _HAS_BIO = False

# ---------- หน้าเว็บ / เมตา ----------
st.set_page_config(page_title="Simple DNA Compare", layout="wide")
st.title("เปรียบเทียบความเหมือนของ DNA (อัตโนมัติ) — Human vs 4 Apes")
st.markdown("แอปจะโหลดลำดับ DNA จากไฟล์ FASTA ในโฟลเดอร์เดียวกันโดยอัตโนมัติ แล้วคำนวณ %Identity และ k-mer Jaccard")

# ---------- Utilities ----------
def clean_seq(text: str) -> str:
    """
    รับข้อความ (อาจเป็น FASTA หรือ raw) คืนเป็นสตริง A/C/G/T ตัวใหญ่เท่านั้น
    - ตัดบรรทัดหัว (ขึ้นต้นด้วย '>') ออก
    - ตัดอักขระที่ไม่ใช่ ACGT ออก
    """
    lines = text.splitlines()
    parts = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        parts.append(line)
    seq = "".join(parts).upper()
    return "".join(ch for ch in seq if ch in "ACGT")

def percent_identity_poswise(s1: str, s2: str):
    """
    เปรียบเทียบตำแหน่งต่อ-ตำแหน่ง: นับตำแหน่งที่ตรงกันในช่วงความยาวขั้นต่ำของ s1,s2
    คืนค่า (pid_fraction, matches, compared_length, len_s1, len_s2)
    """
    n1, n2 = len(s1), len(s2)
    if n1 == 0 or n2 == 0:
        return None, 0, 0, n1, n2
    L = min(n1, n2)
    matches = sum(1 for i in range(L) if s1[i] == s2[i])
    pid = matches / L
    return pid, matches, L, n1, n2

# --- แทนที่ฟังก์ชัน percent_identity_poswise เดิมด้วยเวอร์ชันใหม่ด้านล่าง ---
def percent_identity_global(seq1: str, seq2: str):
    """
    คำนวณ %identity ด้วย global alignment (pairwise2.align.globalxx)
    - คืน (pct, matches, denom) โดย pct = matches / max(len(seq1), len(seq2)) * 100
    - ถ้าไม่มี Biopython จะ fallback เป็นวิธีง่าย ๆ (เทียบซ้อนทับตามตำแหน่ง)
    """
    n1, n2 = len(seq1), len(seq2)
    denom = max(n1, n2)
    if denom == 0:
        return None, 0, 0
    if n1 == 0 or n2 == 0:
        return 0.0, 0, denom  # อีกตัวว่าง → 0%

    if _HAS_BIO:
        aln = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
        aln1, aln2 = aln[0], aln[1]
        matches = sum(a == b for a, b in zip(aln1, aln2))
        pct = matches / denom * 100.0
        return pct, matches, denom
    else:
        # Fallback ถ้าไม่มี Biopython: เทียบตำแหน่งซ้อนทับตรงๆ
        L = min(n1, n2)
        matches = sum(seq1[i] == seq2[i] for i in range(L))
        pct = matches / denom * 100.0
        return pct, matches, denom


def kmer_jaccard(s1: str, s2: str, k: int):
    """
    สร้างชุด k-mer สำหรับแต่ละสตริง แล้วคำนวณ Jaccard index = |A∩B| / |A∪B|
    คืนค่า None ถ้าความยาว < k
    """
    if len(s1) < k or len(s2) < k:
        return None
    set1 = set(s1[i:i+k] for i in range(len(s1)-k+1))
    set2 = set(s2[i:i+k] for i in range(len(s2)-k+1))
    if not set1 and not set2:
        return 1.0
    inter = len(set1 & set2)
    uni = len(set1 | set2)
    return inter / uni if uni > 0 else 0.0

def read_fasta_file(path: Path) -> dict[str, str]:
    """
    อ่าน FASTA จากดิสก์ -> {header: cleaned_sequence(ACGT only)}
    """
    records = {}
    hdr, parts = None, []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    records[hdr] = clean_seq("".join(parts))
                hdr = line[1:].strip()
                parts = []
            else:
                parts.append(line)
    if hdr is not None:
        records[hdr] = clean_seq("".join(parts))
    return records

# ---------- คำหลักช่วยแมป header → species ----------
SPECIES = ["Human", "Chimpanzee", "Bonobo", "Gorilla", "Orangutan"]

SPECIES_KEYWORDS = {
    "human": ["homo sapiens", "human", "nc_012920", "rcrs", "hsapiens"],
    "chimpanzee": ["pan troglodytes", "chimpanzee", "chimp", "nc_001643"],
    "bonobo": ["pan paniscus", "bonobo", "nc_001644"],
    "gorilla": ["gorilla gorilla", "gorilla", "nc_001645", "nc_011120"],
    "orangutan": ["pongo", "orangutan", "pongo abelii", "pongo pygmaeus", "nc_002083", "x97707"],
}

def guess_species_from_header(header: str) -> str | None:
    hl = header.lower()
    for sp in SPECIES:
        key = sp.lower()
        if key in hl:
            return sp
        for kw in SPECIES_KEYWORDS.get(key, []):
            if kw in hl:
                return sp
    return None

# ---------- กำหนดชื่อไฟล์ที่ต้องโหลด (อยู่โฟลเดอร์เดียวกับสคริปต์) ----------
APP_DIR = Path(__file__).parent
FASTA_FILENAME = "coi_sequences.fasta"   # ← เปลี่ยนชื่อไฟล์ที่นี่ถ้าคุณใช้ชื่ออื่น
FASTA_PATH = APP_DIR / FASTA_FILENAME

@st.cache_data(show_spinner=False)
def load_records() -> dict[str, str]:
    if not FASTA_PATH.exists():
        return {}
    return read_fasta_file(FASTA_PATH)

records = load_records()

# ---------- Sidebar: มีเฉพาะตัวเลือก k ----------
st.sidebar.header("Options")
k = st.sidebar.slider("k (k-mer) for Jaccard", min_value=3, max_value=12, value=9)

# ---------- ตรวจไฟล์ / แมปให้ 5 สปีชีส์อัตโนมัติ ----------
if not records:
    st.error(
        f"ไม่พบไฟล์ FASTA: **{FASTA_PATH.name}** ในโฟลเดอร์ `{APP_DIR}`\n"
        f"→ กรุณาวางไฟล์ชื่อเดียวกันไว้ในโฟลเดอร์เดียวกับสคริปต์"
    )
    st.stop()

st.success(f"พบไฟล์ FASTA: **{FASTA_PATH.name}** — โหลด {len(records)} record(s)")

headers = list(records.keys())

# สร้าง mapping แบบ “ไม่ซ้ำ header” และเลือกให้ตรงสายพันธุ์มากที่สุด
unused = headers.copy()
mapping: dict[str, str] = {}

# รอบแรก: พยายามหา header ที่แมตช์โดยตรงตามคำหลัก
for sp in SPECIES:
    pick = None
    for h in list(unused):
        if guess_species_from_header(h) == sp:
            pick = h
            break
    if pick:
        mapping[sp] = pick
        unused.remove(pick)

# รอบสอง: ถ้าเหลือสปีชีส์ไหนยังไม่ได้แมป ให้หยิบ header ที่ยังไม่ใช้มาเติม (เพื่อให้ครบ 5)
for sp in SPECIES:
    if sp not in mapping:
        if not unused:
            break
        mapping[sp] = unused.pop(0)

# เตรียมลำดับตาม mapping
seqs: dict[str, str] = {sp: records.get(mapping.get(sp, ""), "") for sp in SPECIES}

# ---------- Debug: ตรวจจับว่ามีลำดับซ้ำ ๆ หรือไม่ ----------
hash_to_spp = defaultdict(list)
for sp in SPECIES:
    s = seqs[sp]
    md5 = hashlib.md5(s.encode()).hexdigest() if s else None
    hash_to_spp[md5].append(sp)

dups = [grp for grp in hash_to_spp.values() if len(grp) > 1 and None not in grp]
if dups:
    st.warning(
        "พบลำดับเหมือนกันเป๊ะในหลายสปีชีส์: " + ", ".join([str(grp) for grp in dups]) +
        "\nโปรดตรวจไฟล์ FASTA ว่า header ถูกต้อง/ไม่ซ้ำ และเป็นลำดับคนละชนิดจริง"
    )

# แสดง mapping และความยาว
st.markdown(
    "**Mapping อัตโนมัติ (Species → FASTA header):**  " +
    ", ".join(f"{sp} → `{mapping.get(sp, '-')}`" for sp in SPECIES)
)
st.markdown(
    "**Sequence lengths (bp):**  " +
    ", ".join(f"{sp}: {len(seqs[sp])}" for sp in SPECIES)
)

# --- ในส่วน "คำนวณผล" ให้แก้การเรียกใช้เป็นแบบใหม่ ---
human = seqs["Human"]
if len(human) == 0:
    st.warning("กรุณาใส่ลำดับของ Human ก่อน (สามารถใช้ตัวอย่างเพื่อทดลองได้)")
else:
    st.subheader("ผลการเปรียบเทียบ (Human เทียบกับแต่ละสายพันธุ์)")
    rows = []
    for sp in SPECIES[1:]:
        other = seqs[sp]
        pct, matches, denom = percent_identity_global(human, other)
        if pct is None:
            pid_str = "ไม่สามารถคำนวณ (ลำดับว่าง)"
        else:
            pid_str = f"{pct:.2f}%  ({matches}/{denom})"

        jk = kmer_jaccard(human, other, k)
        jk_str = f"{jk*100:.2f}%" if jk is not None else f"ไม่พอความยาว (ต้อง >= {k})"

        rows.append({
            "Species": sp,
            "Len (Human)": len(human),
            "Len (other)": len(other),
            "Pos-wise %ID (global)": pid_str,
            f"k-mer Jaccard (k={k})": jk_str
        })
    st.table(rows)

    # แจ้งเตือนถ้าไม่มี Biopython เพื่อให้ได้ผล global alignment ที่แท้จริง
    if not _HAS_BIO:
        st.info("กำลังใช้ fallback (ไม่มี Biopython). ติดตั้งด้วย:  pip install biopython  เพื่อใช้ global alignment จริง ๆ")


# ---------- คำเตือน / ข้อจำกัด ----------
st.markdown(
    """
**คำเตือน / ข้อจำกัด (สำคัญ):**
- ตอนนี้ใช้ **global alignment** (Needleman–Wunsch ผ่าน `pairwise2.align.globalxx`) โดยให้คะแนนแบบง่าย: **match=1, mismatch=0, gap=0**  
  ⇒ ช่องว่าง (indel) **ไม่ถูกปรับโทษ** ทำให้ %identity อาจ **สูงเกินจริง** เมื่อมี indel ยาว ๆ
- ถ้าไม่มี Biopython โค้ดจะ **fallback เป็นการเทียบตำแหน่งซ้อนทับ** เฉย ๆ (ไม่ใช่การจัดเรียงลำดับจริง)
- **k-mer Jaccard** เป็นวิธี alignment-free ที่ดีสำหรับลำดับยาวไม่เท่ากัน แต่ไม่เห็นการจับคู่ฐานต่อฐาน
- งานวิจัยจริงมักใช้สคอร์ที่สมจริงกว่า (เช่น `pairwise2.align.globalms` พร้อม **affine gap penalties**) หรือทำ **multiple sequence alignment** (MAFFT ฯลฯ) ก่อนวิเคราะห์
    """
)