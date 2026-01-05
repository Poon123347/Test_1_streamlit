# streamlit_simple_dna.py
import streamlit as st

# --- หน้าเว็บ / เมตา ---
st.set_page_config(page_title="Simple DNA Compare", layout="wide")
st.title("เปรียบเทียบความเหมือนของ DNA (ง่าย ๆ) — Human vs 4 Apes")
st.markdown("วางลำดับ DNA (FASTA หรือ raw) ลงในกล่องด้านล่าง แล้วเลือกค่า `k` สำหรับ k-mer Jaccard")

# --- ฟังก์ชันช่วยเหลือ (utilities) ---
def clean_seq(text):
    """
    รับข้อความ (อาจเป็น FASTA หรือ raw) คืนเป็นสตริงของตัวอักษร A/C/G/T ตัวใหญ่เท่านั้น
    - เอา header บรรทัดที่ขึ้นต้นด้วย '>' ออก
    - เอาช่องว่างและตัวอักษรที่ไม่ใช่ ACGT ออก
    """
    lines = text.splitlines()
    parts = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            continue
        parts.append(line)
    seq = "".join(parts).upper()
    seq = "".join(ch for ch in seq if ch in "ACGT")
    return seq

def percent_identity_poswise(s1, s2):
    """
    เปรียบเทียบตำแหน่งต่อ ตำแหน่ง: นับตำแหน่งที่ตรงกันในช่วงความยาวขั้นต่ำของ s1,s2
    คืนค่า (pid_fraction, matches, compared_length, len_s1, len_s2)
    - pid_fraction = matches / compared_length (หรือ None ถ้าไม่สามารถเปรียบเทียบได้)
    """
    n1, n2 = len(s1), len(s2)
    if n1 == 0 or n2 == 0:
        return None, 0, 0, n1, n2
    L = min(n1, n2)
    matches = sum(1 for i in range(L) if s1[i] == s2[i])
    pid = matches / L
    return pid, matches, L, n1, n2

def kmer_jaccard(s1, s2, k):
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

# --- ตัวอย่างสั้น ๆ (เพื่อทดสอบ) ---
EXAMPLES = {
    "Human": "ATGCTAGCTAGCTACGATCGATCGTAGCTAG",
    "Chimpanzee": "ATGCTAGCTAGCTACGATCGATCGTAGCTAG",
    "Bonobo": "ATGCTAGCTAGCTACGATCGATCGTGGCTAG",
    "Gorilla": "ATGCTAGCTAGCTACGATCGATCGTGGCTAA",
    "Orangutan": "ATGCTAGCTAGCTACGATCAATCGTGGCTAG",
}

# --- Sidebar: ตัวเลือก ---
st.sidebar.header("Options")
use_examples = st.sidebar.checkbox("ใช้ตัวอย่างสั้น ๆ (demo)", value=True)
k = st.sidebar.slider("k (k-mer) for Jaccard", min_value=3, max_value=9, value=5)

# --- ช่องให้ผู้ใช้กรอก / วาง sequence ---
st.header("วางหรือพิมพ์ลำดับ DNA (FASTA หรือ raw)")
species = ["Human", "Chimpanzee", "Bonobo", "Gorilla", "Orangutan"]
seqs = {}
for sp in species:
    default = EXAMPLES[sp] if use_examples else ""
    txt = st.text_area(f"{sp}", value=default, height=120, key=sp)
    seqs[sp] = clean_seq(txt)

# --- ตรวจสอบ / คำนวณผล ---
human = seqs["Human"]
if len(human) == 0:
    st.warning("กรุณาใส่ลำดับของ Human ก่อน (สามารถใช้ตัวอย่างเพื่อทดลองได้)")
else:
    st.subheader("ผลการเปรียบเทียบ (Human เทียบกับแต่ละสายพันธุ์)")
    rows = []
    for sp in species[1:]:
        other = seqs[sp]
        pid, matches, compared_len, n1, n2 = percent_identity_poswise(human, other)
        if pid is None:
            pid_str = "ไม่สามารถคำนวณ (ลำดับว่าง)"
        else:
            pid_str = f"{pid*100:.2f}%  ({matches}/{compared_len} ตำแหน่งเทียบได้)"
        jk = kmer_jaccard(human, other, k)
        jk_str = f"{jk*100:.2f}%" if jk is not None else f"ไม่พอความยาว (ต้อง >= {k})"
        rows.append({
            "Species": sp,
            "Len (Human)": len(human),
            "Len (other)": len(other),
            "Pos-wise %ID": pid_str,
            f"k-mer Jaccard (k={k})": jk_str
        })
    st.table(rows)

    # เพิ่มคำอธิบายสั้น ๆ
    st.markdown("""
**คำเตือน / ข้อจำกัด (สำคัญ):**
- Percent identity แบบตำแหน่งต่อ-ตำแหน่งในที่นี้เป็นวิธีง่าย ๆ — มัน **ไม่ใช่** global alignment จริง (เช่น Needleman–Wunsch) และจะนับเฉพาะตำแหน่งที่เปรียบเทียบกันได้ (ช่วงความยาวต่ำสุด)
- k-mer Jaccard เป็นวิธี alignment-free ที่วัดความเหมือนของ **ชุด** ก้อนตัวอักษรยาว k; ได้ประโยชน์เมื่อต้องการวัดความคล้ายโดยไม่ต้องจัดเรียง
- สำหรับลำดับจีโนมยาว ๆ หรือการวิเคราะห์เชิงลึก ควรใช้เครื่องมือเฉพาะทาง (Biopython, mafft, bwa, minimap2 ฯลฯ)
""")
