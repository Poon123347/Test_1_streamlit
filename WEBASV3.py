import streamlit as st
from Bio import Entrez, SeqIO, pairwise2
from io import StringIO
import pandas as pd
Entrez.email = "poonthakorn@gmail.com"

species_accessions = {
    "Human": "NC_012920.1",
    "Chimpanzee": "NC_001643.1",
    "Gorilla": "NC_001645.1",
    "Orangutan": "NC_002083.1"
    # เพิ่มสายพันธุ์ได้ตามต้องการ
}

# ฟังก์ชันดึงลำดับ ND5 gene จาก accession number
def extract_nd5(accession_id):
    with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
        # วนลูปหา feature ประเภท CDS ที่ชื่อ gene เป็น nd5
        for feature in record.features:
            if feature.type == "CDS" and "gene" in feature.qualifiers:
                if "nd5" in feature.qualifiers["gene"][0].lower():
                    # คืนค่า sequence ของ ND5 gene ที่เจอ
                    return feature.extract(record.seq)
    return None  # ถ้าไม่เจอ gene nd5 ให้คืนค่า None

def to_fasta(header, seq) -> str:
    s = str(seq)
    buf = StringIO()
    buf.write(f">{header}\n")
    for i in range(0, len(s), 60):
        buf.write(s[i:i+60] + "\n")
    return buf.getvalue()

def percent_identity(seq1, seq2):
    # ทำ global alignment แบบง่าย ๆ (score นับตำแหน่งที่ตรงกัน)
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    aln1, aln2, score, start, end = alignments[0]
    # นับจำนวนเบสที่ตรงกันใน alignment
    matches = sum(base1 == base2 for base1, base2 in zip(aln1, aln2))
    # คืนค่าเป็นเปอร์เซ็นต์เทียบกับความยาวที่ยาวที่สุดของสองลำดับ
    return matches / max(len(seq1), len(seq2)) * 100

st.set_page_config(page_title="Simple DNA Compare", layout="wide")
st.title("เปรียบเทียบความเหมือนของ DNA ระหว่างมนุษย์และลิง 4 สายพันธุ์")
st.markdown("แอปจะโหลดลำดับ DNA จากเซิร์ฟเวอร์ NCBI ในโฟลเดอร์เดียวกันโดยอัตโนมัติ แล้วคำนวณ %Identity")
st.caption("กดปุ่มด้านล่างเพื่อดึง ND5 ของทุกสายพันธุ์ แล้วดาวน์โหลดเป็น FASTA / ดูลำดับ")

if "nd5_seqs" not in st.session_state:
    st.session_state.nd5_seqs = {}

if st.button("ดึง ND5 และสร้างปุ่มดาวน์โหลด"):
    st.markdown("ปุ่มกดสำเร็จ")
    nd5_seqs = {}
    for sp, acc in species_accessions.items():
        with st.spinner(f"กำลังดึง {sp} ({acc}) ..."):
            seq = extract_nd5(acc)
            if seq:
                nd5_seqs[sp] = seq
    st.session_state.nd5_seqs = nd5_seqs
    print("Download",nd5_seqs)

st.divider()
st.subheader("ไฟล์ดาวน์โหลด & แสดง DNA ที่ดึงมา")
nd5_seqs = st.session_state.nd5_seqs
print("after download",nd5_seqs)
if nd5_seqs:
    for sp, seq in nd5_seqs.items():
        fasta_text = to_fasta(f"{sp}_ND5", seq)

        c1, c2 = st.columns([1, 2])
        with c1:
            st.download_button(
                label=f"⬇️ ดาวน์โหลด {sp} ND5 (FASTA)",
                data=fasta_text,
                file_name=f"{sp}_ND5.fasta",
                mime="text/plain",
                use_container_width=True
            )
        with c2:
            with st.expander(f"แสดงลำดับ {sp} ND5 ({len(seq)} bp)"):
                st.code(str(seq), language="text")
    st.session_state.nd5_seqs = nd5_seqs
else:
    st.info("ยังไม่มีลำดับ ND5 ที่ดึงมาได้")

st.divider()
if st.button("✅ CHECK similarity (%Identity)"):
    if len(nd5_seqs) < 2:
        st.warning("ต้องมีลำดับอย่างน้อย 2 สายพันธุ์เพื่อคำนวณความเหมือน")
    else:
        # สร้างตารางจับคู่ทุกสายพันธุ์
        species = list(nd5_seqs.keys())
        matrix = []
        for i, sp_i in enumerate(species):
            row = []
            for j, sp_j in enumerate(species):
                if i == j:
                    row.append(100.0)
                else:
                    pid = percent_identity(nd5_seqs[sp_i], nd5_seqs[sp_j])
                    row.append(round(pid, 2))
                print(row)
            matrix.append(row)
            print(matrix)
        df = pd.DataFrame(matrix, index=species, columns=species)
        st.subheader("ตาราง %Identity (ND5 CDS)")
        st.dataframe(df, use_container_width=True)