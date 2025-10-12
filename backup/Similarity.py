from Bio import Entrez, SeqIO, pairwise2

# ตั้งค่าอีเมลสำหรับติดต่อ NCBI ตามนโยบายของเขา
Entrez.email = "poonthakorn@gmail.com"

# รายชื่อสายพันธุ์พร้อม accession ของ mitochondrial genome ที่จะดึง ND5 gene
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

# โหลดลำดับ ND5 ของแต่ละสายพันธุ์มาเก็บใน dictionary
nd5_seqs = {}
for sp, acc in species_accessions.items():
    seq = extract_nd5(acc)
    if seq:
        nd5_seqs[sp] = seq

# ฟังก์ชันคำนวณเปอร์เซ็นต์ความเหมือนกัน (percent identity) ระหว่าง 2 sequences
def percent_identity(seq1, seq2):
    # ทำ global alignment แบบง่าย ๆ (score นับตำแหน่งที่ตรงกัน)
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    aln1, aln2, score, start, end = alignments[0]
    # นับจำนวนเบสที่ตรงกันใน alignment
    matches = sum(base1 == base2 for base1, base2 in zip(aln1, aln2))
    # คืนค่าเป็นเปอร์เซ็นต์เทียบกับความยาวที่ยาวที่สุดของสองลำดับ
    return matches / max(len(seq1), len(seq2)) * 100

# สร้างรายการสายพันธุ์เพื่อจัดลำดับในการแสดงผล
species = list(nd5_seqs.keys())

# พิมพ์ header ตารางเว้นช่องว่างมุมบนซ้าย 12 ช่อง (เพื่อจัดรูปแบบ)
print(f"{'Species':19}", end='')
for sp in species:
     print(f"{sp:12}", end='')
print()  # ขึ้นบรรทัดใหม่

# พิมพ์ตารางเปอร์เซ็นต์ความเหมือนกันของ ND5 gene ระหว่างแต่ละคู่สายพันธุ์
for sp1 in species:
    
    print(f"{sp1:12}", end='')  # ชื่อสายพันธุ์แถวแรกของแต่ละแถว
    for sp2 in species:
        # print(f"{sp2:12}", end='')
        pid = percent_identity(nd5_seqs[sp1], nd5_seqs[sp2])
        print(f"{pid:11.2f}%", end=' ')  # ค่าเปอร์เซ็นต์จัดช่องว่าง 11 ตัวอักษร ทศนิยม 2 ตำแหน่ง
    print()  # ขึ้นบรรทัดใหม่สำหรับสายพันธุ์ถัดไป