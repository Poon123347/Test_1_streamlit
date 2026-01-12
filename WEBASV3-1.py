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
# NCBI CONFIG
# ===============================
Entrez.email = "poonthakorn@gmail.com"

# ===============================
# CONSTANT: ND5 LENGTH
# ===============================
ND5_LENGTH = 1812

# ===============================
# SPECIES ACCESSIONS
# ===============================
species_accessions = {
    "Human": "NC_012920.1",
    "Chimpanzee": "NC_001643.1",
    "Bonobo": "NC_001644.1",
    "Western Gorilla": "NC_001645.1",
    "Eastern Gorilla": "NC_011120.1",
    "Sumatran Orangutan": "NC_002083.1",
    "Bornean Orangutan": "NC_001646.1",
    "Common Gibbon": "NC_002082.1",
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
    "Black Colobus": "NC_008219.1",
    "Pig": "NC_000845.1",
    "Cow": "NC_001567.1",
    "Horse": "NC_001640.1",
    "Dog": "NC_002008.4",
    "Cat": "NC_001700.1",
    "Sheep": "NC_001941.1",
    "Goat": "NC_005044.2",
    "Mouse": "NC_005089.1",
    "Rat": "NC_001665.2",
    "Elephant": "NC_003319.1",
    "Chicken": "NC_001323.1",
    "Duck": "NC_009684.1",
    "Zebrafish": "NC_002333.2"
}

# ===============================
# LANGUAGE SUPPORT
# ===============================
LANG = {
    "English": {
        "title": "🧬 ND5 Phylogenetic Explorer",
        "species_select": "Pleas select the following species (min 2)",
        "fetch_nd5": "🔍 Fetch ND5",
        "identity": "📊 Identity Matrix",
        "build_tree": "🌳 Build Phylogenetic Tree",
        "no_data": "No ND5 data yet"
    },
    "ภาษาไทย": {
        "title": "🧬 ตัวสำรวจสายวิวัฒนาการ ND5",
        "species_select": "กรุณาเลือกสายพันธุ์ตัวอย่าง (อย่างน้อย 2 ชนิด)",
        "fetch_nd5": "🔍 ดึง ND5",
        "identity": "📊 ตาราง % ความเหมือน",
        "build_tree": "🌳 สร้างแผนภาพต้นไม้",
        "no_data": "ยังไม่มีข้อมูล ND5"
    }
}

# ===============================
# SIDEBAR
# ===============================
language = st.sidebar.selectbox("🌐 Language / ภาษา", ["English", "ภาษาไทย"])
T = LANG[language]

st.title(T["title"])

# ===============================
# SPECIES SELECTOR
# ===============================
selected_species = st.multiselect(
    T["species_select"],
    options=list(species_accessions.keys())
)

selected_accessions = {
    name: species_accessions[name]
    for name in selected_species
}

# ===============================
# ND5 EXTRACTION (HARD CUT TO 1812 bp)
# ===============================
def extract_nd5(accession_id):
    with Entrez.efetch(
        db="nucleotide",
        id=accession_id,
        rettype="gb",
        retmode="text"
    ) as handle:
        record = SeqIO.read(handle, "genbank")

    for feat in record.features:
        if feat.type == "CDS" and "gene" in feat.qualifiers:
            if feat.qualifiers["gene"][0].lower() == "nd5":
                seq = feat.extract(record.seq)
                return seq[:ND5_LENGTH] if len(seq) >= ND5_LENGTH else None

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
    return score / ND5_LENGTH * 100

# ===============================
# FETCH ND5
# ===============================
if st.button(T["fetch_nd5"]):
    if len(selected_accessions) < 2:
        st.warning(T["species_select"])
    else:
        nd5_seqs = {}
        with st.spinner("Fetching ND5 (1812 bp)..."):
            for name, acc in selected_accessions.items():
                seq = extract_nd5(acc)
                if seq:
                    nd5_seqs[name] = seq

        st.session_state["nd5_seqs"] = nd5_seqs

# ===============================
# SHOW ND5 LENGTH
# ===============================
if "nd5_seqs" in st.session_state:
    st.subheader("ND5 length (bp)")
    for name, seq in st.session_state["nd5_seqs"].items():
        st.text(f"{name} → {len(seq)} bp")

# ===============================
# IDENTITY MATRIX
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
                row.append(100.0 if i == j else round(percent_identity(seqs[a], seqs[b]), 2))
            matrix.append(row)

        df = pd.DataFrame(matrix, index=names, columns=names)
        st.dataframe(df)

# ===============================
# BUILD PHYLOGENETIC TREE
# ===============================
if st.button(T["build_tree"]):
    seqs = st.session_state.get("nd5_seqs", {})

    if not seqs:
        st.warning(T["no_data"])
    else:
        # ===== FASTA =====
        fasta_io = StringIO()
        for name, seq in seqs.items():
            fasta_io.write(f">{name}\n{seq}\n")
        fasta_io.seek(0)

        # ===== Alignment (แบบง่าย) =====
        alignment = AlignIO.read(fasta_io, "fasta")

        # ===== Distance + Tree =====
        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(alignment)

        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        tree.ladderize()

        # ซ่อน Inner node
        for clade in tree.find_clades():
            if clade.name and clade.name.startswith("Inner"):
                clade.name = None

        # ===== Plot =====
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111)

        Phylo.draw(tree, axes=ax, do_show=False)

        st.pyplot(fig)



# ===============================
# DETAILED vs SIMPLE DNA MECHANISMS (2 LANGUAGES)
# ===============================

if language == "English":
    # ===== Detailed Zone =====
    st.markdown(
        """
        ## 🔬 Molecular Mechanisms of DNA (Detailed)

        **DNA Packaging and Chromatin Structure**  
        Eukaryotic DNA is organized into chromatin, with nucleosomes as the fundamental 
        unit, where ~146–147 base pairs of DNA wrap around an octamer of histone proteins.  
        Chromatin exists in forms such as **euchromatin** (loosely packed, transcriptionally 
        active) and **heterochromatin** (tightly packed, transcriptionally repressed). 

        **Histone Modifications and Transcription Regulation**  
        Histone tails can undergo **post-translational modifications** (e.g., acetylation, 
        methylation). Acetylation (such as H3K27ac) is associated with active chromatin, 
        while methylation can either activate or repress depending on context. 

        **Chromatin Remodeling**  
        Chromatin structure is modulated by ATP-dependent remodelers and covalent histone 
        modifiers, allowing nucleosomal dynamics that either expose or occlude DNA from 
        transcription factors and RNA polymerase. 

        **Transcription and Cis-Regulatory Elements**  
        Gene transcription begins at promoters via RNA polymerase II and transcription 
        factors. Enhancers increase transcription from distances, while silencers 
        decrease or repress transcription, collectively shaping gene expression.

        **Epigenetics**  
        Epigenetics refers to inheritable changes in gene activity without altering DNA 
        sequence, typically through mechanisms like DNA methylation and histone modifications 
        that change chromatin accessibility.

        **Non-coding DNA and RNA**  
        Much of the genome does not code for protein but contains regulatory sequences and 
        non-coding RNAs (e.g., microRNAs, lncRNAs), which regulate gene expression at 
        transcriptional or post-transcriptional levels.
        """
    )

    # ===== Simple Zone =====
    st.markdown(
        """
        ## 🧬 Molecular Mechanisms of DNA (Simplified)

        Imagine DNA as a very long instruction manual. To store it in a small space, 
        cells wind DNA around proteins called histones — like winding yarn around spools.
        Tightly wound regions (heterochromatin) are hard to read, while loosely wound 
        regions (euchromatin) are easier to read.

        Some chemical markers act like sticky notes:
        - **Acetyl marks** help open the DNA manual for easier reading.
        - **Methyl marks** can act like brakes or locks, sometimes closing parts of 
          the manual.

        Long-distance switches like **enhancers** act like accelerators, making “reading” 
        faster, while **silencers** act like brakes.

        After DNA is “read” into RNA, some RNA sequences (like microRNAs) can act like 
        editors, deciding whether the instructions are used or thrown away.
        """
    )

elif language == "ภาษาไทย":
    # ===== Detailed Zone =====
    st.markdown(
        """
        ## 🔬 กลไกระดับโมเลกุลของ DNA (ระดับลึก)

        **การจัดบรรจุ DNA และโครงสร้างโครมาติน**  
        DNA ของยูคาริโอตจะถูกจัดเป็นโครมาติน โดยหน่วยพื้นฐานคือ **นิวคลีโอโซม** ซึ่ง DNA 
        ประมาณ 146–147 เบสพันรอบฮิสโตน โปรตีนทั้ง 8 ตัว จากนั้นมีการจัดเรียงเป็นโครมาตินแบบ 
        **ยูโครมาติน** (หลวม ถอดรหัสได้) และ **เฮเทอโรโครมาติน** (แน่น ถอดรหัสยาก)

        **การปรับแต่งฮิสโตนและการควบคุมการถอดรหัส**  
        หางฮิสโตนสามารถถูกปรับแต่งทางเคมี เช่น **อะซิทิล (acetylation)** หรือ **เมทิล (methylation)** 
        ซึ่งส่งผลต่อการเข้าถึง DNA และการถอดรหัส โดยบางป้าย (เช่น H3K27ac) มักสัมพันธ์กับ 
        การเปิดการถอดรหัส

        **โครมาตินรีโมเดลลิง**  
        กลไก dynamic ของโครมาตินเกิดจากเอนไซม์ที่ใช้ ATP และตัวดัดแปลงฮิสโตน ที่สามารถเลื่อน, 
        ถอด หรือเปลี่ยนโครงสร้างนิวคลีโอโซม เพื่อเปิดหรือปิดบริเวณ DNA สำหรับยีนต่าง ๆ

        **การถอดรหัสและ cis-regulatory elements**  
        การถอดรหัสเริ่มที่โปรโมเตอร์ ผ่าน RNA polymerase II และโปรตีน transcription factors 
        enhancer จะช่วยเร่งการถอดรหัสจากระยะไกล ในขณะที่ silencer จะชะลอ/ยับยั้งการถอดรหัส

        **เอพิเจเนติกส์ (Epigenetics)**  
        Epigenetics คือการเปลี่ยนแปลงที่สืบทอดได้ของระดับการแสดงออกของยีนโดยไม่เปลี่ยน 
        ลำดับ DNA เช่น การเมทิลของ DNA หรือการปรับแต่งฮิสโตน ซึ่งเปลี่ยนการเข้าถึงโครมาติน

        **DNA/RNA ที่ไม่เข้ารหัส**  
        ส่วนใหญ่ของจีโนมไม่เข้ารหัสโปรตีน แต่มียีน regulatory sequences และ RNA ที่ไม่เข้ารหัส 
        เช่น microRNA หรือ lncRNA ซึ่งช่วยควบคุมการถอดรหัสและระดับหลังการถอดรหัส
        """
    )

    # ===== Simple Zone =====
    st.markdown(
        """
        ## 🧬 กลไกระดับโมเลกุลของ DNA (ในภาษาง่าย)

        ลองนึกว่า DNA เป็นคู่มือคำสั่งชีวิตที่ยาวมาก เพื่อจัดเก็บมัน เราจะพัน DNA รอบโปรตีน 
        ที่เรียกว่า **ฮิสโตน** เหมือนพันด้ายรอบม้วน การพันแน่น (heterochromatin) เหมือนการ 
        ม้วนหนังสือแน่น ๆ ทำให้อ่านยาก แต่การพันหลวม (euchromatin) เหมือนเปิดหนังสืออ่านง่าย

        เครื่องหมายทางเคมีบางอย่างเปรียบเสมือนกระดาษโน้ต:
        - **อะซิทิล** ช่วยให้ DNA เปิดอ่านได้ง่ายขึ้น
        - **เมทิล** บางครั้งเหมือนล็อกหรือเบรก ทำให้ยากต่อการอ่าน

        ตัวควบคุมระยะไกลอย่าง **enhancer** ทำงานเหมือนสวิตช์เร่งให้การอ่านเร็ว, 
        ส่วน **silencer** เหมือนเบรกที่ชะลอการอ่าน

        หลังจาก DNA ถูกถอดเป็น RNA แล้ว RNA บางชนิด เช่น microRNA ก็เหมือน 
        บรรณาธิการที่ช่วยตัดสินใจว่าจะใช้ข้อความต่อไปหรือไม่
        """
    )
