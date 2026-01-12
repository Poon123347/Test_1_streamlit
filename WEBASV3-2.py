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

elif language == "English":
    st.markdown("""
    ## 🔬 Molecular Mechanisms of DNA (Comprehensive)

    **DNA Packaging and Chromatin Structure**  
    Eukaryotic DNA is packaged into chromatin with nucleosomes as the basic units. Chromatin may be 
    loosely packed (**euchromatin**) or densely packed (**heterochromatin**) which affects accessibility 
    for transcription. :contentReference[oaicite:9]{index=9}

    **Histone Modifications**  
    Post-translational modifications like acetylation (e.g., **H3K27ac**) are associated with active 
    transcription and enhancer regions, while methylation marks vary with context and can activate or 
    repress transcription. :contentReference[oaicite:10]{index=10}

    **DNA Methylation (CpG Islands)**  
    DNA methylation mainly occurs at CpG dinucleotides and hypermethylation at promoter CpG islands 
    can block transcription factor binding and silence gene expression. :contentReference[oaicite:11]{index=11}

    **Chromatin Remodeling**  
    ATP-dependent remodelers reposition or remove nucleosomes to regulate access to DNA.

    **Transcription and Cis-Regulatory Elements**  
    Transcription starts at promoters with RNA polymerase II and transcription factors. Enhancers 
    boost transcription, while silencers repress it.

    **Epigenetics**  
    Epigenetics refers to heritable changes in gene expression without altering DNA sequence. DNA 
    methylation and histone modifications define cell-specific expression profiles. :contentReference[oaicite:12]{index=12}

    **Non-coding DNA and RNA**  
    Non-coding RNAs including **microRNAs** and **lncRNAs** regulate gene expression post-transcriptionally 
    by inhibiting translation or promoting degradation. :contentReference[oaicite:13]{index=13}

    --- 

    ## 🧬 Human DNA Analysis with Bioinformatics
    **Bioinformatics** uses algorithms and computational tools to analyze DNA sequences for function, 
    evolution, and disease association.

    **Study Objectives**  
    1. Analyze human DNA sequences to identify genetic variation.  
    2. Construct **phylogenetic trees** to infer evolutionary relationships.  
    3. Identify mutations related to disease and gene function.

    **Expected Benefits**  
    1. Understanding human lineage and genetic relationships.  
    2. Insight into origins and dispersal of humans.  
    3. Supporting genomic research in medicine.

    **Bioinformatics Pipeline**  
    1. Sample collection & sequence QC  
    2. Sequence assembly  
    3. Multiple sequence alignment  
    4. Evolutionary distance calculation  
    5. Phylogenetic tree construction

    --- 

    ## 🧬 ND5 Gene Function  
    The MT-ND5 gene encodes a subunit of mitochondrial Complex I, which plays a key role in ATP production. :contentReference[oaicite:14]{index=14}

    ## 🚨 Clinical Significance  
    Mutations in ND5 are associated with mitochondrial diseases such as MELAS and Leigh syndrome that impair energy production. :contentReference[oaicite:15]{index=15}
    """)


elif language == "ภาษาไทย":
    # ===== Detailed Zone =====
    st.markdown("""
    ## 🔬 กลไกระดับโมเลกุลของ DNA (ระดับลึกและสรุปครบถ้วน)

    **การจัดบรรจุ DNA และโครมาติน**  
    DNA ของยูคาริโอตถูกบรรจุเป็นโครมาติน โดยหน่วยพื้นฐานคือ **นิวคลีโอโซม** ซึ่ง DNA ประมาณ 146–147 เบสพันรอบฮิสโตน โปรตีนทั้ง 8 ตัว โครมาตินมีสองสภาพหลักคือ **ยูโครมาติน** (หลวม สามารถเข้าถึงได้และรองรับการถอดรหัส) และ **เฮเทอโรโครมาติน** (แน่น ยับยั้งการถอดรหัส). :contentReference[oaicite:1]{index=1}

    **การปรับแต่งฮิสโตน (Histone Modifications)**  
    หางฮิสโตนสามารถถูกปรับแต่ง เช่น **อะซิทิล (acetylation)** และ **เมทิล (methylation)** ซึ่งส่งผลต่อการเข้าถึง DNA และการถอดรหัส  
    - การอะซิทิล เช่น **H3K27ac** มักสัมพันธ์กับการเปิดยีนและ enhancer activity. :contentReference[oaicite:2]{index=2}  
    - การเมทิลที่บางตำแหน่งอาจเปิดหรือปิดยีน ขึ้นอยู่กับบริบท (เช่น H3K4, H3K9, H3K27). :contentReference[oaicite:3]{index=3}

    **การเมทิลของ DNA (CpG islands)**  
    การเมทิลของไซโตซีนที่ CpG dinucleotides โดยเฉพาะใน CpG islands สามารถยับยั้งการจับของ transcription factors และก่อให้เกิดการปิดการแสดงออกของยีน. :contentReference[oaicite:4]{index=4}

    **โครมาตินรีโมเดลลิง (Chromatin Remodeling)**  
    โครมาตินสามารถเปลี่ยนโครงสร้างโดยโปรตีนที่ใช้พลังงาน ATP เพื่อเปิดหรือปิดการเข้าถึง DNA สำหรับเครื่องจักรถอดรหัส.

    **การถอดรหัสและ cis-regulatory elements**  
    การถอดรหัสเริ่มที่โปรโมเตอร์ โดย RNA polymerase II และ transcription factors. **Enhancer** ส่งเสริมการถอดรหัส ขณะที่ **silencer** ยับยั้งการแสดงออกของยีน.

    **เอพิเจเนติกส์ (Epigenetics)**  
    เอพิเจเนติกส์คือการควบคุมการแสดงออกของยีนโดยไม่เปลี่ยนลำดับเบสของ DNA ช่วยกำหนดโปรไฟล์การแสดงออกเฉพาะเซลล์ ซึ่งรวมทั้ง DNA methylation และ histone modifications. :contentReference[oaicite:5]{index=5}

    **DNA และ RNA ที่ไม่เข้ารหัส**  
    แม้ DNA ส่วนใหญ่ของมนุษย์จะไม่เข้ารหัสโปรตีน แต่ยังมีบทบาทควบคุม เช่น **microRNA** และ **lncRNA** ที่ยับยั้งการแปลหรือทำให้ mRNA ถูกย่อยสลาย. :contentReference[oaicite:6]{index=6}

    --- 

    ## 🧬 การถอดรหัสดีเอ็นเอของมนุษย์ด้วยชีวสารสนเทศ (Bioinformatics)
    **ชีวสารสนเทศ (Bioinformatics)** คือการใช้คอมพิวเตอร์และอัลกอริทึมในการวิเคราะห์ข้อมูลลำดับ DNA เพื่อศึกษาเชื้อชาติ หน้าที่ยีน วิวัฒนาการ และการกลายพันธุ์.

    **วัตถุประสงค์ของการศึกษา**  
    1. วิเคราะห์ลำดับ DNA ของมนุษย์เพื่อศึกษาโครงสร้างและความแตกต่างเชิงพันธุกรรม  
    2. สร้าง **แผนภูมิวิวัฒนาการ (phylogenetic tree)** เพื่อทำความเข้าใจความสัมพันธ์เชื้อชาติ  
    3. ระบุการกลายพันธุ์ที่เกี่ยวข้องกับโรคและฟังก์ชันยีน

    **ประโยชน์ที่คาดว่าจะได้รับ**  
    1. เข้าใจสายพันธุ์และความสัมพันธ์ของมนุษย์จากข้อมูล DNA  
    2. เข้าใจต้นกำเนิดและการกระจายของมนุษย์  
    3. สนับสนุนการวิจัยทางการแพทย์ด้วยข้อมูล genomic

    ## 🛠️ วิธีการวิเคราะห์ (Bioinformatics Pipeline)
    1. เก็บตัวอย่าง & QC ของลำดับ  
    2. จัดเรียง/ประกอบลำดับ (assembly)  
    3. จัดเรียงหลายลำดับ (alignment)  
    4. คำนวณระยะทางวิวัฒนาการ  
    5. สร้างและตีความ **phylogenetic tree**

    --- 

    ## 🧬 หน้าที่ทางชีววิทยาของยีน ND5
    ยีน MT-ND5 ให้ข้อมูลในการสร้างหน่วยย่อยของ **Complex I** ในไมโตคอนเดรีย ซึ่งมีบทบาทในกระบวนการสร้างพลังงาน (ATP). :contentReference[oaicite:7]{index=7}

    ## 🚨 ความสำคัญทางการแพทย์
    การกลายพันธุ์ใน ND5 เกี่ยวข้องกับโรค mitochondrial เช่น MELAS และ Leigh syndrome ซึ่งส่งผลต่อการผลิตพลังงาน. :contentReference[oaicite:8]{index=8}
    """)
    # ===== Simple Zone =====
    st.markdown("""
    ## 🧬 กลไกระดับโมเลกุลของ DNA (ภาษาเข้าใจง่าย)

    DNA เหมือนคู่มือชีวิตที่บรรจุอย่างแน่นหนาในนิวเคลียส โดยพันรอบโปรตีนฮิสโตน
    เมื่อ DNA เปิดจะอ่านง่าย แต่พันแน่นจะอ่านยาก เช่น enhancer ช่วยเร่ง ถ้า methylation มาเยอะ
    การเปิดอ่านก็ยากขึ้น และ RNA บางชนิดอย่าง microRNA ก็ช่วยควบคุมอีกที.
    """)

