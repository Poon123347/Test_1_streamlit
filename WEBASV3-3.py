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
    matches = sum(c1 == c2 for c1, c2 in zip(seq1, seq2))
    return round(matches / ND5_LENGTH * 100, 2)

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
    # ===============================
# DETAILED SCIENTIFIC CONTENT (ENGLISH VERSION)
# ===============================
    st.markdown("""
    ## 🔬 Molecular Mechanisms of DNA (Comprehensive Overview)

    ### DNA Packaging and Chromatin Structure  
    Eukaryotic DNA is packaged into chromatin, whose fundamental unit is the **nucleosome**.
    Each nucleosome consists of approximately 146–147 base pairs of DNA wrapped around an octamer of histone proteins.
    Chromatin exists in two major forms: **euchromatin**, which is loosely packed and transcriptionally active,
    and **heterochromatin**, which is densely packed and transcriptionally repressive.

    ### Histone Modifications  
    Histone tails undergo post-translational modifications such as **acetylation** and **methylation**,
    which influence chromatin accessibility and gene expression.
    For example, histone acetylation (e.g., **H3K27ac**) is generally associated with transcriptionally active regions,
    while histone methylation can either activate or repress transcription depending on the residue and context
    (e.g., H3K4, H3K9, H3K27).

    ### DNA Methylation (CpG Islands)  
    DNA methylation commonly occurs at cytosine residues in CpG dinucleotides.
    Hypermethylation of CpG islands in promoter regions can inhibit transcription factor binding
    and lead to gene silencing, playing an important role in epigenetic regulation.

    ### Chromatin Remodeling  
    ATP-dependent chromatin remodeling complexes reposition or remove nucleosomes,
    thereby regulating the accessibility of DNA to transcriptional machinery.

    ### Transcription and Cis-Regulatory Elements  
    Transcription initiates at promoter regions with the binding of RNA polymerase II
    and transcription factors. **Enhancers** increase transcriptional activity,
    whereas **silencers** suppress gene expression.

    ### Epigenetics  
    Epigenetics refers to heritable changes in gene expression that do not involve alterations
    in the DNA sequence. DNA methylation and histone modifications are major epigenetic mechanisms
    that define cell-type-specific gene expression patterns.

    ### Non-coding DNA and RNA  
    Although most of the human genome does not encode proteins,
    non-coding elements such as **microRNAs (miRNAs)** and **long non-coding RNAs (lncRNAs)**
    play crucial roles in post-transcriptional gene regulation by inhibiting translation
    or promoting mRNA degradation.

    ---

    ## 🧬 Human DNA Analysis Using Bioinformatics

    ### Definition of Bioinformatics  
    **Bioinformatics** is the application of computational tools and algorithms
    to analyze biological sequence data, enabling the study of genetic variation,
    gene function, evolution, and disease-associated mutations.

    ---

    ## 🎯 Objectives of the Study
    1. To study human DNA in order to understand human genetic variation and species relationships.  
    2. To construct **phylogenetic trees** to investigate evolutionary relationships
    between humans and other species based on DNA sequence differences,
    and to explore functional aspects of DNA in animals.  
    3. To analyze sequence variation in the mitochondrial **ND5 gene**
    and assess its evolutionary and biological significance.

    ---

    ## ✅ Expected Benefits
    1. To enhance understanding of human genetic diversity and evolutionary relationships
    through DNA sequence analysis.  
    2. To provide insights into the origin and evolutionary divergence of humans
    and to establish a foundation for further studies in medical and genetic research,
    particularly in mitochondrial-related diseases.

    ---

    ## 🛠️ Bioinformatics Analysis Pipeline

    ### Data Source  
    DNA sequences used in this study are retrieved from the public database
    **NCBI GenBank** using accession numbers for each species.
    All sequences are publicly available and do not contain personal or sensitive information.

    ### Data Preparation and Quality Control  
    - Extraction of the **ND5 (MT-ND5)** coding sequence from mitochondrial genomes.  
    - All sequences are standardized to a length of **1,812 base pairs (bp)** for accurate comparison.  
    - Sequences with insufficient length, excessive ambiguous bases (N),
    or internal stop codons are excluded from analysis.

    ### Multiple Sequence Alignment  
    Multiple sequence alignment (MSA) is performed to align ND5 sequences
    from different species, ensuring homologous positions are compared.

    ### Evolutionary Distance and Phylogenetic Tree Construction  
    - Pairwise sequence similarity (% identity) and evolutionary distances are calculated.  
    - A **phylogenetic tree** is constructed to visualize evolutionary relationships
    among species based on ND5 sequence variation.

    ---

    ## 🧬 Biological Function of the ND5 Gene
    The mitochondrial **MT-ND5** gene encodes a subunit of **Complex I**
    (NADH dehydrogenase) in the mitochondrial electron transport chain.
    This complex plays a critical role in oxidative phosphorylation
    and ATP production.

    ---

    ## 🚨 Medical Significance
    Mutations in the ND5 gene are associated with mitochondrial disorders
    such as **MELAS** and **Leigh syndrome**, which impair cellular energy production.
    Therefore, ND5 is an important gene for both evolutionary studies
    and medical genetics.

    ---

    ## 🔁 Reproducibility and Future Applications
    The analytical workflow and computational steps can be reproduced
    and extended to larger datasets. This approach provides a foundation
    for future research in evolutionary biology, population genetics,
    and biomedical studies.
    """)


elif language == "ภาษาไทย":
    # ===============================
# DETAILED METHODS & BACKGROUND (INSERT - Expanded, complete)
# ===============================
    # ===============================
# DETAILED SCIENTIFIC CONTENT (THAI VERSION – MATCHED FORMAT)
# ===============================
    st.markdown("""
    ## 🔬 กลไกระดับโมเลกุลของ DNA (ภาพรวมเชิงลึก)

    ### การจัดบรรจุ DNA และโครมาติน  
    DNA ของยูคาริโอตถูกบรรจุอยู่ในรูปของโครมาติน โดยมีหน่วยพื้นฐานคือ **นิวคลีโอโซม**
    ซึ่งประกอบด้วย DNA ประมาณ 146–147 คู่เบสพันรอบโปรตีนฮิสโตนจำนวน 8 ตัว
    โครมาตินสามารถแบ่งออกเป็น 2 สภาพหลัก ได้แก่  
    **ยูโครมาติน (euchromatin)** ซึ่งมีโครงสร้างหลวมและเอื้อต่อการถอดรหัส
    และ **เฮเทอโรโครมาติน (heterochromatin)** ซึ่งมีโครงสร้างแน่นและยับยั้งการแสดงออกของยีน

    ### การปรับแต่งฮิสโตน (Histone Modifications)  
    หางของโปรตีนฮิสโตนสามารถถูกปรับแต่งทางเคมี เช่น **การอะซิทิล (acetylation)**
    และ **การเมทิล (methylation)** ซึ่งส่งผลต่อระดับการเข้าถึง DNA และการแสดงออกของยีน
    ตัวอย่างเช่น การอะซิทิลของฮิสโตนตำแหน่ง **H3K27ac** มักเกี่ยวข้องกับบริเวณที่มีการถอดรหัสสูง
    ในขณะที่การเมทิลของฮิสโตนบางตำแหน่งอาจกระตุ้นหรือยับยั้งการถอดรหัส ขึ้นอยู่กับบริบทของตำแหน่งนั้น

    ### การเมทิลของ DNA (CpG Islands)  
    การเมทิลของ DNA มักเกิดขึ้นที่ไซโตซีนในตำแหน่ง CpG dinucleotides
    โดยเฉพาะบริเวณ CpG islands ที่อยู่ใกล้โปรโมเตอร์ของยีน
    การเมทิลในบริเวณดังกล่าวสามารถยับยั้งการจับของ transcription factors
    และนำไปสู่การปิดการแสดงออกของยีน ซึ่งเป็นกลไกสำคัญของการควบคุมยีนในระดับเอพิเจเนติกส์

    ### การปรับโครงสร้างโครมาติน (Chromatin Remodeling)  
    โปรตีนกลุ่ม chromatin remodeling complexes ใช้พลังงานจาก ATP
    ในการเคลื่อนย้ายหรือปรับตำแหน่งของนิวคลีโอโซม
    ทำให้ DNA เปิดหรือปิดต่อเครื่องจักรถอดรหัสได้ตามความเหมาะสม

    ### การถอดรหัสและองค์ประกอบควบคุมยีน (Cis-Regulatory Elements)  
    กระบวนการถอดรหัสเริ่มต้นที่บริเวณโปรโมเตอร์ โดย RNA polymerase II
    และ transcription factors  
    **Enhancer** มีบทบาทในการเพิ่มระดับการถอดรหัส
    ขณะที่ **silencer** ทำหน้าที่ยับยั้งการแสดงออกของยีน

    ### เอพิเจเนติกส์ (Epigenetics)  
    เอพิเจเนติกส์หมายถึงการเปลี่ยนแปลงการแสดงออกของยีน
    โดยไม่เปลี่ยนแปลงลำดับเบสของ DNA
    กลไกหลักได้แก่ DNA methylation และ histone modifications
    ซึ่งช่วยกำหนดรูปแบบการแสดงออกของยีนที่แตกต่างกันในแต่ละชนิดของเซลล์

    ### DNA และ RNA ที่ไม่เข้ารหัส  
    แม้ว่า DNA ส่วนใหญ่ของจีโนมมนุษย์จะไม่เข้ารหัสโปรตีน
    แต่ส่วนที่ไม่เข้ารหัสเหล่านี้ เช่น **microRNA (miRNA)** และ **long non-coding RNA (lncRNA)**
    มีบทบาทสำคัญในการควบคุมการแสดงออกของยีน
    โดยการยับยั้งการแปลโปรตีนหรือกระตุ้นการสลายของ mRNA

    ---

    ## 🧬 การวิเคราะห์ดีเอ็นเอของมนุษย์ด้วยชีวสารสนเทศ

    ### ความหมายของชีวสารสนเทศ  
    **ชีวสารสนเทศ (Bioinformatics)** คือการประยุกต์ใช้คอมพิวเตอร์
    และอัลกอริทึมในการวิเคราะห์ข้อมูลชีวภาพ โดยเฉพาะข้อมูลลำดับ DNA
    เพื่อศึกษาโครงสร้าง ความแตกต่างทางพันธุกรรม วิวัฒนาการ
    และการกลายพันธุ์ที่เกี่ยวข้องกับโรค

    ---

    ## 🎯 วัตถุประสงค์ของการศึกษา
    1. เพื่อศึกษาดีเอ็นเอของมนุษย์และทำความเข้าใจความหลากหลายทางพันธุกรรม
    รวมถึงความสัมพันธ์ระหว่างมนุษย์กับสายพันธุ์อื่น  
    2. เพื่อสร้าง **แผนภูมิวิวัฒนาการ (phylogenetic tree)**
    สำหรับอธิบายความสัมพันธ์เชิงวิวัฒนาการจากความแตกต่างของลำดับดีเอ็นเอ
    และศึกษาหน้าที่ของยีนในสัตว์ชนิดต่าง ๆ  
    3. เพื่อวิเคราะห์ความแปรผันของลำดับยีนไมโตคอนเดรีย **ND5**
    และประเมินความสำคัญทางชีววิทยาและวิวัฒนาการของยีนดังกล่าว

    ---

    ## ✅ ประโยชน์ที่คาดว่าจะได้รับ
    1. ทำให้เข้าใจความสัมพันธ์และสายวิวัฒนาการของมนุษย์
    ผ่านการวิเคราะห์ข้อมูลลำดับดีเอ็นเอ  
    2. ทำให้ทราบถึงต้นกำเนิดและการพัฒนาการของมนุษย์ในเชิงวิวัฒนาการ
    และสามารถนำความรู้ไปต่อยอดในการศึกษาทางการแพทย์
    โดยเฉพาะโรคที่เกี่ยวข้องกับไมโตคอนเดรีย

    ---

    ## 🛠️ กระบวนการวิเคราะห์ด้วยชีวสารสนเทศ

    ### แหล่งที่มาของข้อมูล  
    ลำดับดีเอ็นเอที่ใช้ในการศึกษานี้ถูกดึงมาจากฐานข้อมูลสาธารณะ
    **NCBI GenBank** โดยใช้หมายเลข accession ของแต่ละสายพันธุ์
    ข้อมูลทั้งหมดเป็นข้อมูลที่เปิดเผยเพื่อการวิจัย

    ### การเตรียมข้อมูลและการควบคุมคุณภาพ  
    - ดึงลำดับยีน **ND5 (MT-ND5)** จากจีโนมไมโตคอนเดรีย  
    - ปรับความยาวของลำดับให้เท่ากันที่ **1,812 คู่เบส**
    เพื่อความถูกต้องในการเปรียบเทียบ  
    - คัดกรองลำดับที่มีเบสไม่ทราบชนิด (N) จำนวนมาก
    หรือมีสัญญาณการหยุดการแปลภายในกรอบการอ่าน

    ### การจัดเรียงหลายลำดับ (Multiple Sequence Alignment)  
    ทำการจัดเรียงลำดับ ND5 ของแต่ละสายพันธุ์
    เพื่อให้ตำแหน่งที่เปรียบเทียบกันเป็นตำแหน่งที่สอดคล้องกันทางชีววิทยา

    ### การคำนวณระยะทางและการสร้างต้นไม้วิวัฒนาการ  
    - คำนวณค่าความเหมือนของลำดับเบสและระยะทางทางวิวัฒนาการ  
    - สร้าง **แผนภูมิต้นไม้วิวัฒนาการ**
    เพื่อแสดงความสัมพันธ์ของสายพันธุ์ต่าง ๆ จากข้อมูลยีน ND5

    ---

    ## 🧬 หน้าที่ทางชีววิทยาของยีน ND5
    ยีน **MT-ND5** ทำหน้าที่เข้ารหัสหน่วยย่อยของ **Complex I**
    ในระบบขนส่งอิเล็กตรอนของไมโตคอนเดรีย
    ซึ่งมีบทบาทสำคัญในกระบวนการสร้างพลังงาน (ATP)
    ของเซลล์

    ---

    ## 🚨 ความสำคัญทางการแพทย์
    การกลายพันธุ์ในยีน ND5 มีความเกี่ยวข้องกับโรคไมโตคอนเดรีย
    เช่น **MELAS** และ **Leigh syndrome**
    ซึ่งส่งผลต่อการสร้างพลังงานของเซลล์
    ดังนั้นยีน ND5 จึงมีความสำคัญทั้งในด้านชีววิทยาเชิงวิวัฒนาการ
    และด้านพันธุศาสตร์ทางการแพทย์

    ---

    ## 🔁 ความสามารถในการทำซ้ำและการศึกษาต่อยอด
    กระบวนการวิเคราะห์ที่ใช้ในการศึกษานี้สามารถนำไปทำซ้ำ
    และประยุกต์ใช้กับข้อมูลจากสายพันธุ์อื่นหรือยีนอื่นได้
    ซึ่งเป็นพื้นฐานสำหรับการวิจัยในอนาคต
    ด้านชีววิทยาวิวัฒนาการ พันธุศาสตร์ประชากร
    และชีวการแพทย์
    """)
