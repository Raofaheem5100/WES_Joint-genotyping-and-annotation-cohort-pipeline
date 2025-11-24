# 🧬 WES Joint Genotyping & Annotation Cohort Pipeline  
### **A Complete Multi-Sample Whole Exome Sequencing (WES) Pipeline for Variant Discovery, Annotation & Prioritization**

---

## 📌 Overview  
This repository provides a fully automated **multi-sample Whole Exome Sequencing (WES) analysis pipeline** integrating:

- **GATK Best Practices**
- **Joint Genotyping of Multiple Samples**
- **ANNOVAR Functional Annotation**
- **Multi-step Biological Filtering**
- **CMD/CM Gene-Panels Prioritization**
- **Comprehensive Excel/TSV Outputs**

Designed for **rare pathogenic variant discovery** in large cohorts with high accuracy and efficiency.

---

## 🚀 Key Features  

### 🔹 1. FASTQ → GVCF (Per Sample)
- Fastp QC + adapter trimming  
- BWA-MEM alignment to GRCh38  
- Samtools sorting & read-group assignment  
- Picard MarkDuplicates  
- GATK BaseRecalibrator & ApplyBQSR  
- GATK HaplotypeCaller (GVCF mode)

### 🔹 2. Joint Genotyping
- GATK CombineGVCFs  
- GATK GenotypeGVCFs  
- Cohort-level raw VCF generation  

### 🔹 3. Artifact & Hard Filters
- QD, FS, MQ, SOR, RankSum filters  
- PASS filtering  
- High-complexity region (HCR) removal  
- DP/AD/VAF genotype-level filtering  
- Remove sites with no ALT allele  

### 🔹 4. ANNOVAR Annotation  
Integrated annotation from:  
- refGene  
- ClinVar  
- gnomAD  
- avsnp150  
- dbNSFP 4.2a  

### 🔹 5. Biological Filtering Framework  
- Remove intronic variants  
- Keep only exonic, splicing, UTRs  
- Non-synonymous variant selection  
- Population AF < 1%  
- ≥3 damaging predictors (SIFT, PolyPhen, PROVEAN, MT, LRT)  
- CMD/CM gene-panel prioritization  
- Excel/TSV final reports  

---

## 📁 Folder Structure  

```
WES_Joint-genotyping-and-annotation-cohort-pipeline/
│
├── multi_wes_cmd_cm4.sh          # Full pipeline script
├── CMD_genes2.txt                # CMD gene-panel list
├── CM_genes2.txt                 # CM gene-panel list
├── runs.txt                      # FASTQ run accession list
├── README.md                     # Documentation
└── OUTPUT/                       # Results after pipeline run
```

---

## 🛠 Requirements  

### **Software used**
| Tool | Version |
|------|---------|
| GATK | 4.x |
| BWA-MEM | 0.7.x |
| Samtools | 1.15+ |
| Picard | 3.x |
| Fastp | 0.23+ |
| ANNOVAR | Latest |
| bcftools | 1.15+ |
| bedtools | 2.30+ |
| R (data.table, openxlsx) | latest |
| Python3 | with pandas |

### **Reference Genome**
- **GRCh38 / hg38 FASTA**  
- dbSNP, Mills + 1000G indels  
- High-complexity regions BED (optional)

---

## 📥 Input Files  

### 1. FASTQ files (paired-end)  
Automatically downloaded using **runs.txt**

Format:

```
run_accession    fastq_ftp_urls
SRR3671535       ftp://...;ftp://...
```

---

## 🔄 Pipeline Workflow  

```
FASTQ
  ↓
QC + Trimming (fastp)
  ↓
BWA Alignment → BAM Cleanup
  ↓
GATK BaseRecalibrator (BQSR)
  ↓
HaplotypeCaller (GVCF)
  ↓
Joint Genotyping (CombineGVCFs + GenotypeGVCFs)
  ↓
Artifact & Hard Filters
  ↓
ANNOVAR Annotation
  ↓
Biological Prioritization (CMD/CM)
  ↓
Final Excel + TSV Reports
```

---

## 👨‍💻 Author  
Rao Faheem 
Bioinformatics | NGS | Variant Calling | scRNA-seq | Machine Learning  

---

## 📄 License  
This project is open-source under the **MIT License**.

---

## ⭐ If you find this useful, please star the repository!  
