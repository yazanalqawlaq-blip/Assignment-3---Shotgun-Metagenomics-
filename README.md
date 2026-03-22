# Assignment 3: Shotgun Metagenomics Analysis of Vegan and Omnivore Gut Microbiomes

---
 
## Introduction
 
The human gut microbiome is a complex community of microorganisms that plays a central role in host metabolism, immune function, and overall health (Turnbaugh et al., 2006). Among the many factors that shape gut microbial composition, diet is one of the most well-documented drivers. Long-term dietary patterns have been shown to influence both the taxonomic structure and functional capacity of the gut microbiome (Wu et al., 2011; David et al., 2014). In particular, plant-based diets tend to enrich for fiber-degrading taxa and short-chain fatty acid (SCFA) producers, while omnivorous diets that include red meat have been associated with protein-fermenting and bile-tolerant microorganisms (De Filippo et al., 2010; David et al., 2014). A recent large-scale study across 21,561 individuals found that microbial profiles could distinguish vegan, vegetarian, and omnivore diets with high accuracy, and that omnivore-associated taxa such as *Ruminococcus torques* and *Bilophila wadsworthia* correlated with less favourable cardiometabolic health markers (Fackelmann et al., 2025).
 
The present analysis uses shotgun metagenomic data from De Filippis et al. (2019), who sequenced the gut metagenomes of 97 healthy Italian adults following vegan (n=36), vegetarian (n=38), or omnivore (n=23) diets. Their study focused on strain-level diversity of *Prevotella copri* across dietary groups, and the raw sequencing data are publicly available on NCBI under accession SRP126540. For this assignment, I selected 4 vegan and 4 omnivore samples with comparable sequencing depth (~9-11 GB each) to compare gut microbial community composition, diversity, and differential abundance between the two dietary groups.
 
Shotgun metagenomics was the sequencing approach used in this study, and it offers several advantages over 16S rRNA amplicon sequencing (metabarcoding) for this type of analysis. While 16S sequencing is cheaper and benefits from more complete reference databases, it is limited to taxonomic profiling at genus-level resolution, is subject to PCR amplification biases, and cannot capture functional information directly from the reads (Knight et al., 2018). Shotgun metagenomics sequences all DNA in a sample, allowing for species- and strain-level resolution, detection of functional genes, and avoidance of primer bias. The trade-off is that shotgun approaches require deeper sequencing and larger reference databases for adequate classification.
 
For taxonomic classification, I used Kraken2 (Wood et al., 2019), a kmer-based classifier that assigns taxonomic labels by matching exact subsequences against a reference database. Kraken2 was chosen over alignment-based alternatives such as MetaPhlAn (Beghini et al., 2021), which relies on clade-specific marker genes and therefore has a smaller database footprint but can miss taxa not represented by markers. Centrifuge (Kim et al., 2016) was also considered, but Kraken2 is better established in the course material and is well-suited to paired-end short read data. One limitation of Kraken2 is its tendency to produce false positives at lower confidence thresholds, which I addressed by setting `--confidence 0.15` as recommended in the course lectures. Abundance re-estimation was performed with Bracken (Lu et al., 2017), which redistributes reads assigned at higher taxonomic levels back down to species level while accounting for genome size differences. For downstream analysis in R, I used phyloseq (McMurdie and Holmes, 2013) for data handling and visualization, the vegan package (Oksanen et al., 2024) for diversity metrics and PERMANOVA, ALDEx2 (Fernandes et al., 2014) for differential abundance analysis, and ANCOM-BC2 (Lin and Peddada, 2023) as a second differential abundance method. Both ALDEx2 and ANCOM-BC2 are designed to handle the compositionality problem inherent to metagenomic count data, where counts represent relative rather than absolute abundances and are therefore not independent (Gloor et al., 2017). Methods like DESeq2 or edgeR, which were originally designed for RNA-seq, have been shown to produce inflated false positive rates when applied to compositional microbiome data (Nearing et al., 2022), which is why I opted for tools specifically developed for this purpose.
 
---
 
## Methods
 
### Data Acquisition
 
Raw paired-end shotgun metagenomic reads were downloaded from the NCBI Sequence Read Archive (SRA) under study accession SRP126540 (De Filippis et al., 2019). I selected 4 vegan and 4 omnivore samples (Table 1) based on comparable sequencing depth to minimize bias in downstream comparisons. All samples were sequenced on the Illumina NextSeq 500 platform, generating 2x150 bp paired-end reads. Data was downloaded using the SRA Toolkit's `prefetch` and `fasterq-dump` commands on the Narval HPC cluster (Compute Canada).
 
**Table 1. Samples used in this analysis.**
 
| SRR Accession | Diet | Subject ID | City | Raw Reads (M) |
|:---|:---|:---|:---|:---|
| SRR8146974 | Vegan | VOV113 | Turin | 70.8 |
| SRR8146973 | Vegan | VOV114 | Turin | 69.2 |
| SRR8146965 | Vegan | 11BA | Bari | 65.7 |
| SRR8146961 | Vegan | 17PR | Parma | 63.8 |
| SRR8146969 | Omnivore | VOV39 | Turin | 73.5 |
| SRR8146975 | Omnivore | VOV77 | Turin | 71.2 |
| SRR8146970 | Omnivore | VOV36 | Turin | 62.6 |
| SRR8146976 | Omnivore | VOV70 | Turin | 57.2 |
 
### Quality Control
 
Raw reads were trimmed and quality-filtered using fastp v1.0.1 (Chen et al., 2018). Parameters were set as follows: minimum base quality of Phred 20 (`--qualified_quality_phred 20`), minimum read length of 50 bp (`--length_required 50`), sliding window trimming from both ends with a window size of 4 and a mean quality threshold of 20 (`--cut_front`, `--cut_tail`, `--cut_window_size 4`, `--cut_mean_quality 20`), automatic adapter detection for paired-end reads (`--detect_adapter_for_pe`), and base correction in overlapping regions (`--correction`). A Phred score of 20 corresponds to 99% base call accuracy and is a standard cutoff for Illumina data. The minimum read length of 50 bp was chosen to ensure reads are long enough for reliable kmer-based classification with Kraken2 (default kmer length = 35).
 
### Taxonomic Classification
 
Quality-filtered reads were classified using Kraken2 v2.1.6 (Wood et al., 2019) against the Standard-8 pre-built database (downloaded October 2025). This database includes RefSeq genomes for bacteria, archaea, viruses, plasmids, human, and UniVec_Core sequences, and is designed to run in 8 GB of RAM. A confidence threshold of 0.15 was used (`--confidence 0.15`), which requires that at least 15% of the kmers in a read must map to a given taxon for the classification to be accepted. The default confidence of 0 accepts any classification regardless of how few kmers match, which tends to inflate false positive rates, particularly with reduced databases. The course instructor recommended a confidence of 0.15 or higher to mitigate this. The `--paired` and `--gzip-compressed` flags were used to specify the input format.
 
### Abundance Re-estimation
 
Bracken v3.0 (Lu et al., 2017) was used to redistribute reads classified at higher taxonomic levels (genus or above) down to the species level. This step accounts for the fact that closely related species share large portions of their genomes, so many reads cannot be uniquely assigned to a single species by Kraken2. Bracken uses the kmer distribution information in the Kraken2 database to probabilistically re-assign these ambiguous reads. I ran Bracken at species level (`-l S`) with a read length of 150 bp (`-r 150`) matching our data, and a minimum read threshold of 10 (`-t 10`) to filter out very low-abundance taxa that are likely noise.
 
### BIOM File Construction
 
Individual Bracken species-level reports for all 8 samples were combined into a single BIOM-format abundance table using kraken-biom (Dabdoub, 2016). The resulting BIOM file was imported into R using the biomformat package.
 
### Diversity and Statistical Analysis
 
All downstream analysis was performed in R v4.5.0 using RStudio. The BIOM file was imported into a phyloseq object (McMurdie and Holmes, 2013), and sample metadata (diet group, city of origin, subject ID) was added manually. Taxonomy column names were standardized from Rank1-Rank7 to Kingdom, Phylum, Class, Order, Family, Genus, and Species.
 
**Alpha diversity** was measured using three complementary metrics: Observed species richness (counts the number of species present), Shannon index (accounts for both richness and evenness), and Simpson index (measures dominance, giving more weight to abundant species). These three metrics cover different aspects of alpha diversity as categorized by Cassol et al. (2024): richness, information-based, and dominance-based measures, respectively. Differences between vegan and omnivore groups were tested using Wilcoxon rank-sum tests.
 
**Beta diversity** was assessed using Bray-Curtis dissimilarity, which considers both species identity and abundance when comparing samples. This was chosen over Jaccard distance, which only uses presence/absence information. Both PCoA and NMDS ordinations were generated to visualize community-level differences. PERMANOVA (adonis2 from the vegan package) was used to test whether diet group explained a significant portion of the variation in community composition, with 999 permutations (Anderson, 2001).
 
**Differential abundance** was assessed using two methods. First, ALDEx2 v1.38.0 (Fernandes et al., 2014), which uses a Bayesian approach with Monte Carlo sampling (128 instances) to estimate technical variation, followed by Welch's t-test with Benjamini-Hochberg correction. ALDEx2 applies a centered log-ratio (CLR) transformation to address compositionality. Second, ANCOM-BC2 v2.8.0 (Lin and Peddada, 2023), which estimates and corrects for sample-specific biases in microbial absolute abundance. ANCOM-BC2 was run at genus level with Holm correction for multiple testing, a prevalence filter of 10% (`prv_cut = 0.10`), a library size cutoff of 1000 reads (`lib_cut = 1000`), and structural zero detection enabled (`struc_zero = TRUE`). Both methods were applied because they are considered conservative tools that minimize false positives in compositional data (Nearing et al., 2022).
 
### Software Versions
 
**Table 2. Software and versions used.**
 
| Software | Version | Purpose | Reference |
|:---|:---|:---|:---|
| SRA Toolkit | 3.0.9 | Data download | NCBI |
| fastp | 1.0.1 | Quality control and trimming | Chen et al., 2018 |
| Kraken2 | 2.1.6 | Taxonomic classification | Wood et al., 2019 |
| Bracken | 3.0 | Abundance re-estimation | Lu et al., 2017 |
| kraken-biom | 1.2.0 | BIOM file generation | Dabdoub, 2016 |
| R | 4.5.0 | Statistical analysis | R Core Team, 2024 |
| phyloseq | 1.50.0 | Microbiome data handling | McMurdie and Holmes, 2013 |
| vegan | 2.7-3 | Diversity and PERMANOVA | Oksanen et al., 2024 |
| ALDEx2 | 1.38.0 | Differential abundance | Fernandes et al., 2014 |
| ANCOM-BC2 | 2.8.0 | Differential abundance | Lin and Peddada, 2023 |
 
---
