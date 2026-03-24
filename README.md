# Assignment 3: Shotgun Metagenomics Analysis of Vegan and Omnivore Gut Microbiomes

**Course:** BINF*6110 - Bioinformatics

**Author:** Yazan Alqawlaq

**Date:** March 23, 2026

**Repository:** https://github.com/yazanalqawlaq-blip/Assignment-3---Shotgun-Metagenomics-

---

## Table of Contents
- [Introduction](#introduction)
- [Project Structure](#project-structure)
- [Methods](#methods)
- [Results](#results)
- [Discussion](#discussion)
- [Limitations](#limitations)
- [Conclusion](#conclusion)
- [References](#references)


---

## Introduction

The human gut harbors a dense and diverse community of microorganisms, collectively referred to as the gut microbiome, that plays a well-documented role in immune function and metabolic health (Clemente et al., 2012; Turnbaugh et al., 2006). Among the environmental factors that shape microbial community composition, habitual diet is one of the most consistent and reproducible drivers. Long-term dietary patterns influence which microbial taxa colonize the gut and what metabolic functions they perform (Wu et al., 2011; David et al., 2014). Plant-based diets, which are rich in dietary fiber and complex carbohydrates, tend to select for saccharolytic and fiber-degrading taxa that produce short-chain fatty acids like butyrate and propionate, which are known to support gut barrier integrity and reduce inflammation (De Filippo et al., 2010). In contrast, omnivorous diets that include red meat and animal fat have been associated with enrichment of protein-fermenting and bile-tolerant taxa such as *Bilophila wadsworthia* and *Alistipes putredinis*, some of which have been linked to poorer cardiometabolic outcomes (David et al., 2014; Fackelmann et al., 2025). A recent cross-sectional study of 21,561 individuals across five cohorts from three countries found that microbial profiles could distinguish vegan, vegetarian, and omnivore diets with a mean AUC of 0.85, and that red meat consumption was the strongest dietary driver of omnivore-associated microbiome signatures (Fackelmann et al., 2025).

The present analysis uses shotgun metagenomic data from De Filippis et al. (2019), who sequenced the gut metagenomes of 97 healthy Italian adults following vegan (n=36), vegetarian (n=38), or omnivore (n=23) diets. Their study focused specifically on strain-level diversity of *Prevotella copri* across dietary groups and showed that fiber-rich diets were associated with *P. copri* strains carrying enhanced carbohydrate catabolism genes, while omnivore-associated strains had a higher prevalence of genes involved in branched-chain amino acid biosynthesis, a risk factor for glucose intolerance (De Filippis et al., 2019). The raw sequencing data are publicly available on NCBI under accession SRP126540. For this assignment, I selected 4 vegan and 4 omnivore samples with comparable sequencing depth (~9-11 GB each) to compare gut microbial community composition, diversity, and differential abundance between the two dietary groups.

I chose to work with shotgun metagenomic data rather than 16S rRNA amplicon data because shotgun sequencing captures all DNA in a sample, allowing for species- and strain-level taxonomic resolution without the primer biases and limited phylogenetic resolution inherent to amplicon-based approaches (Knight et al., 2018). While 16S rRNA metabarcoding is cheaper, benefits from more complete reference databases, and is sufficient for many community profiling questions, it is restricted to taxonomic classification only and typically resolves taxa to genus level at best. Shotgun metagenomics also has the potential for functional annotation of gene content, though I did not perform functional profiling in this analysis. That said, shotgun approaches require much deeper sequencing and are more dependent on database completeness for accurate classification.

For taxonomic classification, I used Kraken2 (Wood et al., 2019), a kmer-based classifier that assigns taxonomic labels by matching minimizers derived from 35-bp kmers against a pre-built hash table reference database. I chose Kraken2 over alignment-based alternatives such as MetaPhlAn (Beghini et al., 2021), which uses clade-specific marker genes to produce species-level profiles with a much smaller database footprint but can miss taxa not represented in its marker gene set. Centrifuge (Kim et al., 2016) was also considered as a kmer-based alternative, but Kraken2 is more widely adopted, better supported in the course material, and well-suited to paired-end Illumina short read data. One recognized limitation of Kraken2 is its tendency toward false positive classifications at lower confidence thresholds, especially with reduced databases where reads from unrepresented species may be spuriously assigned to related taxa in the database (Wood et al., 2019). I addressed this by setting `--confidence 0.15`, which requires at least 15% of the kmers in a read to support a classification before it is accepted, as recommended by the lecture notes. Species-level abundance re-estimation was performed with Bracken (Lu et al., 2017), which redistributes reads that Kraken2 assigned at genus or higher levels back down to species level by using kmer distribution data to account for the proportion of each genome that is unique versus shared with closely related species.

For downstream statistical analysis in R, I used phyloseq (McMurdie and Holmes, 2013) and the vegan package (Oksanen et al., 2024) for data handling, visualization, and diversity analyses, alongside two complementary differential abundance tools: ALDEx2 (Fernandes et al., 2014) and ANCOM-BC2 (Lin and Peddada, 2024). The choice of differential abundance method is particularly relevant for microbiome data because count tables generated from sequencing are inherently compositional: they represent relative rather than absolute abundances, meaning that an increase in one taxon's proportion necessarily decreases others (Gloor et al., 2017). Standard RNA-seq differential expression tools such as DESeq2 or edgeR do not account for this property and have been shown to produce inflated false positive rates when applied to microbiome abundance data (Nearing et al., 2022). I therefore used ALDEx2, which applies a centered log-ratio transformation to address compositionality and uses Monte Carlo sampling to estimate technical variation, and ANCOM-BC2, which models and corrects for sample-specific biases in absolute abundance estimation. Both methods are considered conservative but reliable, and running both let me check whether findings hold up under both approaches.

---

## Project Structure

```
Assignment-3---Shotgun-Metagenomics-/
├── README.md
├── .gitignore
├── scripts/
│   ├── 01_download_data.sh
│   ├── 02_qc_fastp.sh
│   ├── 03_download_kraken2_db.sh
│   ├── 04_kraken2_classify.sh
│   ├── 05_bracken_reestimate.sh
│   ├── 06_prep_biom.sh
│   └── 07_R_analysis.R
├── output_files/
│   ├── qc/                      
│   ├── kraken2/                
│   ├── bracken/                 
│   └── R_analysis/
│       └── combined_bracken.biom
└── input_data/ # gitignored (large files not included)                  
```

---


## Methods

### Data Acquisition

Raw paired-end shotgun metagenomic reads were downloaded from the NCBI Sequence Read Archive (SRA) under study accession SRP126540 (De Filippis et al., 2019). I selected 4 vegan and 4 omnivore samples (Table 1) with comparable sequencing depth to minimize systematic bias in downstream diversity and abundance comparisons. All samples were sequenced on the Illumina NextSeq 500 platform, generating 2x150 bp paired-end reads. I downloaded the data using SRA Toolkit's `prefetch` and `fasterq-dump` commands on the Narval HPC cluster login node (Compute Canada), as compute nodes do not have internet access. Files were compressed with `gzip` immediately after extraction to conserve disk space.

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

Raw reads were trimmed and quality-filtered using fastp v1.0.1 (Chen et al., 2018) with the following parameters: minimum base quality of Phred 20 (`--qualified_quality_phred 20`), minimum read length of 50 bp (`--length_required 50`), sliding window trimming from both ends (`--cut_front`, `--cut_tail`) with a window size of 4 and a mean quality threshold of 20 (`--cut_window_size 4`, `--cut_mean_quality 20`), automatic adapter detection for paired-end reads (`--detect_adapter_for_pe`), and base correction in overlapping regions of paired reads (`--correction`). I set the minimum read length to 50 bp because Kraken2 uses a default kmer length of 35 bp, meaning that reads shorter than approximately 50 bp would contain too few informative kmers for reliable classification. The `--correction` flag was included because it takes advantage of the overlapping region in paired-end reads to correct mismatched bases, which can reduce error-driven misclassifications during kmer matching.

### Taxonomic Classification

Quality-filtered reads were classified using Kraken2 v2.1.6 (Wood et al., 2019) against the Standard-8 pre-built database (released October 2025). This database includes RefSeq genomes for bacteria, archaea, viruses, plasmids, human, and UniVec_Core sequences, and is designed to run within 8 GB of RAM. I chose this database because it is the same one used in the course tutorials and provides a reasonable balance between taxonomic coverage and computational requirements. The full Kraken2 Standard database would provide more comprehensive classification but requires over 70 GB of RAM, while more specialized databases such as GTDB (Parks et al., 2022) offer more thorough bacterial and archaeal taxonomy but exclude viral sequences entirely. A confidence threshold of 0.15 was applied (`--confidence 0.15`) as described in the Introduction and recommended by the course notes, to reduce false positive classifications that occur when unrepresented taxa are spuriously assigned to their closest database relative. The `--paired` and `--gzip-compressed` flags were used to specify the input data format.

### Abundance Re-estimation

Bracken v3.0 (Lu et al., 2017) was used to redistribute reads classified at higher taxonomic levels (genus or above) down to species level. This step is necessary because many shotgun reads map to genomic regions shared between closely related species, meaning Kraken2 assigns them conservatively at the genus or family level rather than committing to a species. Bracken uses the kmer distribution information in the Kraken2 database to probabilistically reassign these ambiguous reads based on the proportion of each genome that is uniquely identifiable at the kmer level. I ran Bracken at species level (`-l S`) with a read length parameter of 150 bp (`-r 150`), matching the actual read length of the NextSeq 500 data. Bracken uses this parameter to look up the appropriate precomputed kmer distribution file (`database150mers.kmer_distrib`), which was included in the Standard-8 database download. The minimum read threshold was set to 10 (`-t 10`) to filter out taxa supported by very few reads, which are more likely to be spurious classifications than genuine low-abundance community members.

### BIOM File Construction

Individual Bracken species-level reports for all 8 samples were combined into a single BIOM-format abundance table using kraken-biom v1.2.0 (Dabdoub, 2016). BIOM is a standardized format for storing biological observation matrices and is the expected input format for phyloseq and other microbiome analysis packages in R. The resulting BIOM file contained 225 species across 8 samples.

### Diversity and Statistical Analysis

All downstream analysis was performed in R v4.5.0 using RStudio. The BIOM file was imported into a phyloseq object using the biomformat package (McMurdie and Holmes, 2013), and sample metadata including diet group, city of origin, and subject ID were added manually, since kraken-biom does not carry sample-level metadata. Taxonomy column names generated by kraken-biom (Rank1 through Rank7) were renamed to standard levels (Kingdom, Phylum, Class, Order, Family, Genus, Species).

**Alpha diversity** was measured using three metrics chosen to capture different aspects of within-sample diversity following the framework described by Cassol et al. (2025): Observed species richness, Shannon index, and Simpson index. These three were chosen to cover three of the four major alpha diversity categories described by Cassol et al. (2025), specifically richness, information, and dominance. I did not include phylogenetic diversity (Faith's PD) because the Kraken2/Bracken pipeline does not produce a phylogenetic tree, which is required for this metric. Differences between vegan and omnivore groups were tested using Wilcoxon rank-sum tests, which are appropriate for small sample sizes and do not assume normally distributed data.

**Beta diversity** was assessed using Bray-Curtis dissimilarity, which quantifies the difference between two samples based on both species identity and abundance. I chose Bray-Curtis over Jaccard because Jaccard only considers presence/absence and would discard the abundance information in the Bracken count data. Both PCoA and NMDS ordinations were generated to visualize community-level differences between diet groups. I included both ordination approaches because the course lectures present both, and comparing them helps confirm that observed patterns are not artifacts of the ordination method. For NMDS, stress values below 0.1 were considered indicative of a good fit. PERMANOVA was performed using `adonis2()` from the vegan package (Anderson, 2001) with 999 permutations to test whether diet group explained a significant portion of the variation in Bray-Curtis community dissimilarity.

**Differential abundance** was assessed using two complementary approaches. First, ALDEx2 v1.38.0 (Fernandes et al., 2014) was applied at genus level with 128 Monte Carlo instances (`mc.samples = 128`) to estimate the technical variation introduced by random sampling. ALDEx2 applies a centered log-ratio (CLR) transformation to address compositionality before testing. Between-group differences were evaluated using Welch's t-test with Benjamini-Hochberg correction for multiple comparisons. I used 128 Monte Carlo instances because this is the recommended default for exploratory analysis; higher values (e.g., 1000) provide more precise estimates but at substantially increased computation time, and the results are generally stable at 128 instances (Fernandes et al., 2014). Second, ANCOM-BC2 v2.8.0 (Lin and Peddada, 2024) was run at genus level with Holm correction for multiple testing, a prevalence filter of 10% (`prv_cut = 0.10`, meaning taxa must be present in at least 10% of samples to be tested), a library size cutoff of 1000 reads (`lib_cut = 1000`, removing samples with very low total counts), and structural zero detection enabled (`struc_zero = TRUE`). The prevalence filter was set at 10% rather than a higher threshold because with only 8 samples, even a taxon present in a single sample represents 12.5% prevalence, and I did not want to exclude taxa that might be genuinely absent from one diet group. Structural zeros are taxa that are completely absent from one experimental group but present in the other, and ANCOM-BC2 identifies and reports these separately from the main differential abundance test because they violate the assumptions of the log-linear model used for estimation (Lin and Peddada, 2024).

### Software Versions

**Table 2. Software and versions used.**
| Software | Version | Purpose |
|:---|:---|:---|
| SRA Toolkit | 3.0.9 | Data download and FASTQ conversion |
| fastp | 1.0.1 | Quality control and adapter trimming |
| Kraken2 | 2.1.6 | Taxonomic classification |
| Bracken | 3.0 | Species-level abundance re-estimation |
| kraken-biom | 1.2.0 | BIOM file generation |
| R | 4.5.0 | Statistical computing environment |
| phyloseq | 1.50.0 | Microbiome data handling and visualization |
| vegan | 2.7-3 | Diversity metrics and PERMANOVA |
| ALDEx2 | 1.38.0 | Compositional differential abundance |
| ANCOM-BC2 | 2.8.0 | Bias-corrected differential abundance |

---

## Results

### Quality Control

Quality filtering with fastp retained a high proportion of reads across all samples (Table 3). Seven of eight samples had pass rates above 97.7%, with Q20 scores exceeding 96% and Q30 scores above 91%. SRR8146976 had somewhat lower quality metrics (Q30 = 83.3%, pass rate = 96.0%), likely due to a less optimal sequencing run, but still retained 54.9 million reads after filtering, which I considered sufficient for downstream classification.

**Table 3. Quality control summary.**

| Sample | Diet | Raw Reads (M) | Passed (M) | Pass % | Q20 % | Q30 % |
|:---|:---|:---|:---|:---|:---|:---|
| SRR8146974 | Vegan | 70.8 | 69.7 | 98.4 | 96.7 | 91.4 |
| SRR8146973 | Vegan | 69.2 | 68.1 | 98.4 | 96.7 | 91.5 |
| SRR8146965 | Vegan | 65.7 | 64.7 | 98.6 | 96.9 | 92.1 |
| SRR8146961 | Vegan | 63.8 | 62.4 | 97.8 | 96.7 | 91.7 |
| SRR8146969 | Omnivore | 73.5 | 72.5 | 98.7 | 97.2 | 92.6 |
| SRR8146975 | Omnivore | 71.2 | 69.8 | 98.1 | 96.5 | 91.1 |
| SRR8146970 | Omnivore | 62.6 | 61.4 | 98.1 | 96.8 | 91.8 |
| SRR8146976 | Omnivore | 57.2 | 54.9 | 96.0 | 93.3 | 83.3 |

### Taxonomic Classification

Classification rates with Kraken2 were low across all samples, with the large majority of reads remaining unclassified (over 96% in all samples). For example, SRR8146974 had approximately 34.9 million quality-filtered reads, of which only around 265,000 were classified at the species level and subsequently processed by Bracken. This low classification rate is expected with the Standard-8 database, which is a reduced version of the full Kraken2 standard database and contains only a fraction of the reference genomes available in larger databases. After Bracken re-estimation, 225 species were identified across the 8 samples, with a mean of 107 species per sample (range: 91 to 123). Rarefaction curves for all eight samples reached an asymptote well before the total number of classified reads, confirming that the sequencing depth was adequate to capture the species diversity present in each sample despite the low overall classification rate (Figure 1).

<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/1a49dff8-85b9-4311-be85-ded4b2f4bf4f" />

**Figure 1. Rarefaction curves for all eight samples.** Each curve represents one sample, showing the number of species detected as a function of sequencing depth. All curves plateau well before reaching the total number of classified reads, indicating that sequencing depth was sufficient to capture the majority of species diversity in each sample.

### Taxonomic Composition

At the phylum level, both vegan and omnivore gut microbiomes were dominated by Bacillota (formerly Firmicutes) and Bacteroidota (formerly Bacteroidetes), which together accounted for the majority of classified reads in all samples (Figure 2). This is expected, as these two phyla consistently dominate the healthy human gut (Eckburg et al., 2005). Actinomycetota was visibly more abundant in the vegan samples, driven primarily by *Bifidobacterium* species. Other phyla detected at lower relative abundances included Verrucomicrobiota (represented almost entirely by *Akkermansia muciniphila*), Pseudomonadota, and Thermodesulfobacteriota.


<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/08797263-3319-4c46-b68c-7a640f9e3002" />

**Figure 2. Phylum-level taxonomic composition of vegan and omnivore gut microbiomes.** Stacked bar plot showing relative abundance of major phyla (those exceeding 1% relative abundance in any sample) in each sample, faceted by diet group. Both groups are dominated by Bacillota and Bacteroidota. Actinomycetota is visibly enriched in vegan samples.



At the family level (Figure 3), differences between groups became more apparent. Bacteroidaceae, Lachnospiraceae, and Oscillospiraceae were consistently abundant across both groups. Prevotellaceae showed high inter-individual variability, being particularly abundant in certain samples from both groups (e.g., SRR8146961, a vegan subject from Parma, and SRR8146969, an omnivore from Turin). Bifidobacteriaceae was noticeably more abundant in the vegan samples, consistent with the phylum-level pattern.


<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/c0c5976b-8059-4cf1-b482-c67ae5106fcb" />

**Figure 3. Family-level taxonomic composition.** Stacked bar plot showing relative abundance of the top 15 families by mean abundance, with remaining families grouped as "Other." Faceted by diet group. Bifidobacteriaceae is enriched in vegans. Prevotellaceae abundance is highly variable across individuals in both groups.


### Alpha Diversity

Alpha diversity did not differ significantly between vegan and omnivore groups for any of the three metrics tested (Figure 4, Table 4). Observed species richness was nearly identical between groups (Vegan mean = 108.8, Omnivore mean = 105.0; Wilcoxon p = 1.00). Shannon diversity was slightly higher in vegans (3.01 vs. 2.90; p = 0.49), and Simpson index values were also comparable (0.891 vs. 0.876; p = 0.69). None of these comparisons reached statistical significance, indicating that vegan and omnivore gut communities in this subset of the De Filippis et al. cohort harbor similar levels of species richness and evenness.

<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/3e9a1d8b-e1b6-4ccc-b37a-28ae81134198" />

**Figure 4. Alpha diversity comparison between vegan and omnivore groups.** Boxplots of Observed species richness, Shannon index, and Simpson index, with individual sample points overlaid. No significant differences were detected for any metric (Wilcoxon rank-sum test, all p > 0.05).



**Table 4. Alpha diversity summary statistics.**

| Metric | Vegan Mean | Omnivore Mean | Wilcoxon p-value |
|:---|:---|:---|:---|
| Observed | 108.8 | 105.0 | 1.000 |
| Shannon | 3.01 | 2.90 | 0.486 |
| Simpson | 0.891 | 0.876 | 0.686 |

### Beta Diversity

PCoA ordination based on Bray-Curtis dissimilarity showed partial separation between vegan and omnivore samples along the first two principal coordinate axes, though with overlap between groups (Figure 5). NMDS ordination produced a similar pattern (Figure 6), with a stress value of 0.048, indicating an excellent fit. The consistency between PCoA and NMDS confirms that the partial separation between diet groups is not an artifact of either ordination method.

PERMANOVA revealed that diet explained 21.2% of the total variation in community composition (R² = 0.212, F = 1.61), though this result did not reach statistical significance (p = 0.108, 999 permutations; Table 5). The R² value is a meaningful effect size for a microbiome study, where many factors beyond diet (age, genetics, medication use, individual history) contribute to community variation, but the small sample size (n = 4 per group) limits the ability to detect this effect as statistically significant.

<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/0e4172f3-a153-4668-8bd9-2199950612b5" />

**Figure 5. PCoA ordination of gut microbial communities based on Bray-Curtis dissimilarity.** Points represent individual samples, colored by diet group. Ellipses show 95% confidence intervals. Partial separation between vegan and omnivore groups is visible. PERMANOVA: R² = 0.212, F = 1.61, p = 0.108 (999 permutations).

<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/0a3c5cde-3207-4e22-9d3e-c60c87fc6e11" />

**Figure 6. NMDS ordination based on Bray-Curtis dissimilarity.** Stress = 0.048, indicating excellent ordination fit. Diet groups show partial separation consistent with the PCoA result.



**Table 5. PERMANOVA results (adonis2, 999 permutations).**

| Source | Df | Sum of Squares | R² | F | p-value |
|:---|:---|:---|:---|:---|:---|
| Diet | 1 | 0.354 | 0.212 | 1.611 | 0.108 |
| Residual | 6 | 1.319 | 0.788 | | |
| Total | 7 | 1.673 | 1.000 | | |

### Differential Abundance

Neither ALDEx2 nor ANCOM-BC2 identified any taxa as significantly differentially abundant between vegan and omnivore groups after multiple testing correction (ALDEx2: all BH-adjusted p > 0.05; ANCOM-BC2: all Holm-adjusted q > 0.05). Both tools issued warnings about the small sample size: ANCOM-BC2 specifically noted that "variance estimation would be unstable when the sample size is < 5 per group."

Despite the lack of statistically significant results from the main differential abundance tests, ANCOM-BC2 detected 15 structural zeroes (Table 6). These are taxa that were completely absent from all samples in one diet group while being present in one or more samples in the other group. Several of these taxa have well-documented associations with dietary patterns.

**Table 6. Selected structural zeroes detected by ANCOM-BC2.** Taxa present in one diet group but absent from all samples in the other.

| Taxon | Present in Omnivore | Present in Vegan |
|:---|:---|:---|
| *Escherichia* | Yes | No |
| *Desulfovibrio* | Yes | No |
| *Streptococcus* | Yes | No |
| *Cutibacterium* | Yes | No |
| *Oxalobacter* | Yes | No |
| *Amedibacillus* | Yes | No |
| *Christensenella* | No | Yes |
| *Eisenbergiella* | No | Yes |
| *Megamonas* | No | Yes |
| *Wansuia* | No | Yes |
| *Haemophilus* | No | Yes |

ALDEx2 effect size analysis showed trends in taxon abundance between groups even though none reached statistical significance. The genera with the largest negative effect sizes (indicating trends toward higher abundance in omnivores) included *Phocaeicola* (effect = -0.65) and *Parabacteroides* (effect = -0.63). On the vegan side, *Bifidobacterium* showed a strong trend toward enrichment (ANCOM-BC2 log fold change = 2.98, q = 0.68), and the raw count data confirmed this: *Bifidobacterium adolescentis* read counts ranged from 13,683 to 44,505 in vegan samples compared to 0 to 1,246 in omnivore samples.

<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/9f602bc5-1469-40f5-b0be-8008c28de73a" />

**Figure 7. ANCOM-BC2 differential abundance results at genus level.** Forest plot showing log fold change (Vegan vs. Omnivore) for the top 20 genera ranked by Holm-adjusted q-value. Error bars represent standard error. No taxa reached q < 0.05. Positive values indicate higher abundance in vegans.

<img width="600" height="500" alt="Image" src="https://github.com/user-attachments/assets/8c5ac854-98d9-4a24-85d4-9250e1dc3443" />

**Figure 8. ALDEx2 effect size plot.** Top 10 genera with the largest positive and negative ALDEx2 effect sizes between omnivore and vegan groups. Positive effect size indicates higher abundance in omnivores. No genera reached Benjamini-Hochberg corrected significance at p < 0.05.



---

## Discussion

### Sulfate-Reducing and Bile-Tolerant Taxa in Omnivores

The most biologically notable pattern among the structural zeroes was the exclusive presence of *Desulfovibrio* in omnivore samples. *Desulfovibrio* species are sulfate-reducing bacteria (SRB) that produce hydrogen sulfide (H2S) as a metabolic byproduct of dissimilatory sulfate reduction. H2S is cytotoxic to colonic epithelial cells and inhibits butyrate oxidation in colonocytes, and has been implicated in the pathogenesis of ulcerative colitis (Loubinoux et al., 2002; Rowan et al., 2010). The enrichment of SRB in omnivores is consistent with higher dietary intake of sulfur-containing amino acids (methionine and cysteine) from animal protein, which serve as substrates for sulfate reduction in the colon (David et al., 2014). This is consistent with Fackelmann et al. (2025), who identified *Bilophila wadsworthia*, a closely related sulfite-reducing bacterium, as one of the key signature taxa of omnivore microbiomes and showed that its abundance correlated negatively with cardiometabolic health markers. In the present dataset, *B. wadsworthia* was detected in all samples but showed numerically higher read counts in several omnivore samples, consistent with its known association with dietary saturated fat and taurine-conjugated bile acid availability (Devkota et al., 2012).

*Escherichia* was also detected exclusively in omnivore samples. While the genus includes commensal strains that are normal inhabitants of the human gut, certain *Escherichia* species, particularly adherent-invasive strains of *E. coli*, have been associated with gut inflammation and are more frequently isolated from individuals consuming Western-style diets high in animal protein (Agus et al., 2016). The complete absence of *Escherichia* from all four vegan samples, while present in all four omnivore samples, is consistent with the broader pattern of reduced Proteobacteria in plant-based diets reported by multiple studies (De Filippo et al., 2010; David et al., 2014).

### Fiber-Associated and SCFA-Producing Taxa in Vegans

On the vegan side, the strongest signal was the trend toward higher *Bifidobacterium* abundance. The raw count data was particularly striking for *Bifidobacterium adolescentis*, which had counts ranging from 13,683 to 44,505 in the vegan samples compared to 0 to 1,246 in omnivores. *Bifidobacterium* species are saccharolytic organisms that ferment dietary fiber and complex plant-derived carbohydrates to produce acetate and lactate, which are then converted to butyrate by cross-feeding with other gut bacteria such as *Faecalibacterium prausnitzii* and *Roseburia* species (Riviere et al., 2016). The enrichment of *Bifidobacterium* in plant-based diets is one of the most consistently reported findings in diet-microbiome studies (De Filippo et al., 2010), and it reflects the higher availability of fermentable substrates like inulin and resistant starch in vegan diets.

*Christensenella* was present exclusively in vegan samples. This genus is notable as one of the most heritable taxa in the human gut microbiome and has been consistently associated with lean body mass across multiple populations (Goodrich et al., 2014; Waters and Ley, 2019). Experimental work has shown that *Christensenella minuta* can reduce adiposity in gnotobiotic mice and modulates community composition when introduced into an obese-associated microbiome (Goodrich et al., 2014). Its exclusive detection in vegans in the present analysis is consistent with the association between plant-based diets and leaner body composition, though why *Christensenella* is more common in plant-based diets is not well understood.

*Megamonas*, a member of the Selenomonadaceae family, was also found only in vegan samples. *Megamonas* species are carbohydrate fermenters known to produce propionate and acetate as major end products. Their presence in vegans is consistent with the higher dietary fiber intake associated with plant-based diets, though *Megamonas* is less well studied than other SCFA producers in the human gut.

### Community-Level Patterns and Comparison with Published Data

The lack of statistically significant differences in alpha diversity between groups is consistent with the findings of De Filippis et al. (2019), who did not report major differences in species richness between dietary groups in their full cohort of 97 individuals. Fackelmann et al. (2025) did observe slightly greater richness in omnivores compared to vegans across their much larger cohort of 21,561 individuals, but noted that the differences were small and that richness alone is not a reliable indicator of microbiome health.

The PERMANOVA R² of 0.212 indicates that diet explains approximately one-fifth of the variation in community composition, which is a substantial effect size for a microbiome study where many other factors contribute to individual variation. For comparison, De Filippis et al. (2019) reported that dietary variables explained a modest but significant portion of Bray-Curtis variation in their full cohort. The failure to reach p < 0.05 in my analysis is almost certainly a power issue: with only 4 samples per group and 999 permutations, the test has limited ability to distinguish a real dietary effect from random variation.

The high inter-individual variability in *Segatella copri* (formerly *Prevotella copri*) abundance across both diet groups is worth noting. Some individuals had very high read counts (e.g., SRR8146961 with 166,246 reads, SRR8146969 with 169,596 reads) while others in the same diet group had minimal counts. This variability is consistent with the findings of De Filippis et al. (2019), who showed that *P. copri* presence was not significantly different between dietary groups at the species level but that functionally distinct strains were associated with different diets. This strain-level variation cannot be captured by species-level shotgun classification, which is a recognized limitation of the approach used here.

---

## Limitations

Several limitations should be noted when interpreting these results:

1. **Sample size.** With only 4 samples per group, statistical power for detecting differential abundance or significant beta diversity differences is limited. Both ALDEx2 and ANCOM-BC2 issued explicit warnings about unstable variance estimates at this sample size. A minimum of 10-20 samples per group is generally recommended for robust microbiome comparisons (Knight et al., 2018).

2. **Database completeness.** The Standard-8 Kraken2 database classified fewer than 4% of reads in most samples. This means the vast majority of sequencing data was not used in the analysis. A larger database such as the full Kraken2 Standard database or GTDB would likely classify substantially more reads, though at higher computational cost. The low classification rate also means that some taxa present in the samples may have been missed entirely.

3. **Compositionality.** As noted in the Introduction, metagenomic count data is compositional. While both ALDEx2 and ANCOM-BC2 account for this, the structural zero results should be interpreted cautiously: a taxon being absent from all 4 samples in one group could reflect true biological absence, but it could also result from insufficient sequencing depth or database coverage in combination with natural low abundance.

4. **No functional analysis.** This analysis is limited to taxonomic profiling. Functional annotation with tools such as HUMAnN (Beghini et al., 2021) could reveal metabolic pathway differences between diet groups that are not detectable from taxonomy alone. This would be particularly relevant for capturing the strain-level functional variation in *Segatella copri* that De Filippis et al. (2019) identified as the main diet-dependent signal in this cohort.

5. **Geographic confounding.** All four omnivore samples are from Turin, while the vegan samples span Turin, Bari, and Parma. Regional dietary traditions in Italy differ substantially, and some of the observed variation may reflect geographic differences in diet composition beyond the vegan/omnivore distinction itself.

6. **Host DNA.** *Homo sapiens* reads were detected in all samples at low levels, which is expected for fecal metagenomes and does not substantially affect the analysis, but indicates that no host depletion step was applied during sample preparation.

---

## Conclusion

This analysis compared vegan and omnivore gut microbial communities from the De Filippis et al. (2019) Italian cohort using shotgun metagenomics with Kraken2/Bracken classification and compositional differential abundance testing in R. While no taxa reached statistical significance in the differential abundance tests after multiple testing correction, the structural zero analysis and abundance trends point toward biologically coherent differences between diet groups. Omnivore-specific taxa including *Desulfovibrio* and *Escherichia* are associated with sulfate reduction and protein fermentation, while vegan-associated taxa including *Christensenella* and *Megamonas* are linked to lean body composition and carbohydrate fermentation. The strong trend toward *Bifidobacterium* enrichment in vegans is one of the most consistently reported findings in the diet-microbiome literature. These patterns are in line with published data from larger cohorts, and a study with more samples and a more complete reference database would be better positioned to detect these trends as statistically significant.

---

## References

Agus, A., Denizot, J., Thevenot, J., Martinez-Medina, M., Massier, S., Sauvanet, P., Bernalier-Donadille, A., Denis, S., Hofman, P., Bonnet, R., Billard, E., and Barnich, N. (2016). Western diet induces a shift in microbiota composition enhancing susceptibility to adherent-invasive *E. coli* infection and intestinal inflammation. *Scientific Reports*, 6, 19032. https://doi.org/10.1038/srep19032

Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32-46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x

Beghini, F., McIver, L. J., Blanco-Miguez, A., Dubois, L., Asnicar, F., Maharjan, S., Mailyan, A., Manghi, P., Scholz, M., Thomas, A. M., Valles-Colomer, M., Weingart, G., Zhang, Y., Zolfo, M., Huttenhower, C., Franzosa, E. A., and Segata, N. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife*, 10, e65088. https://doi.org/10.7554/eLife.65088

Cassol, I., Ibañez, M., and Bustamante, J. P. (2025). Key features and guidelines for the application of microbial alpha diversity metrics. *Scientific Reports*, 15, 622. https://doi.org/10.1038/s41598-024-77864-y

Chen, S., Zhou, Y., Chen, Y., and Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560

Clemente, J. C., Ursell, L. K., Parfrey, L. W., and Knight, R. (2012). The impact of the gut microbiota on human health: an integrative view. *Cell*, 148(6), 1258-1270. https://doi.org/10.1016/j.cell.2012.01.035

Dabdoub, S. (2016). kraken-biom: Enabling interoperative format conversion for Kraken results. https://github.com/smdabdoub/kraken-biom. GitHub.

David, L. A., Maurice, C. F., Carmody, R. N., Gootenberg, D. B., Button, J. E., Wolfe, B. E., Ling, A. V., Devlin, A. S., Varma, Y., Fischbach, M. A., Biddinger, S. B., Dutton, R. J., and Turnbaugh, P. J. (2014). Diet rapidly and reproducibly alters the human gut microbiome. *Nature*, 505(7484), 559-563. https://doi.org/10.1038/nature12820

De Filippo, C., Cavalieri, D., Di Paola, M., Ramazzotti, M., Poullet, J. B., Massart, S., Collini, S., Pieraccini, G., and Lionetti, P. (2010). Impact of diet in shaping gut microbiota revealed by a comparative study in children from Europe and rural Africa. *Proceedings of the National Academy of Sciences*, 107(33), 14691-14696. https://doi.org/10.1073/pnas.1005963107

De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., and Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. *Cell Host & Microbe*, 25(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004

Devkota, S., Wang, Y., Musch, M. W., Leone, V., Fehlner-Peach, H., Nadimpalli, A., Antonopoulos, D. A., Bhatt, B., and Chang, E. B. (2012). Dietary-fat-induced taurocholic acid promotes pathobiont expansion and colitis in Il10-/- mice. *Nature*, 487(7405), 104-108. https://doi.org/10.1038/nature11225

Eckburg, P. B., Bik, E. M., Bernstein, C. N., Purdom, E., Dethlefsen, L., Sargent, M., Gill, S. R., Nelson, K. E., and Relman, D. A. (2005). Diversity of the human intestinal microbial flora. *Science*, 308(5728), 1635-1638. https://doi.org/10.1126/science.1110591

Fackelmann, G., Manghi, P., Carlino, N., Heidrich, V., Piccinno, G., Ricci, L., Piperni, E., Arrè, A., Bakker, E., Creedon, A. C., Francis, L., Capdevila Pujol, J., Davies, R., Wolf, J., Bermingham, K. M., Berry, S. E., Spector, T. D., Asnicar, F., and Segata, N. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. *Nature Microbiology*, 10(1), 41-52. https://doi.org/10.1038/s41564-024-01870-z

Fernandes, A. D., Reid, J. N., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., and Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome*, 2, 15. https://doi.org/10.1186/2049-2618-2-15

Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., and Egozcue, J. J. (2017). Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology*, 8, 2224. https://doi.org/10.3389/fmicb.2017.02224

Goodrich, J. K., Waters, J. L., Poole, A. C., Sutter, J. L., Koren, O., Blekhman, R., Beaumont, M., Van Treuren, W., Knight, R., Bell, J. T., Spector, T. D., Clark, A. G., and Ley, R. E. (2014). Human genetics shape the gut microbiome. *Cell*, 159(4), 789-799. https://doi.org/10.1016/j.cell.2014.09.053

Kim, D., Song, L., Breitwieser, F. P., and Salzberg, S. L. (2016). Centrifuge: rapid and sensitive classification of metagenomic sequences. *Genome Research*, 26(12), 1721-1729. https://doi.org/10.1101/gr.210641.116

Knight, R., Vrbanac, A., Taylor, B. C., Aksenov, A., Callewaert, C., Debelius, J., Gonzalez, A., Kosciolek, T., McCall, L. I., McDonald, D., Melnik, A. V., Morton, J. T., Navas, J., Quinn, R. A., Sanders, J. G., Swafford, A. D., Thompson, L. R., Tripathi, A., Xu, Z. Z., Zaneveld, J. R., Zhu, Q., Caporaso, J. G., and Dorrestein, P. C. (2018). Best practices for analysing microbiomes. *Nature Reviews Microbiology*, 16(7), 410-422. https://doi.org/10.1038/s41579-018-0029-9

Lin, H. and Peddada, S. D. (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. *Nature Methods*, 21(1), 83-91. https://doi.org/10.1038/s41592-023-02092-7

Loubinoux, J., Bronowicki, J. P., Pereira, I. A., Mougenel, J. L., and Faou, A. E. (2002). Sulfate-reducing bacteria in human feces and their association with inflammatory bowel diseases. *FEMS Microbiology Ecology*, 40(2), 107-112. https://doi.org/10.1111/j.1574-6941.2002.tb00942.x

Lu, J., Breitwieser, F. P., Thielen, P., and Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104. https://doi.org/10.7717/peerj-cs.104

McMurdie, P. J. and Holmes, S. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLoS ONE*, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217

Nearing, J. T., Douglas, G. M., Hayes, M. G., MacDonald, J., Desai, D. K., Allward, N., Jones, C. M. A., Wright, R. J., Dhanani, A. S., Comeau, A. M., and Langille, M. G. I. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications*, 13(1), 342. https://doi.org/10.1038/s41467-022-28034-z

Oksanen, J., Simpson, G. L., Blanchet, F. G., Kindt, R., Legendre, P., Minchin, P. R., O'Hara, R. B., Solymos, P., Stevens, M. H. H., Szoecs, E., Wagner, H., Barbour, M., Bedward, M., Bolker, B., Borcard, D., Carvalho, G., Chirico, M., De Caceres, M., Durand, S., Evangelista, H., FitzJohn, R., Friendly, M., Furneaux, B., Hannigan, G., Hill, M. O., Lahti, L., McGlinn, D., Ouellette, M., Ribeiro Cunha, E., Smith, T., Stier, A., Ter Braak, C., and Weedon, J. (2024). vegan: Community Ecology Package. R package version 2.7-3. https://CRAN.R-project.org/package=vegan

Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P. A., and Hugenholtz, P. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*, 50(D1), D199-D207. https://doi.org/10.1093/nar/gkab776

Riviere, A., Selak, M., Lantin, D., Leroy, F., and De Vuyst, L. (2016). Bifidobacteria and butyrate-producing colon bacteria: importance and strategies for their stimulation in the human gut. *Frontiers in Microbiology*, 7, 979. https://doi.org/10.3389/fmicb.2016.00979

Rowan, F. E., Docherty, N. G., Murphy, M., Murphy, B., Calvin Coffey, J., and O'Connell, P. R. (2010). Desulfovibrio bacterial species are increased in ulcerative colitis. *Diseases of the Colon & Rectum*, 53(11), 1530-1536. https://doi.org/10.1007/DCR.0b013e3181f1e620


Turnbaugh, P. J., Ley, R. E., Mahowald, M. A., Magrini, V., Mardis, E. R., and Gordon, J. I. (2006). An obesity-associated gut microbiome with increased capacity for energy harvest. *Nature*, 444(7122), 1027-1031. https://doi.org/10.1038/nature05414

Waters, J. L. and Ley, R. E. (2019). The human gut bacteria Christensenellaceae are widespread, heritable, and associated with health. *BMC Biology*, 17, 83. https://doi.org/10.1186/s12915-019-0699-4

Wood, D. E., Lu, J., and Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257. https://doi.org/10.1186/s13059-019-1891-0

Wu, G. D., Chen, J., Hoffmann, C., Bittinger, K., Chen, Y. Y., Keilbaugh, S. A., Bewtra, M., Knights, D., Walters, W. A., Knight, R., Sinha, R., Gilroy, E., Gupta, K., Baldassano, R., Nessel, L., Li, H., Bushman, F. D., and Lewis, J. D. (2011). Linking long-term dietary patterns with gut microbial enterotypes. *Science*, 334(6052), 105-108. https://doi.org/10.1126/science.1208344

---
