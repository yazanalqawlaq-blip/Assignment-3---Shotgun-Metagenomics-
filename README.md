# Assignment 3: Shotgun Metagenomics Analysis of Vegan and Omnivore Gut Microbiomes
 
## Table of Contents
- [Introduction](#introduction)
- [Methods](#methods)
- [Results](#results)
- [Discussion](#discussion)
- [Limitations](#limitations)
- [Conclusion](#conclusion)
- [References](#references)
- [Project Structure](#project-structure)
 
---
 
## Introduction
 
The human gut microbiome is a complex community of microorganisms that plays a central role in host metabolism, immune function, and overall health (Turnbaugh et al., 2006). Among the many factors that shape gut microbial composition, diet is one of the most well-documented drivers. Long-term dietary patterns have been shown to influence both the taxonomic structure and functional capacity of the gut microbiome (Wu et al., 2011; David et al., 2014). In particular, plant-based diets tend to enrich for fiber-degrading taxa and short-chain fatty acid (SCFA) producers, while omnivorous diets that include red meat have been associated with protein-fermenting and bile-tolerant microorganisms (De Filippo et al., 2010; David et al., 2014). A recent large-scale study across 21,561 individuals found that microbial profiles could distinguish vegan, vegetarian, and omnivore diets with high accuracy, and that omnivore-associated taxa such as *Ruminococcus torques* and *Bilophila wadsworthia* correlated with less favourable cardiometabolic health markers (Fackelmann et al., 2025).
 
The present analysis uses shotgun metagenomic data from De Filippis et al. (2019), who sequenced the gut metagenomes of 97 healthy Italian adults following vegan (n=36), vegetarian (n=38), or omnivore (n=23) diets. Their study focused on strain-level diversity of *Prevotella copri* across dietary groups, and the raw sequencing data are publicly available on NCBI under accession SRP126540. For this assignment, I selected 4 vegan and 4 omnivore samples with comparable sequencing depth (~9-11 GB each) to compare gut microbial community composition, diversity, and differential abundance between the two dietary groups.
 
Shotgun metagenomics was the sequencing approach used in this study, and it offers several advantages over 16S rRNA amplicon sequencing (metabarcoding) for this type of analysis. While 16S sequencing is cheaper and benefits from more complete reference databases, it is limited to taxonomic profiling at genus-level resolution, is subject to PCR amplification biases, and cannot capture functional information directly from the reads (Knight et al., 2018). Shotgun metagenomics sequences all DNA in a sample, allowing for species- and strain-level resolution, detection of functional genes, and avoidance of primer bias. The trade-off is that shotgun approaches require deeper sequencing and larger reference databases for adequate classification.
 
For taxonomic classification, I used Kraken2 (Wood et al., 2019), a kmer-based classifier that assigns taxonomic labels by matching exact subsequences against a reference database. Kraken2 was chosen over alignment-based alternatives such as MetaPhlAn (Beghini et al., 2021), which relies on clade-specific marker genes and therefore has a smaller database footprint but can miss taxa not represented by markers. Centrifuge (Kim et al., 2016) was also considered, but Kraken2 is better established in the course material and is well-suited to paired-end short read data. One limitation of Kraken2 is its tendency to produce false positives at lower confidence thresholds, which I addressed by setting `--confidence 0.15` as recommended in the course lectures. Abundance re-estimation was performed with Bracken (Lu et al., 2017), which redistributes reads assigned at higher taxonomic levels back down to species level while accounting for genome size differences. For downstream analysis in R, I used phyloseq (McMurdie and Holmes, 2013) for data handling and visualization, the vegan package (Oksanen et al., 2024) for diversity metrics and PERMANOVA, ALDEx2 (Fernandes et al., 2014) for differential abundance analysis, and ANCOM-BC2 (Lin and Peddada, 2023) as a second differential abundance method. Both ALDEx2 and ANCOM-BC2 are designed to handle the compositionality problem inherent to metagenomic count data, where counts represent relative rather than absolute abundances and are therefore not independent (Gloor et al., 2017). Methods like DESeq2 or edgeR, which were originally designed for RNA-seq, have been shown to produce inflated false positive rates when applied to compositional microbiome data (Nearing et al., 2022), which is why I opted for tools specifically developed for this purpose.
 
---
 
## Methods
