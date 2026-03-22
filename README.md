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
 
## Results
 
### Quality Control
 
Quality filtering with fastp retained a high proportion of reads across all samples (Table 3). Seven of eight samples had pass rates above 97.7%, with Q20 scores exceeding 96% and Q30 above 91%. SRR8146976 had slightly lower quality (Q30 = 83.3%), likely due to a lower quality sequencing run, but still retained 54.9 million reads after filtering, which is sufficient for downstream classification.
 
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
 
Kraken2 classification rates were low across all samples, with the majority of reads remaining unclassified (>96% in all samples). For example, SRR8146974 had 34.9 million reads but only approximately 265,000 were classified at the species level. This is expected when using the Standard-8 database, which is a reduced version of the full Kraken2 standard database. After Bracken re-estimation, 225 species were identified across the 8 samples, with an average of 107 species per sample.
 
### Taxonomic Composition
 
At the phylum level, both vegan and omnivore gut microbiomes were dominated by Bacillota (formerly Firmicutes) and Bacteroidota (formerly Bacteroidetes), which together accounted for the majority of classified reads in all samples (Figure 1). Actinomycetota was visibly more abundant in vegan samples, driven primarily by *Bifidobacterium* species. Other phyla detected included Verrucomicrobiota (represented almost entirely by *Akkermansia muciniphila*), Pseudomonadota, and Thermodesulfobacteriota.
 
**Figure 1. Phylum-level taxonomic composition of vegan and omnivore gut microbiomes.** Stacked bar plot showing relative abundance of major phyla (>1% in any sample) in each sample, faceted by diet group. Both groups are dominated by Bacillota and Bacteroidota. Actinomycetota appears enriched in vegan samples.
 
![Phylum Barplot](output_files/R_analysis/phylum_barplot.png)
 
At the family level (Figure 2), composition differences between groups became more apparent. The Prevotellaceae family was prominent in several samples from both groups, though its distribution was highly variable across individuals. Bacteroidaceae, Lachnospiraceae, and Oscillospiraceae were consistently abundant in both groups. Bifidobacteriaceae was notably more abundant in the vegan samples.
 
**Figure 2. Family-level taxonomic composition.** Stacked bar plot showing relative abundance of the top 15 families, with remaining families grouped as "Other." Faceted by diet group. Notable differences include higher Bifidobacteriaceae in vegans and variable Prevotellaceae representation across individuals in both groups.
 
![Family Barplot](output_files/R_analysis/family_barplot.png)
 
### Alpha Diversity
 
Alpha diversity did not differ significantly between vegan and omnivore groups for any of the three metrics tested (Figure 3, Table 4). Observed species richness was nearly identical between groups (Vegan mean = 108.8, Omnivore mean = 105.0; Wilcoxon p = 1.00). Shannon diversity was slightly higher in vegans (3.01 vs 2.90; p = 0.49), and Simpson index was also comparable (0.89 vs 0.88; p = 0.69). None of these differences reached statistical significance.
 
**Figure 3. Alpha diversity comparison between vegan and omnivore groups.** Boxplots of Observed species richness, Shannon index, and Simpson index. Points represent individual samples. No significant differences were detected for any metric (Wilcoxon rank-sum test, all p > 0.05).
 
![Alpha Diversity](output_files/R_analysis/alpha_diversity.png)
 
**Table 4. Alpha diversity summary statistics.**
 
| Metric | Vegan Mean | Omnivore Mean | Wilcoxon p-value |
|:---|:---|:---|:---|
| Observed | 108.8 | 105.0 | 1.000 |
| Shannon | 3.01 | 2.90 | 0.486 |
| Simpson | 0.891 | 0.876 | 0.686 |
 
### Beta Diversity
 
PCoA ordination based on Bray-Curtis dissimilarity showed partial separation between vegan and omnivore samples, though with considerable overlap (Figure 4). NMDS ordination produced a similar pattern with a stress value of 0.048, indicating an excellent fit (stress < 0.1 is generally considered good). PERMANOVA revealed that diet explained 21.2% of the variation in community composition (R² = 0.212, F = 1.61), though this result was not statistically significant at the conventional threshold (p = 0.108, 999 permutations; Table 5).
 
**Figure 4. PCoA ordination of gut microbial communities based on Bray-Curtis dissimilarity.** Points represent individual samples, colored by diet group. Ellipses show 95% confidence intervals. Partial separation is visible between groups. PERMANOVA: R² = 0.212, F = 1.61, p = 0.108.
 
![PCoA](output_files/R_analysis/pcoa_bray.png)
 
**Figure 5. NMDS ordination based on Bray-Curtis dissimilarity.** Stress = 0.048, indicating excellent ordination fit. Diet groups show partial separation consistent with the PCoA result.
 
![NMDS](output_files/R_analysis/nmds_bray.png)
 
**Table 5. PERMANOVA results (adonis2, 999 permutations).**
 
| Source | Df | Sum of Squares | R² | F | p-value |
|:---|:---|:---|:---|:---|:---|
| Diet | 1 | 0.354 | 0.212 | 1.611 | 0.108 |
| Residual | 6 | 1.319 | 0.788 | | |
| Total | 7 | 1.673 | 1.000 | | |
 
### Differential Abundance
 
Neither ALDEx2 nor ANCOM-BC2 identified any taxa as significantly differentially abundant between vegan and omnivore groups after multiple testing correction (ALDEx2: all BH-adjusted p > 0.05; ANCOM-BC2: all Holm-adjusted q > 0.05). Both tools warned that the small sample size (n = 4 per group) limits statistical power and makes variance estimation unstable.
 
However, ANCOM-BC2 detected 15 structural zeroes, which are taxa present in one group but completely absent from the other (Table 6). Several of these have well-documented associations with diet. Taxa found exclusively in omnivore samples included *Escherichia*, *Desulfovibrio*, *Streptococcus*, and *Cutibacterium*. Taxa exclusive to vegan samples included *Christensenella*, *Eisenbergiella*, *Megamonas*, and *Wansuia*.
 
**Table 6. Selected structural zeroes detected by ANCOM-BC2.** These taxa were present in one diet group but absent from all samples in the other.
 
| Taxon | Present in Omnivore | Present in Vegan |
|:---|:---|:---|
| *Escherichia* | Yes | No |
| *Desulfovibrio* | Yes | No |
| *Streptococcus* | Yes | No |
| *Cutibacterium* | Yes | No |
| *Oxalobacter* | Yes | No |
| *Haemophilus* | No | Yes |
| *Christensenella* | No | Yes |
| *Eisenbergiella* | No | Yes |
| *Megamonas* | No | Yes |
| *Wansuia* | No | Yes |
 
ALDEx2 effect size analysis showed trends in taxa abundance between groups even though none reached significance. Among the genera with the largest effect sizes, *Phocaeicola* (effect = -0.65) and *Parabacteroides* (effect = -0.63) trended higher in omnivores, while *Bifidobacterium* showed the strongest trend toward enrichment in vegans (ANCOM-BC2 log fold change = 2.98, q = 0.68).
 
**Figure 6. ANCOM-BC2 differential abundance results.** Lollipop plot of log fold change (Vegan vs Omnivore) for genera with q < 0.25. Error bars represent standard error. No taxa reached q < 0.05 (all shown in grey). Positive values indicate higher abundance in vegans.
 
![ANCOM-BC2 Lollipop](output_files/R_analysis/ancombc2_lollipop.png)
 
**Figure 7. ALDEx2 effect size plot.** Top 10 genera with the largest positive and negative effect sizes between diet groups. Positive effect indicates higher abundance in omnivores. None reached BH-corrected significance at p < 0.05.
 
![ALDEx2 Effect](output_files/R_analysis/aldex2_effect.png)
 
---
 
## Discussion
 
Although no taxa reached statistical significance in the differential abundance analysis, the compositional trends and structural zeroes observed in this analysis are consistent with established patterns in the diet-microbiome literature. The lack of statistical significance is most likely attributable to the small sample size (n = 4 per group), which limits power for both diversity comparisons and differential abundance testing. ANCOM-BC2 itself flagged that variance estimation is unstable with fewer than 5 samples per group.
 
The structural zero analysis provided the most interpretable results. *Desulfovibrio*, which was present in omnivore samples but absent from all vegan samples, is a sulfate-reducing bacterium that produces hydrogen sulfide (H2S). H2S is cytotoxic to colonic epithelial cells and has been linked to gut inflammation (Loubinoux et al., 2002). The enrichment of *Desulfovibrio* in omnivores is consistent with previous findings that animal-based diets promote sulfate-reducing bacteria due to higher intake of sulfur-containing amino acids from meat (David et al., 2014). Similarly, *Bilophila wadsworthia*, a sulfite-reducing bacterium found in all samples but at higher counts in several omnivore samples, has been shown to expand in response to dietary saturated fat through increased taurine-conjugated bile acid availability (Devkota et al., 2012).
 
*Escherichia* was detected exclusively in omnivore samples. While *Escherichia* includes commensal strains, certain strains are associated with gut inflammation, and its increased abundance in omnivore microbiomes has been reported in studies comparing dietary patterns (Fackelmann et al., 2025).
 
On the vegan side, *Christensenella* was found only in vegan samples and was absent from all omnivore samples. This genus is of particular interest because it has been associated with lean body mass and is one of the most heritable taxa in the human gut (Goodrich et al., 2014). Its presence in vegans is consistent with the broader association between plant-based diets and lean-associated microbial profiles. *Megamonas*, also exclusive to vegans in this analysis, is a member of the Selenomonadaceae family and is involved in carbohydrate fermentation and SCFA production.
 
The trend toward higher *Bifidobacterium* abundance in vegan samples (ANCOM-BC2 lfc = 2.98) is consistent with previous findings. *Bifidobacterium* species are well-known fiber fermenters that produce acetate and lactate, and they are commonly enriched in individuals consuming plant-based diets rich in dietary fiber and complex carbohydrates (De Filippo et al., 2010). The enrichment of *Bifidobacterium adolescentis* in particular has been reported in vegan populations (Wu et al., 2011).
 
The PERMANOVA result, while not reaching the p < 0.05 threshold (p = 0.108), did show that diet explained 21.2% of the variation in community composition. This is a meaningful effect size for a microbiome study, where individual variation is typically high and many factors beyond diet contribute to microbial community structure. A study with more samples would likely detect this as significant. Alpha diversity was comparable between groups, which aligns with the findings of Fackelmann et al. (2025), who reported that omnivores had slightly greater richness than vegans across larger cohorts, though the differences were small.
 
The De Filippis et al. (2019) source paper focused on *Prevotella copri* (now reclassified as *Segatella copri*) strain-level diversity rather than community-level analysis. In my data, *Segatella copri* was present in samples from both groups, with particularly high read counts in some individuals (e.g., 166,246 reads in SRR8146961, a vegan, and 169,596 in SRR8146969, an omnivore). This high inter-individual variability in *S. copri* abundance is consistent with the original paper's observation that *P. copri* presence was not significantly different between dietary groups at the species level, but that functionally distinct strains were associated with different diets. Species-level shotgun analysis, as performed here, cannot capture this strain-level variation.
 
---
 
## Limitations
 
Several limitations should be noted when interpreting these results:
 
1. **Sample size.** With only 4 samples per group, the statistical power to detect differential abundance or significant PERMANOVA results is limited. Both ALDEx2 and ANCOM-BC2 warned about unstable variance estimates at this sample size.
 
2. **Database completeness.** The Standard-8 Kraken2 database is a reduced database designed to fit in 8 GB of RAM. The classification rate was very low (<4% of reads classified in most samples). A larger database such as the full Kraken2 Standard database or GTDB would classify more reads but requires substantially more memory and storage.
 
3. **Compositionality.** Metagenomic count data is compositional, meaning that an increase in one taxon's relative abundance necessarily decreases others. Both ALDEx2 and ANCOM-BC2 are designed to handle this, but it remains a fundamental constraint of the data (Gloor et al., 2017).
 
4. **No functional analysis.** This analysis is limited to taxonomic composition. Functional profiling with tools such as HUMAnN (Beghini et al., 2021) or analysis of functional gene content could reveal metabolic differences between groups that are not apparent from taxonomy alone.
 
5. **Geographic bias.** All omnivore samples are from Turin, while the vegan samples span Turin, Bari, and Parma. Regional dietary differences beyond the vegan/omnivore distinction could contribute to the observed variation.
 
6. **Host DNA contamination.** *Homo sapiens* reads were detected in all samples at low levels. While the Kraken2 database includes human sequences and these reads are classified rather than lost, they do occupy a small fraction of the abundance table.
 
---
 
## Conclusion
 
This analysis compared the gut microbial communities of vegan and omnivore individuals from the De Filippis et al. (2019) cohort using a shotgun metagenomics pipeline consisting of fastp for quality control, Kraken2 and Bracken for taxonomic classification and abundance estimation, and phyloseq, ALDEx2, and ANCOM-BC2 for statistical analysis in R. While no taxa reached statistical significance after multiple testing correction, the structural zero analysis and effect size trends point toward biologically meaningful differences. Omnivore-specific taxa such as *Desulfovibrio* and *Escherichia* are associated with sulfate reduction and inflammation, while vegan-specific taxa such as *Christensenella* and *Megamonas* are linked to lean phenotypes and carbohydrate fermentation. These findings are broadly consistent with the published literature on diet-microbiome interactions, though a larger sample size and a more comprehensive classification database would be needed to detect statistically significant differential abundance.
 
---
 
## References
 
Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32-46. https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x
 
Beghini, F., McIver, L. J., Blanco-Miguez, A., Dubois, L., Asnicar, F., Maharjan, S., Mailyan, A., Manghi, P., Scholz, M., Thomas, A. M., Valles-Colomer, M., Weingart, G., Zhang, Y., Zolfo, M., Huttenhower, C., Franzosa, E. A., and Segata, N. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife*, 10, e65088. https://doi.org/10.7554/eLife.65088
 
Cassol, I., Ibarra, A., Singh, J., Subedi, S., Zeballos, J., Kaviani, S., and Nearing, J. T. (2024). Comprehensive evaluation of taxonomic diversity metrics for microbiome data. *Scientific Reports*, 14, 25230. https://doi.org/10.1038/s41598-024-77864-y
 
Chen, S., Zhou, Y., Chen, Y., and Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560
 
Dabdoub, S. (2016). kraken-biom: Enabling interoperative format conversion for Kraken results. https://github.com/smdabdoub/kraken-biom
 
David, L. A., Maurice, C. F., Carmody, R. N., Gootenberg, D. B., Button, J. E., Wolfe, B. E., Ling, A. V., Devlin, A. S., Varma, Y., Fischbach, M. A., Biddinger, S. B., Dutton, R. J., and Turnbaugh, P. J. (2014). Diet rapidly and reproducibly alters the human gut microbiome. *Nature*, 505(7484), 559-563. https://doi.org/10.1038/nature12820
 
De Filippo, C., Cavalieri, D., Di Paola, M., Ramazzotti, M., Poullet, J. B., Massart, S., Collini, S., Pieraccini, G., and Lionetti, P. (2010). Impact of diet in shaping gut microbiota revealed by a comparative study in children from Europe and rural Africa. *Proceedings of the National Academy of Sciences*, 107(33), 14691-14696. https://doi.org/10.1073/pnas.1005963107
 
De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., and Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets. *Cell Host & Microbe*, 25(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004
 
Devkota, S., Wang, Y., Musch, M. W., Leone, V., Fehlner-Peach, H., Nadimpalli, A., Antonopoulos, D. A., Bhatt, B., and Chang, E. B. (2012). Dietary-fat-induced taurocholic acid promotes pathobiont expansion and colitis in Il10-/- mice. *Nature*, 487(7405), 104-108. https://doi.org/10.1038/nature11225
 
Fackelmann, G., Gimenez-Bastida, J. A., Tett, A., Segata, N., and Acharjee, A. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. *Nature Microbiology*, 10, 42-52. https://doi.org/10.1038/s41564-024-01870-z
 
Fernandes, A. D., Reid, J. N., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., and Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome*, 2, 15. https://doi.org/10.1186/2049-2618-2-15
 
Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., and Egozcue, J. J. (2017). Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology*, 8, 2224. https://doi.org/10.3389/fmicb.2017.02224
 
Goodrich, J. K., Waters, J. L., Poole, A. C., Sutter, J. L., Koren, O., Blekhman, R., Beaumont, M., Van Treuren, W., Knight, R., Bell, J. T., Spector, T. D., Clark, A. G., and Ley, R. E. (2014). Human genetics shape the gut microbiome. *Cell*, 159(4), 789-799. https://doi.org/10.1016/j.cell.2014.09.053
 
Kim, D., Song, L., Breitwieser, F. P., and Salzberg, S. L. (2016). Centrifuge: rapid and sensitive classification of metagenomic sequences. *Genome Research*, 26(12), 1721-1729. https://doi.org/10.1101/gr.210641.116
 
Knight, R., Vrbanac, A., Taylor, B. C., Aksenov, A., Callewaert, C., Debelius, J., Gonzalez, A., Kosciolek, T., McCall, L. I., McDonald, D., Melnik, A. V., Morton, J. T., Navas, J., Quinn, R. A., Sanders, J. G., Swafford, A. D., Thompson, L. R., Tripathi, A., Xu, Z. Z., Zaneveld, J. R., Zhu, Q., Caporaso, J. G., and Dorrestein, P. C. (2018). Best practices for analysing microbiomes. *Nature Reviews Microbiology*, 16(7), 410-422. https://doi.org/10.1038/s41579-018-0029-9
 
Lin, H. and Peddada, S. D. (2023). Analysis of compositions of microbiomes with bias correction. *Nature Methods*, 20(2), 188-191. https://doi.org/10.1038/s41592-023-02092-7
 
Loubinoux, J., Bronowicki, J. P., Pereira, I. A., Mougenel, J. L., and Faou, A. E. (2002). Sulfate-reducing bacteria in human feces and their association with inflammatory bowel diseases. *FEMS Microbiology Ecology*, 40(2), 107-112. https://doi.org/10.1111/j.1574-6941.2002.tb00942.x
 
Lu, J., Breitwieser, F. P., Thielen, P., and Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104. https://doi.org/10.7717/peerj-cs.104
 
McMurdie, P. J. and Holmes, S. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLoS ONE*, 8(4), e61217. https://doi.org/10.1371/journal.pone.0061217
 
Nearing, J. T., Douglas, G. M., Hayes, M. G., MacDonald, J., Desai, D. K., Allward, N., Jones, C. M. A., Wright, R. J., Dhanani, A. S., Comeau, A. M., and Langille, M. G. I. (2022). Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications*, 13(1), 342. https://doi.org/10.1038/s41467-022-28034-z
 
Oksanen, J., Simpson, G. L., Blanchet, F. G., Kindt, R., Legendre, P., Minchin, P. R., O'Hara, R. B., Solymos, P., Stevens, M. H. H., Szoecs, E., Wagner, H., Barbour, M., Bedward, M., Bolker, B., Borcard, D., Carvalho, G., Chirico, M., De Caceres, M., Durand, S., Evangelista, H., FitzJohn, R., Friendly, M., Furneaux, B., Hannigan, G., Hill, M. O., Lahti, L., McGlinn, D., Ouellette, M., Ribeiro Cunha, E., Smith, T., Stier, A., Ter Braak, C., and Weedon, J. (2024). vegan: Community Ecology Package. R package version 2.7-3. https://CRAN.R-project.org/package=vegan
 
R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
 
Turnbaugh, P. J., Ley, R. E., Mahowald, M. A., Magrini, V., Mardis, E. R., and Gordon, J. I. (2006). An obesity-associated gut microbiome with increased capacity for energy harvest. *Nature*, 444(7122), 1027-1031. https://doi.org/10.1038/nature05414
 
Wood, D. E., Lu, J., and Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20(1), 257. https://doi.org/10.1186/s13059-019-1891-0
 
Wu, G. D., Chen, J., Hoffmann, C., Bittinger, K., Chen, Y. Y., Keilbaugh, S. A., Bewtra, M., Knights, D., Walters, W. A., Knight, R., Sinha, R., Gilroy, E., Gupta, K., Baldassano, R., Nessel, L., Li, H., Bushman, F. D., and Lewis, J. D. (2011). Linking long-term dietary patterns with gut microbial enterotypes. *Science*, 334(6052), 105-108. https://doi.org/10.1126/science.1208344
 
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
├── input_data/
│   ├── raw_reads/          # gitignored (large files)
│   └── kraken2_db/         # gitignored (large files)
├── output_files/
│   ├── qc/                 # fastp reports
│   ├── kraken2/            # classification reports
│   ├── bracken/            # abundance estimates
│   └── R_analysis/         # BIOM file, figures, stats
└── slurm/                  # job logs (gitignored)
```
 
