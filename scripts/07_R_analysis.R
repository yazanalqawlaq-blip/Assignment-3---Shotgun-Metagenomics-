## 07_R_analysis.R
## Diversity analysis and differential abundance of vegan vs omnivore gut metagenomes
## Dataset: De Filippis et al. (2019) Cell Host & Microbe - SRP126540
## Input: combined_bracken.biom (from kraken-biom)
## Output: taxonomy barplots, alpha/beta diversity plots, ALDEx2 and ANCOM-BC2 results

# ---- Load packages ----
library(phyloseq)
library(biomformat)
library(vegan)
library(ALDEx2)
library(ANCOMBC)
library(ggplot2)

# ---- Set paths ----
base_dir <- "C:/Users/yazan/Desktop/Assignment-3---Shotgun-Metagenomics-"
biom_path <- file.path(base_dir, "output_files", "R_analysis", "combined_bracken.biom")
output_dir <- file.path(base_dir, "output_files", "R_analysis")

# ---- Import BIOM into phyloseq ----
biom_data <- read_biom(biom_path)
physeq <- import_biom(biom_data)

physeq
head(tax_table(physeq))
sample_names(physeq)

# fix sample names, kraken-biom appends "_bracken_species"
sample_names(physeq) <- gsub("_bracken_species", "", sample_names(physeq))

# ---- Add metadata ----
sample_df <- data.frame(
    SampleID = c("SRR8146974", "SRR8146973", "SRR8146965", "SRR8146961",
                 "SRR8146969", "SRR8146975", "SRR8146970", "SRR8146976"),
    Diet = c("Vegan", "Vegan", "Vegan", "Vegan",
             "Omnivore", "Omnivore", "Omnivore", "Omnivore"),
    City = c("Turin", "Turin", "Bari", "Parma",
             "Turin", "Turin", "Turin", "Turin"),
    SubjectID = c("VOV113", "VOV114", "11BA", "17PR",
                  "VOV39", "VOV77", "VOV36", "VOV70")
)
rownames(sample_df) <- sample_df$SampleID
sample_data(physeq) <- sample_data(sample_df)

# ---- Fix taxonomy column names ----
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "Species")
head(tax_table(physeq))

# ---- Rarefaction curves ----
otu_df <- as.data.frame(t(otu_table(physeq)))
rarecurve(otu_df, step = 50, label = TRUE,
          main = "Rarefaction Curves", xlab = "Reads", ylab = "Species")

# ---- Taxonomic abundance barplot (Phylum level) ----
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))
physeq_phy <- tax_glom(physeq_rel, taxrank = "Phylum")

df_phy <- psmelt(physeq_phy)

major_phyla <- unique(df_phy$Phylum[df_phy$Abundance >= 0.01])
df_phy$Phylum_plot <- ifelse(df_phy$Phylum %in% major_phyla, df_phy$Phylum, "Other")

ggplot(df_phy, aes(x = Sample, y = Abundance, fill = Phylum_plot)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Diet, scales = "free_x") +
    labs(y = "Relative Abundance", x = "Sample", fill = "Phylum",
         title = "Phylum-Level Taxonomic Composition") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# ---- Taxonomic abundance barplot (Family level) ----
physeq_fam <- tax_glom(physeq_rel, taxrank = "Family")
df_fam <- psmelt(physeq_fam)

top_fam <- names(sort(tapply(df_fam$Abundance, df_fam$Family, mean),
                      decreasing = TRUE))[1:15]
df_fam$Family_plot <- ifelse(df_fam$Family %in% top_fam, df_fam$Family, "Other")

ggplot(df_fam, aes(x = Sample, y = Abundance, fill = Family_plot)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Diet, scales = "free_x") +
    labs(y = "Relative Abundance", x = "Sample", fill = "Family",
         title = "Family-Level Taxonomic Composition") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.text = element_text(size = 7))

# ---- Alpha diversity ----
plot_richness(physeq, x = "Diet", color = "Diet",
              measures = c("Observed", "Shannon", "Simpson")) +
    geom_boxplot(alpha = 0.3) +
    theme_bw() +
    labs(title = "Alpha Diversity by Diet Group")

alpha_div <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Diet <- sample_data(physeq)$Diet

wilcox_observed <- wilcox.test(Observed ~ Diet, data = alpha_div)
wilcox_shannon <- wilcox.test(Shannon ~ Diet, data = alpha_div)
wilcox_simpson <- wilcox.test(Simpson ~ Diet, data = alpha_div)

cat("Wilcoxon tests (Vegan vs Omnivore):\n")
cat("Observed:", wilcox_observed$p.value, "\n")
cat("Shannon:", wilcox_shannon$p.value, "\n")
cat("Simpson:", wilcox_simpson$p.value, "\n")

alpha_stats <- data.frame(
    Metric = c("Observed", "Shannon", "Simpson"),
    Vegan_mean = c(mean(alpha_div$Observed[alpha_div$Diet == "Vegan"]),
                   mean(alpha_div$Shannon[alpha_div$Diet == "Vegan"]),
                   mean(alpha_div$Simpson[alpha_div$Diet == "Vegan"])),
    Omnivore_mean = c(mean(alpha_div$Observed[alpha_div$Diet == "Omnivore"]),
                      mean(alpha_div$Shannon[alpha_div$Diet == "Omnivore"]),
                      mean(alpha_div$Simpson[alpha_div$Diet == "Omnivore"])),
    Wilcoxon_p = c(wilcox_observed$p.value, wilcox_shannon$p.value, wilcox_simpson$p.value)
)
write.csv(alpha_stats, file.path(output_dir, "alpha_diversity_stats.csv"), row.names = FALSE)

# ---- Beta diversity: PCoA with Bray-Curtis ----
ord_pcoa_bray <- ordinate(physeq, method = "PCoA", distance = "bray")

plot_ordination(physeq, ord_pcoa_bray, color = "Diet") +
    geom_point(size = 4) +
    stat_ellipse(type = "t", level = 0.95) +
    theme_bw() +
    labs(title = "PCoA (Bray-Curtis Dissimilarity)")

# ---- Beta diversity: NMDS with Bray-Curtis ----
ord_nmds_bray <- ordinate(physeq, method = "NMDS", distance = "bray")

plot_ordination(physeq, ord_nmds_bray, color = "Diet") +
    geom_point(size = 4) +
    stat_ellipse(type = "t", level = 0.95) +
    theme_bw() +
    labs(title = "NMDS (Bray-Curtis Dissimilarity)")

# ---- PERMANOVA ----
metadata <- as(sample_data(physeq), "data.frame")
permanova_result <- adonis2(phyloseq::distance(physeq, method = "bray") ~ Diet,
                            data = metadata, permutations = 999)
cat("\nPERMANOVA result:\n")
print(permanova_result)

capture.output(permanova_result, file = file.path(output_dir, "permanova_result.txt"))

# ---- ALDEx2: Differential abundance at Genus level ----
physeq_genus <- tax_glom(physeq, taxrank = "Genus")

genus_counts <- as.data.frame(otu_table(physeq_genus))
if (!taxa_are_rows(physeq_genus)) {
    genus_counts <- t(genus_counts)
}

tax_info <- as.data.frame(tax_table(physeq_genus))
rownames(genus_counts) <- tax_info$Genus

conditions <- as.character(sample_data(physeq_genus)$Diet)

aldex_results <- aldex(genus_counts, conditions, mc.samples = 128,
                       test = "t", effect = TRUE, denom = "all")

aldex_sig <- subset(aldex_results, we.eBH < 0.05)
cat("\nALDEx2 significant taxa (BH-corrected p < 0.05):\n")
print(aldex_sig[, c("diff.btw", "diff.win", "effect", "we.eBH")])

write.csv(aldex_results, file.path(output_dir, "aldex2_results.csv"))
write.csv(aldex_sig, file.path(output_dir, "aldex2_significant.csv"))

# ALDEx2 MA plot
aldex_results$taxon <- rownames(aldex_results)
aldex_results$significant <- aldex_results$we.eBH < 0.05

aldex.plot(aldex_results, type = "MW", test = "welch")

# ALDEx2 effect size plot (top 10 each direction)
aldex_plot_df <- aldex_results[order(aldex_results$effect), ]
aldex_plot_df <- aldex_plot_df[c(1:10, (nrow(aldex_plot_df)-9):nrow(aldex_plot_df)), ]

ggplot(aldex_plot_df, aes(x = effect, y = reorder(taxon, effect))) +
    geom_point(aes(color = significant), size = 3) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    scale_color_manual(values = c("grey50", "steelblue"),
                       labels = c("Not significant", "BH p < 0.05"),
                       name = "Significance") +
    labs(x = "Effect Size (Omnivore vs Vegan)", y = "Genus",
         title = "Differential Abundance (ALDEx2)") +
    theme_bw()

# ---- ANCOM-BC2: Differential abundance at Genus level ----
ancombc_out <- ancombc2(data = physeq, tax_level = "Genus",
                        fix_formula = "Diet", rand_formula = NULL,
                        p_adj_method = "holm", pseudo_sens = TRUE,
                        prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                        group = "Diet", struc_zero = TRUE, neg_lb = TRUE)

cat("\nStructural zeroes:\n")
struc_zeros <- ancombc_out$zero_ind
struc_diff <- subset(struc_zeros,
                     `structural_zero (Diet = Omnivore)` != `structural_zero (Diet = Vegan)`)
print(struc_diff)

ancombc_res <- ancombc_out$res
ancombc_sig <- subset(ancombc_res, q_DietVegan < 0.05)
cat("\nANCOM-BC2 significant taxa (q < 0.05):\n")
print(ancombc_sig[, c("taxon", "lfc_DietVegan", "q_DietVegan")])

write.csv(ancombc_res, file.path(output_dir, "ancombc2_results.csv"), row.names = FALSE)
write.csv(ancombc_sig, file.path(output_dir, "ancombc2_significant.csv"), row.names = FALSE)

# ANCOM-BC2 plot
ancombc_plot_df <- ancombc_res[order(ancombc_res$q_DietVegan), ][1:20, ]

ggplot(ancombc_plot_df,
       aes(x = lfc_DietVegan, y = reorder(taxon, lfc_DietVegan))) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(xmin = lfc_DietVegan - se_DietVegan,
                    xmax = lfc_DietVegan + se_DietVegan), width = 0.3) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Log Fold Change (Vegan vs Omnivore)", y = "Genus",
       title = "ANCOM-BC2 Differential Abundance") +
  theme_bw()
