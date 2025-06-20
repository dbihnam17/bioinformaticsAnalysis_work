### This script is a visual QC check for the quality of our VCF files returned
### We had a batch of sequencing with way more called variants than expected
### and determined there to be a significant difference in the distribution of
### AD and VAF compared to previous batches

## Load libraries --------------------------------------------------------------
library(ggplot2)
library(vcfR)

## Data wrangling --------------------------------------------------------------
vcf <- read.vcfR('/path/to/sample.muTect2.somatic.snv.vcf.gz')

# Extract the genotype matrix
gt <- extract.gt(vcf, element = "GT")

# Extract the AF (allele frequency) values
af <- extract.gt(vcf, element = "AF", as.numeric = TRUE)

# Extract the AD (allelic depth) values
ad <- extract.gt(vcf, element = "AD")

# Split each value in the matrix by the comma
ad_split <- strsplit(ad, ",")

# Extract the second element from each split result, which corresponds to the variant allele depth
variant_reads <- sapply(ad_split, function(x) as.numeric(x[2]))

# Reshape the result into the original matrix shape (if necessary)
variant_reads_matrix <- matrix(variant_reads, nrow = nrow(ad), ncol = ncol(ad), byrow = FALSE)

# Calculate the total number of variants (rows) for normalization
total_variants <- nrow(variant_reads_matrix)

# Calculate the frequency of each unique alternative read count
alt_read_counts <- table(as.vector(variant_reads_matrix))

# Calculate the proportion of each alternative read count
alt_read_proportions <- alt_read_counts / total_variants

# Create a data frame for plotting
plot_data <- data.frame(
  AltReads = as.numeric(names(alt_read_proportions)),
  Proportion = as.numeric(alt_read_proportions)
)

# Extract VAF values from AF
vaf_values <- af[, 1]

# Calculate the proportion of variants for each VAF value
vaf_proportions <- table(vaf_values) / length(vaf_values)

# Create a data frame for plotting
vaf_data <- data.frame(VAF = as.numeric(names(vaf_proportions)),
                       Proportion = as.numeric(vaf_proportions))

## Plotting --------------------------------------------------------------------
# Plot AD frequency as barplot
ggplot(plot_data, aes(x = AltReads, y = Proportion)) +
  geom_histogram(stat = "identity", binwidth = 1, fill = "red") +
  labs(title = "<sample>: Alternative Read Count Proportion",
       x = "Number of Alternative Reads",
       y = "Proportion of Mutations") +
  xlim(0, 30) +
  ylim(0, max(plot_data$Proportion)) +
  theme(plot.title = element_text(size = 8))

# Plot the VAF frequency as a line
ggplot(vaf_data, aes(x = VAF, y = Proportion)) +
  geom_line(color = "red") +
  labs(title = "<sample>: Variant Allelic Frequency (VAF) Proportion",
       x = "VAF",
       y = "Proportion of Variants") +
  xlim(0, 1) +
  ylim(0, max(vaf_data$Proportion) * 1.1) +
  theme(plot.title = element_text(size = 8))


