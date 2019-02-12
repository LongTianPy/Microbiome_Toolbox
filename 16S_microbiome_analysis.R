# Title     : 16s Microbiome analysis data visualization
# Objective : Data visualization
# Created by: longtian
# Created on: 2/12/19

# Import required packages here
source('http://bioconductor.org/biocLite.R')
# Upgrade BioConductor when necessary
# biocLite("BiocUpgrade")
biocLite('phyloseq')
library("phyloseq")
library("scales")
library("ggplot2")

# Set working directory and import files
# Working directory, replace with where most of the input files are
setwd("/Users/longtian/Desktop/VinatzerLab/Dianthus")
# Import mapping file
##########################################
# Make sure the mapping file has 8 columns
##########################################
map <- import_qiime_sample_data("Fasting_mapVTP64.txt")
# Import OTU table (BIOM file)
otu <- import_biom(BIOMfilename="otu_table_even23675.biom")
# Merge mapping file and otu table into a phyloseq object
run <- merge_phyloseq(otu,map)
# Rename tax ranks with the names we are familiar with
colnames(tax_table(run)) <- c("Kingdom", "Phylum", "Class",
                                    "Order", "Family", "Genus","Species","Rank8","Rank9",
                                    "Rank10","Rank11","Rank12","Rank13","Rank14","Rank15")
# Abundance barplot
library(dplyr)
# Select those genus with abundance over 2%
# Can be used to filter specific taxonomic rank at any abundance cutoff
run_genus <- run %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  psmelt() %>%
  filter(Abundance > 0.02) %>%
  arrange(Genus)
# Import color palette
library(RColorBrewer)
# Use distinctive colors
n <- dim(run_genus)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# Plot barplot
ggplot(run_genus,aes(x=X.SampleID,y=Abundance,fill=Genus)) +
  geom_bar(position="fill",stat="identity") +
  scale_fill_manual(values = col_vector) +
  guides(fill=guide_legend(reverse=T,keywidth = 1,keyheight = 1)) +
  ylab("Relative Abundance (Genera > 2%) \n") +
  xlab("Sample") +
  # ggtitle("Genus Composition")+
  scale_y_continuous(labels=percent_format(),expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90,hjust=1), panel.background = element_blank())

# Alpha diversity
# Calculate rarefaction curve data
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures) # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none')) # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}
# Remember to change the rarefaction steps acoording to the depth
rarefaction_curve_data <- calculate_rarefaction_curves(run, c('Observed'), rep(seq(1,23675,by=1000), each = 10))
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(run)), by.x = 'Sample', by.y = 'row.names')
rarefaction_curve_data_summary_verbose$X.SampleID <- factor(rarefaction_curve_data_summary_verbose$X.SampleID,levels = map$X.SampleID)
n <- dim(rarefaction_curve_data_summary_verbose)[1]
display.brewer.all()
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot( data = rarefaction_curve_data_summary_verbose,
        mapping = aes( x = Depth, y = Alpha_diversity_mean,
                       ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                       ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                       colour = X.SampleID,
                       group = X.SampleID)
        ) +
  scale_fill_manual(values=col_vector) +
  geom_line( ) +
  geom_point(size=0.5 ) +
  guides(fill=guide_legend(title="Sample")) +
  xlab("Sequences per sample") +
  ylab("Observed OTU")

ggsave('rain_VTP064_alpha_diversity.pdf',device = 'pdf',units = 'mm',width = 297,height = 210,dpi = 300)



