#########################
####### LIBRARIES #######
#########################

library("optparse")
library("stringr")
library("data.table")
library("riboWaltz")


##########################
####### PARAMETERS #######
##########################

option_list = list(
  make_option(c("-w", "--work_dir"), type="character", 
              help="Path to working directory"),
  make_option(c("-f", "--feature_type"), type="character", 
              help="Path to designed R functions"),
  make_option(c("-g", "--gtf"), type="character", 
              help="Path and name of the input GTF file"),
  make_option(c("-p", "--path_prefix"), type="character", 
              help="Prefix for RESULTS path") 
)

opt = parse_args(OptionParser(option_list=option_list))

# Path to working directory
local_path <- opt$w

# Path prefix for RESULTS_ALL
results_all_path <- opt$p  # New argument for the base path prefix

# Path to remastered psite() function from riboWaltz
function_psite <- paste0(opt$f, "ribowaltz_psite_with_NA_control.R")

# Name of the gtf file
gtf_file <- opt$g

# Read the config file
params <- scan(file = paste0(local_path, "config.yaml"),
               what = "character",
               sep = ":"
)

readsLength_min <- gsub(" ", "", params[which(params == "readsLength_min") + 1], fixed = TRUE)
readsLength_max <- gsub(" ", "", params[which(params == "readsLength_max") + 1], fixed = TRUE)

# Create directories for riboWaltz and BAM results using the path_prefix
ribowaltz_folder <- paste0(results_all_path, "/riboWaltz.", readsLength_min, "-", readsLength_max, "/")
dir.create(ribowaltz_folder, recursive = TRUE, showWarnings = FALSE)

bam_folder <- paste0(results_all_path, "/BAM_transcriptome.", readsLength_min, "-", readsLength_max, "/")


#########################
####### EXECUTION #######
#########################

# Creates annotation table by transcript names
annotation_db <- riboWaltz::create_annotation(gtf_file)
annotation_db_transcript_with_cds0l <- data.table(annotation_db)
annotation_db_transcript <- annotation_db_transcript_with_cds0l[annotation_db_transcript_with_cds0l$l_cds > 0, ]
# Free unused memory
rm(list = c("annotation_db", "annotation_db_transcript_with_cds0l"))
gc()


# Bam files to be computed
bam_list <- list.files(bam_folder, pattern = "\\.bam$")

samples <- str_replace(bam_list, "\\.[0-9]{1,3}-[0-9]{1,3}\\.bam", "")
names(samples) <- str_remove(bam_list, "\\.bam")

# List of reads coordinates
reads_list <- try(riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples))
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)

# Changing dashes ("-") in names to underscores ("_") as ggplot uses parse which throws an error with digits followed by dashes
samples_renamed <- gsub("-", "_", samples)
names(samples_renamed) <- gsub("-", "_", names(samples_renamed))
names(reads_list) <- gsub("-", "_", names(reads_list))

# Reads lengths distribution
length_dist <- rlength_distr(reads_list, sample = samples_renamed)
for (i in 1:length(samples_renamed)) {
  dir.create(paste0(ribowaltz_folder, samples[i], "/"), showWarnings = FALSE)
  col_selec <- c(1, (i * 2), (i * 2) + 1)
  write.table(length_dist$dt[, ..col_selec],
              paste0(ribowaltz_folder, samples[i], "/reads_distribution_", samples[i], ".csv"),
              quote = FALSE, row.names = FALSE, sep = "\t")
  tiff(file = paste0(ribowaltz_folder, samples[i], "/reads_distribution_", samples[i], ".tiff"))
  plot(length_dist[[paste0("plot_", samples_renamed[i])]])
  dev.off()
}
rm(length_dist)
gc()

# p-site offset calculation
source(function_psite)
psite_offset <- psite_ribowaltz(reads_list,
                                flanking = 6,
                                start = TRUE,
                                extremity = "auto",
                                plot = TRUE,
                                plot_dir = ribowaltz_folder,
                                plot_format = "tiff",
                                cl = 100,
                                txt = TRUE,
                                txt_file = paste0(ribowaltz_folder, "best_offset.csv")
)

reads_psite_list <- riboWaltz::psite_info(reads_list, psite_offset)
write.table(psite_offset, paste0(ribowaltz_folder, "psite_offset.csv"), quote = FALSE, row.names = FALSE, sep = "\t")
rm(list = c("psite_offset", "reads_list"))
gc()

# Proportion of region covered by P-sites
psite_region <- region_psite(reads_psite_list, annotation_db_transcript, sample = samples_renamed)
write.table(psite_region$dt, paste0(ribowaltz_folder, "region_psite.csv"), quote = FALSE, row.names = FALSE, sep = "\t")

tiff(file = paste0(ribowaltz_folder, "region_psite.tiff"))
psite_region[["plot"]]
dev.off()

rm(psite_region)
gc()

# Phasing by read length
frames_stratified <- frame_psite_length(reads_psite_list, sample = samples_renamed, region = "all")
write.table(frames_stratified$dt, paste0(ribowaltz_folder, "frame_psite_length.csv"), quote = FALSE, row.names = FALSE, sep = "\t")

tiff(file = paste0(ribowaltz_folder, "frame_psite_length.tiff"))
frames_stratified[["plot"]]
dev.off()

rm(frames_stratified)
gc()

# Global phasing
frames <- frame_psite(reads_psite_list, sample = samples_renamed, region = "all")
write.table(frames$dt, paste0(ribowaltz_folder, "frame_psite.csv"), quote = FALSE, row.names = FALSE, sep = "\t")

tiff(file = paste0(ribowaltz_folder, "frame_psite.tiff"))
frames[["plot"]]
dev.off()

rm(frames)
gc()
