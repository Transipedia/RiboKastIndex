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
  make_option(c("-w", "--work_dir"), type = "character", help = "Path to working directory"),
  make_option(c("-f", "--feature_type"), type = "character", help = "Path to designed R functions"),
  make_option(c("-g", "--gtf"), type = "character", help = "Path and name of the input GTF file"),
  make_option(c("--path_prefix"), type = "character", help = "Prefix for RESULTS path"),
  make_option(c("-b", "--bam"), type = "character", help = "Path to single BAM file"),
  make_option(c("-n", "--sample_name"), type = "character", help = "Sample name for output labeling")
)

opt = parse_args(OptionParser(option_list = option_list))

# Parameters
local_path <- opt$work_dir
results_all_path <- opt$path_prefix
function_psite <- paste0(opt$feature_type, "ribowaltz_psite_with_NA_control.R")
gtf_file <- opt$gtf
bam_file <- opt$bam
sample_name <- opt$sample_name

# Read config.yaml
params <- scan(file = paste0(local_path, "config.yaml"), what = "character", sep = ":")
readsLength_min <- gsub(" ", "", params[which(params == "readsLength_min") + 1], fixed = TRUE)
readsLength_max <- gsub(" ", "", params[which(params == "readsLength_max") + 1], fixed = TRUE)

# Output directory per sample
ribowaltz_folder <- file.path(results_all_path,
                              paste0("riboWaltz.", readsLength_min, "-", readsLength_max),
                              sample_name)
dir.create(ribowaltz_folder, recursive = TRUE, showWarnings = FALSE)

#########################
####### EXECUTION #######
#########################

# Load annotation
annotation_db <- riboWaltz::create_annotation(gtf_file)
annotation_db_transcript <- data.table(annotation_db)[l_cds > 0, ]

# Handle single BAM file
bam_basename <- basename(bam_file)
samples <- c(sample_name)
names(samples) <- str_remove(bam_basename, "\\.bam")

# Temporary directory for bam
tmp_bam_dir <- tempfile(pattern = "bam_input_")
dir.create(tmp_bam_dir)
file.copy(from = bam_file, to = file.path(tmp_bam_dir, bam_basename))

# Read BAM
reads_list <- riboWaltz::bamtolist(
  bamfolder = tmp_bam_dir,
  annotation = annotation_db_transcript,
  name_samples = samples
)

# Clean sample names
samples_renamed <- gsub("-", "_", samples)
names(samples_renamed) <- gsub("-", "_", names(samples_renamed))
names(reads_list) <- gsub("-", "_", names(reads_list))

# Reads length distribution
length_dist <- rlength_distr(reads_list, sample = samples_renamed)
col_selec <- c(1, 2, 3)
write.table(length_dist$dt[, ..col_selec],
            file.path(ribowaltz_folder, paste0("reads_distribution_", sample_name, ".csv")),
            quote = FALSE, row.names = FALSE, sep = "\t")
tiff(file = file.path(ribowaltz_folder, paste0("reads_distribution_", sample_name, ".tiff")))
plot(length_dist[[paste0("plot_", samples_renamed)]])
dev.off()

# P-site offset
source(function_psite)
psite_offset <- psite_ribowaltz(
  reads_list,
  flanking = 6,
  start = TRUE,
  extremity = "auto",
  plot = TRUE,
  plot_dir = ribowaltz_folder,
  plot_format = "tiff",
  cl = 100,
  txt = TRUE,
  txt_file = file.path(ribowaltz_folder, "best_offset.csv")
)

# Add p-site info
reads_psite_list <- riboWaltz::psite_info(reads_list, psite_offset)
write.table(psite_offset,
            file.path(ribowaltz_folder, "psite_offset.csv"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# Region p-site
psite_region <- region_psite(reads_psite_list, annotation_db_transcript, sample = samples_renamed)
write.table(psite_region$dt,
            file.path(ribowaltz_folder, "region_psite.csv"),
            quote = FALSE, row.names = FALSE, sep = "\t")
tiff(file = file.path(ribowaltz_folder, "region_psite.tiff"))
psite_region[["plot"]]
dev.off()

# Frame p-site length
frames_stratified <- frame_psite_length(reads_psite_list, sample = samples_renamed, region = "all")
write.table(frames_stratified$dt,
            file.path(ribowaltz_folder, "frame_psite_length.csv"),
            quote = FALSE, row.names = FALSE, sep = "\t")
tiff(file = file.path(ribowaltz_folder, "frame_psite_length.tiff"))
frames_stratified[["plot"]]
dev.off()

# Global frame p-site
frames <- frame_psite(reads_psite_list, sample = samples_renamed, region = "all")
write.table(frames$dt,
            file.path(ribowaltz_folder, "frame_psite.csv"),
            quote = FALSE, row.names = FALSE, sep = "\t")
tiff(file = file.path(ribowaltz_folder, "frame_psite.tiff"))
frames[["plot"]]
dev.off()

