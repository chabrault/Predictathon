here::here()


library(BrAPI)
library(data.table)

suppressPackageStartupMessages(library(here)) # root directory
suppressPackageStartupMessages(library(data.table)) # load data
suppressPackageStartupMessages(library(apercu)) # ap
suppressPackageStartupMessages(library(dplyr)) # formatting
suppressPackageStartupMessages(library(purrr)) # handle lists
## devtools::install_github("chabrault/HRSWGS") # custom package instal
#suppressPackageStartupMessages(library(HRSWGS)) # custom package
## devtools::install_github("TriticeaeToolbox/BrAPI.R")
suppressPackageStartupMessages(library(BrAPI)) # T3 connection


here::i_am("analysis/pullGenoTrial.Rmd")





# Create T3 connection

wheat <- createBrAPIConnection("wheat.triticeaetoolbox.org", is_breedbase = TRUE)


#One trial
trial_name <- "2025_AYT_Aurora"

acc <- wheat$wizard("accessions", list(trials = trial_name), verbose = FALSE)

acc_ids   <- acc$data$ids
acc_names <- acc$data$names

length(acc_ids)
head(acc_names)

##Automatically find which protocols actually have a downloadable VCF for Aurora

prot <- wheat$wizard("genotyping_protocols", list(accessions = acc_ids), verbose = FALSE)
prot_ids   <- prot$data$ids
prot_names <- prot$data$names

proj_ids_chr <- as.character(proj_ids)

vcf_hits_list <- lapply(seq_along(prot_ids), function(i) {
  pid   <- prot_ids[i]
  pname <- prot_names[i]
  
  vcf_list <- wheat$vcf_archived_list(
    genotyping_protocol_id = pid,
    verbose = FALSE
  )
  
  if (is.null(vcf_list) || nrow(vcf_list) == 0) return(NULL)
  
  hits <- vcf_list[vcf_list$project_id %in% proj_ids_chr, , drop = FALSE]
  if (nrow(hits) == 0) return(NULL)
  
  hits$protocol_name <- pname
  hits
})

vcf_hits_all <- rbindlist(vcf_hits_list, fill = TRUE)

vcf_hits_all

#Download the VCF

trial_name <- "2025_AYT_Aurora"

chosen <- vcf_hits_all[protocol_id == "301" & project_id == "11050", ]
chosen


#Output folder 

outdir <- here("data", "genotypic")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

out_vcf <- file.path(outdir, paste0(trial_name, "_GBS_SDSU_2025.vcf.gz"))
out_vcf


wheat$vcf_archived(
  output = out_vcf,
  genotyping_protocol_id = as.integer(chosen$protocol_id),
  genotyping_project_id  = as.integer(chosen$project_id),
  file_name              = chosen$file_name,
  verbose = TRUE
)

file.exists(out_vcf)
file.info(out_vcf)$size


####Confirm wheather the VCF file has Aurora lines
#Sample names to a txt file

samples_file <- file.path(outdir, "Aurora_samples.txt")
writeLines(acc_names, samples_file)
samples_file



out_vcf
file.exists(out_vcf)
file.info(out_vcf)


library(vcfR)
library(rrBLUP)

#Read the QC VCF and extract genotype matrix

vcf_path <- "/Users/marcos/Documents/GitHub/Predictathon/data/genotypic/2025_AYT_Aurora_GBS_SDSU_2025_AuroraOnly.QC.vcf.gz"

vcf <- read.vcfR(vcf_path, verbose = FALSE)

# Extract genotype calls
gt <- extract.gt(vcf, element = "GT")   # matrix markers x samples, values like "0/0", "0/1", "1/1", "./."

# Convert to numeric dosage: 0,1,2
dosage <- matrix(NA_real_, nrow = nrow(gt), ncol = ncol(gt),
                 dimnames = dimnames(gt))

dosage[gt %in% c("0/0","0|0")] <- 0
dosage[gt %in% c("0/1","1/0","0|1","1|0")] <- 1
dosage[gt %in% c("1/1","1|1")] <- 2

# transpose to individuals x markers (rrBLUP convention often uses individuals in rows)
M <- t(dosage)

dim(M)        # should be 125 x 11738
M[1:5,1:5]

#Simple imputation and build a genomic relationship matrix
# Mean-impute missing genotypes (fast baseline)
for (j in 1:ncol(M)) {
  mj <- M[, j]
  M[is.na(mj), j] <- mean(mj, na.rm = TRUE)
}

G <- A.mat(M)   # VanRaden G matrix
dim(G)

