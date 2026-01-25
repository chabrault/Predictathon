# ---- Setup ----
suppressPackageStartupMessages({
  library(here)
  library(BrAPI)
  library(data.table)
  library(vcfR)
  library(rrBLUP)
})

here::i_am("code/Template_marcos.r")


here::here()

wheat <- createBrAPIConnection("wheat.triticeaetoolbox.org", is_breedbase = TRUE)

# ---- USER INPUTS (change only these) ----
trial_name <- "AWY1_DVPWA_2024"

# ---- Output folders ----
outdir <- here("data", "genotypic", trial_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Accessions in trial ----
acc <- wheat$wizard("accessions", list(trials = trial_name), verbose = FALSE)
acc_ids   <- acc$data$ids
acc_names <- acc$data$names

length(acc_ids)
head(acc_names)


samples_file <- file.path(outdir, paste0(trial_name, "_samples.txt"))
writeLines(acc_names, samples_file)
samples_file

# ---- Find genotyping projects + protocols for these accessions ----
prot <- wheat$wizard("genotyping_protocols", list(accessions = acc_ids), verbose = FALSE)
prot_ids   <- prot$data$ids
prot_names <- prot$data$names

proj <- wheat$wizard("genotyping_projects", list(accessions = acc_ids), verbose = FALSE)
proj_ids_chr <- as.character(proj$data$ids)

length(prot_ids)
length(proj_ids_chr)

data.frame(protocol_id = prot_ids, protocol_name = prot_names)



# ---- Find archived VCFs that match BOTH protocol and project ----
vcf_hits_list <- lapply(seq_along(prot_ids), function(i) {
  pid   <- prot_ids[i]
  pname <- prot_names[i]
  
  vcf_list <- wheat$vcf_archived_list(genotyping_protocol_id = pid, verbose = FALSE)
  if (is.null(vcf_list) || nrow(vcf_list) == 0) return(NULL)
  
  hits <- vcf_list[vcf_list$project_id %in% proj_ids_chr, , drop = FALSE]
  if (nrow(hits) == 0) return(NULL)
  
  hits$protocol_name <- pname
  hits
})

vcf_hits_all <- rbindlist(vcf_hits_list, fill = TRUE)
vcf_hits_all

# ---- Pick a VCF to download (must exist) ----
chosen <- vcf_hits_all[file_name == "2025-12-10_17:57:39_fileRm9S", ]
stopifnot(nrow(chosen) == 1)

label <- "ThermoFisher_AgriSeq_5K"

raw_vcf_gz <- file.path(outdir, paste0(trial_name, "_", label, ".vcf.gz"))

wheat$vcf_archived(
  output = raw_vcf_gz,
  genotyping_protocol_id = as.integer(chosen$protocol_id),
  genotyping_project_id  = as.integer(chosen$project_id),
  file_name              = chosen$file_name,
  verbose = TRUE
)

raw_vcf_gz
file.exists(raw_vcf_gz)
file.info(raw_vcf_gz)$size
samples_file




acc <- wheat$wizard("accessions", list(trials = "AWY1_DVPWA_2024"), verbose = FALSE)
acc_ids <- acc$data$ids



prot <- wheat$wizard("genotyping_protocols",
                     list(accessions = acc_ids),
                     verbose = FALSE)






######Terminal part


TRIAL="AWY1_DVPWA_2024"
LABEL="ThermoFisher_AgriSeq_5K"
BASE="/Users/marcos/Documents/GitHub/Predictathon/data/genotypic/${TRIAL}"

RAW="${BASE}/${TRIAL}_${LABEL}.vcf.gz"
VCF="${BASE}/${TRIAL}_${LABEL}.vcf"
BGZ="${BASE}/${TRIAL}_${LABEL}.bgz.vcf.gz"
SAMPLES="${BASE}/${TRIAL}_samples.txt"

TRIAL_ONLY="${BASE}/${TRIAL}_${LABEL}.TrialOnly.vcf.gz"
QC="${BASE}/${TRIAL}_${LABEL}.TrialOnly.QC.vcf.gz"

file "$RAW"
mv "$RAW" "$VCF"
bgzip -c "$VCF" > "$BGZ"
bcftools index -c "$BGZ"

bcftools view --force-samples -S "$SAMPLES" -Oz -o "$TRIAL_ONLY" "$BGZ"
bcftools index -c "$TRIAL_ONLY"

bcftools +fill-tags "$TRIAL_ONLY" -- -t AF,MAF,F_MISSING \
| bcftools view -i 'MAF>=0.05 && F_MISSING<=0.20' -Oz -o "$QC"
bcftools index -c "$QC"


wc -l "/Users/marcos/Documents/GitHub/Predictathon/data/genotypic/AWY1_DVPWA_2024/AWY1_DVPWA_2024_samples.txt"

bcftools query -l "/Users/marcos/Documents/GitHub/Predictathon/data/genotypic/AWY1_DVPWA_2024/AWY1_DVPWA_2024_ThermoFisher_AgriSeq_5K.TrialOnly.vcf.gz" | wc -l

bcftools view -H "/Users/marcos/Documents/GitHub/Predictathon/data/genotypic/AWY1_DVPWA_2024/AWY1_DVPWA_2024_ThermoFisher_AgriSeq_5K.TrialOnly.vcf.gz" | wc -l

bcftools view -H "/Users/marcos/Documents/GitHub/Predictathon/data/genotypic/AWY1_DVPWA_2024/AWY1_DVPWA_2024_ThermoFisher_AgriSeq_5K.TrialOnly.QC.vcf.gz" | wc -l


### Back to R code to convert QC VCF → M and save