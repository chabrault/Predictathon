
#' Format and curate vcf file
#'
#' @param vcf.p2f path to the vcf file
#' @param matrix.gt matrix of genotypes, must contains CHROM and POS columns, default is NULL
#' @param p2f.export.vcf path to export curated vcf file, default is NULL (no export), useful for further Beagle imputation.
#' @param IDnum logical, remove the prefix idxx. from the genotype name, default is FALSE
#' @param mrk.info data frame with marker information pulled from the vcf, with columns CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, default is NULL. Needed if `vcf.p2f` is not provided and `p2f.export.vcf` is not null.
#' @param corresp.geno.name table of correspondence between genotype names, with the first column being the name to be matched and the second one being the updated name.
#' @param remove.chrUkn logical, remove markers from unknown chromosome, default is TRUE
#' @param check.mrk.dups logical, check for duplicated markers, default is TRUE
#' @param check.geno.dups logical, check for duplicated genotypes, default is TRUE
#' @param remove.nonPolyMrk logical, remove non-segregating markers, default is TRUE
#' @param tresh.heterozygous numeric, threshold of the maximum proportion of heterozygous genotypes, default is NULL
#' @param thresh.NA.ind numeric, threshold of the maximum proportion missing values for individuals, default is 0.5
#' @param thresh.NA.mrk numeric, threshold of the maximum proportion of missing values for markers, default is 0.5
#' @param format.names logical, use a custom function to transform genotype names, default is FALSE
#' @param concordance.function logical, use a custom function to match genotype names, default is FALSE
#' @param corresp.geno.list named list of correspondence between genotype names, with correct spelling as name and possible synonym as character vector, default is NULL.
#' @param imputation character, impute missing values with kNNI or Beagle or don't impute, default is NULL (no imputation)
#' @param p2f.beagle path to export file for Beagle imputation, default is NULL
#' @param thresh.MAF numeric, threshold of minor allele frequency, default is NULL
#' @param verbose numeric, level of verbosity, default is 1
#'
#' @return data frame with 0/1/2 values from the curated vcf file, with genotypes in row and markers in columns.
#' @author Charlotte Brault
#' @seealso [format_names()], [getGenoTas_to_DF()], [concordance_match_name()]
#' @importFrom vcfR read.vcfR extract.gt getFIX
#' @importFrom methods new
#' @importFrom dplyr arrange filter mutate select distinct
#' @importFrom tibble as_tibble

#' @export

format_curate_vcf <- function(vcf.p2f=NULL,
                              matrix.gt=NULL,
                              mrk.info=NULL,
                              corresp.geno.name=NULL,
                              p2f.export.vcf=NULL,
                              IDnum=FALSE,
                              remove.chrUkn=TRUE,
                              check.mrk.dups=TRUE,
                              remove.nonPolyMrk=TRUE,
                              tresh.heterozygous=NULL,
                              check.geno.dups=TRUE,
                              thresh.NA.ind=0.5,
                              thresh.NA.mrk=0.5,
                              format.names=FALSE,
                              concordance.function=FALSE,
                              corresp.geno.list=NULL,
                              imputation=NULL,
                              p2f.beagle=NULL,
                              thresh.MAF=NULL,
                              verbose=1){

  ### Verifications
  stopifnot(thresh.NA.mrk >= 0, thresh.NA.mrk <= 1,
            thresh.NA.ind >=0, thresh.NA.ind <= 1,
            is.logical(check.geno.dups), is.logical(check.mrk.dups),
            is.logical(remove.chrUkn), is.logical(format.names),
            is.logical(IDnum),is.logical(remove.nonPolyMrk),
            ## check if vcf.p2f or matrix.gt is provided
            xor(is.null(vcf.p2f), is.null(matrix.gt)),
            imputation %in% c(NULL,"kNNI","Beagle"),
            length(imputation) <=1)
  #vcf.file <- as.data.frame(vcf.file)


  ### 1. Load VCF file
  if(!is.null(vcf.p2f)){
    stopifnot(file.exists(vcf.p2f))

    if(verbose > 0) print("Load the vcf file...")
    ## load vcf file and convert to gt
    vcf <- vcfR::read.vcfR(vcf.p2f, verbose = FALSE)
    gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric=FALSE, IDtoRowNames = TRUE)
    meta <- vcf@meta
    mrk.info <- tibble::as_tibble(vcfR::getFIX(vcf))
    mrk.info$POS <- as.numeric(mrk.info$POS)
    mrk.info <- dplyr::arrange(mrk.info, CHROM,POS)
    ## combine marker information with genotype
    vcf.file <- cbind(mrk.info[,c("CHROM","POS")],gt[mrk.info$ID,])
    vcf.file$POS <- as.numeric(vcf.file$POS)

    rm(vcf, gt)

  } else {
    stopifnot(!is.null(matrix.gt), all(c("CHROM","POS") %in% colnames(matrix.gt)[c(1,2)]))
    vcf.file <- matrix.gt
    meta <- NULL
  }

  if(verbose >0){
    print(paste0("Initial dimensions: ", nrow(vcf.file)," markers and ",
                 (ncol(vcf.file) - 2), " genotypes"))
  }

  ## remove the prefix idxx. from the genotype name
  if(IDnum){
    colnames(vcf.file) <- gsub(x=colnames(vcf.file),pattern="^id[0-9]+\\.",replacement="")
  }
  ## replace patterns "./." by NA
  vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],2, function(x)
    gsub(pattern="\\.",x=x,replacement=NA)))
  ### replace multi-allelic marker by mono-allelic, i.e., "2/2"" by "1/1"
  # vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],2, function(x)
  #   gsub(pattern="2/2",x=x,replacement="1/1")))
  vcf.file[,3:ncol(vcf.file)] <- as.data.frame(apply(vcf.file[,3:ncol(vcf.file)],
                                                     2, function(x){
                                                       x=gsub(pattern="2/2",x=x,replacement="1/1")
                                                       x=gsub(pattern="0/2",x=x,replacement="0/1")
                                                       x=gsub(pattern="2/0",x=x,replacement="1/0")
                                                     }))
  ## set position as numeric
  vcf.file$POS <- as.numeric(stringr::str_trim(vcf.file$POS))

  ## update genotype names
  if(!is.null(corresp.geno.name)){
    stopifnot(ncol(corresp.geno.name) >=2)
    if(verbose >0){
      print("Update genotype names..")
    }

    colnames(vcf.file) <- plyr::mapvalues(colnames(vcf.file),
                                          from=corresp.geno.name[,1],
                                          to=corresp.geno.name[,2], warn_missing=FALSE)
  }

  ## handle unknown chromosome, remove markers if few
  if(remove.chrUkn){
    chr.unk <- c("CHRUnknown","U","ChrUnknown","Unknown","Un","ChrUn","chrU","ChrU")
    if(verbose >0){
      print(paste0("Removed ",nrow(vcf.file[vcf.file$CHROM %in% chr.unk,]),
                   " markers from unknown chromosome"))
    }
    vcf.file <- vcf.file[!vcf.file$CHROM %in% chr.unk,]

  }

  ## check for duplicated markers
  # extract map information
  vcf.file$chr_pos <- paste0(vcf.file$CHROM, "_",vcf.file$POS)
  if(check.mrk.dups){
    dup_snp <- unique(vcf.file$chr_pos[duplicated(vcf.file$chr_pos)])
    if(verbose >0){
      print(paste0("Found ",length(dup_snp)," duplicated markers"))
    }
    ## for duplicated loci, keep the one with the lowest missing value
    idx2rem <- c()
    for(mrk in dup_snp){
      ## get row id for duplicated markers
      idx.mrkdup <- grep(mrk,vcf.file$chr_pos)
      ## detect rows with less missing values
      idx2rem <- c(idx2rem,idx.mrkdup[-which.min(rowSums(is.na(vcf.file[vcf.file$chr_pos == mrk,])))])
    }
    ## remove duplicated columns with more missing values
    if(length(idx2rem)>0){
      vcf.file <- vcf.file[-idx2rem,]
    }
  }
  rownames(vcf.file) <- vcf.file$chr_pos
  ## if there is no chr before chr_pos, add it
  if(all(grepl("[1-7]|U",substr(rownames(vcf.file),1,1)))){
    rownames(vcf.file) <- paste0("chr",rownames(vcf.file))
  }
  cols2rem <- c("CHROM","POS", "chr_pos")
  vcf.file[,cols2rem]<- NULL

  ## 2. remove duplicated genotypes

  ## format genotype name with custom function format_names
  if(format.names){
    colnames(vcf.file) <- format_names(colnames(vcf.file))$corrected
  } else if(any(grep(" ", colnames(vcf.file))) & !is.null(p2f.export.vcf)){
    ## need to remove white spaces for synbreed export, replace it by underscore
    colnames(vcf.file) <- gsub(x=colnames(vcf.file),pattern=" ",replacement = "_")
  }
  if(concordance.function){
    stopifnot(!is.null(corresp.geno.list))
    conc <- concordance_match_name(match.name=colnames(vcf.file), corresp.list=corresp.geno.list)
    colnames(vcf.file) <- conc$new.names
  }

  if(check.geno.dups){
    ## classic detection of duplicated genotypes
    dup.genos1 <- unique(colnames(vcf.file)[duplicated(colnames(vcf.file))])

    ### detect duplicated genotypes with .1/.2/... appended at the end of the name
    dup.genos2 <- grep(pattern="\\.[1-9]{1,2}$", colnames(vcf.file), value=TRUE)
    dup.genos2 <- unique(stringr::str_remove(dup.genos2, pattern="\\.[1-9]{1,2}$"))
    dup.genos <- unique(c(dup.genos1,dup.genos2))
    if(verbose >0){
      print(paste0("Found ",length(dup.genos)," duplicated genotypes"))
    }
    ## for duplicated genotypes, fill if missing data and keep the one with less missing data
    for(gen in dup.genos){
      idx.genodup <- grep(paste0("^",gen,"[.1-9]{0,2}"),colnames(vcf.file))
      #colnames(vcf.file.sel)[idx.genodup]
      if(length(idx.genodup) > 1){
        tmp <- vcf.file[,idx.genodup]

        ## find the column with the least NA values
        idxLessNa <- which.min(colSums(is.na(tmp)))
        ## fill NA values with data from other columns if available
        idxna <- which(is.na(tmp[,idxLessNa]))
        for(i in 1:ncol(tmp)){
          tmp[idxna,idxLessNa] <- dplyr::coalesce(tmp[idxna,idxLessNa],tmp[idxna,i])
        }
        ## suppress duplicated columns, keep the one with less NA values and filled
        vcf.file[,idx.genodup[idxLessNa]] <- tmp[,idxLessNa]
        colnames(vcf.file)[idx.genodup[idxLessNa]] <- gen
        idx2rem <- setdiff(idx.genodup,idx.genodup[idxLessNa])
        vcf.file = vcf.file[,-idx2rem]
      }
    }
    stopifnot(length(unique(colnames(vcf.file))) == ncol(vcf.file))
  }


  ## 3. Curate: remove SNPs/ genos with too many missing data

  ### remove markers with too many NA
  snp.miss <- apply(vcf.file, 1, function(x) sum(is.na(x)))
  snp2rem <- names(snp.miss)[snp.miss > thresh.NA.mrk*ncol(vcf.file)]
  if(verbose > 0){
    print(paste0("Removed ",length(snp2rem)," markers with too many missing data"))
  }
  vcf.file <- vcf.file[!rownames(vcf.file) %in% snp2rem,]
  if(nrow(vcf.file) == 0){
    stop("No markers left after removing markers with too many missing data")
  }

  ### remove individuals with too many NA
  gen.miss <- apply(vcf.file, 2, function(x) sum(is.na(x)))
  geno2rem <- names(gen.miss)[gen.miss > thresh.NA.ind*nrow(vcf.file)]
  if(verbose > 0){
    print(paste0("Removed ",length(geno2rem)," individuals with too many missing data"))
  }
  vcf.file <- vcf.file[,!colnames(vcf.file) %in% geno2rem]
  if(ncol(vcf.file) < 4){
    stop("No genotypes left after removing genotypes with too many missing data")
  }



  ## Remove non-polymorphic markers
  if(remove.nonPolyMrk){
    nb.allele <- apply(vcf.file, 1, function(x) length(table(t(x))))
    mrk2rem <- names(nb.allele)[nb.allele <= 1]
    if(verbose > 0){
      print(paste0("Removed ",length(mrk2rem)," non-polymorphic markers"))
    }
    if(length(mrk2rem) >0) vcf.file <- vcf.file[!rownames(vcf.file) %in% mrk2rem,]
    rm(mrk2rem)
  }

  ## Remove markers with too many heterozygous genotypes
  if(!is.null(tresh.heterozygous)){
    hetero <- apply(vcf.file, 2, function(x) {sum(x %in% c("0/1","1/0","0|1","1|0"), na.rm = TRUE)})
    geno2rem <- names(hetero)[hetero > tresh.heterozygous * nrow(vcf.file)]
    if(verbose > 0){
      print(paste0("Removed ",length(geno2rem)," genotypes with too many heterozygous markers"))
    }
    vcf.file <- vcf.file[,!colnames(vcf.file) %in% geno2rem,]
  }

  rownames(vcf.file) <- gsub(x = rownames(vcf.file),pattern="Chr",replacement = "chr")
  ## recode numerically: replace 0/0=0, 0/1 or 1/0=1 and 1/1=2
  rowN <- rownames(vcf.file)

  vcf.num <- apply(vcf.file, 2, function(x) {
    as.numeric(stringr::str_replace_all(gsub("[[:punct:]]", "", x),
                                        c("00"="0","01"="1","10"="1","11"="2")))})
  vcf.num <- as.data.frame(vcf.num)
  rownames(vcf.num) <- rowN

  ### 4. (optional) export to a formatted vcf file for Beagle imputation

  if(!is.null(p2f.export.vcf)){
    ## if no extension or vcf extension provided, extend to vcf.gz extension
    if(!grepl("vcf.gz$",p2f.export.vcf)){
      p2f.export.vcf <- paste0(gsub(".vcf","",p2f.export.vcf),".vcf.gz")
    }

    ## recreate vcfR object
    vcf2exp <- methods::new(Class = "vcfR")
    ### meta element
    if(!is.null(meta)){
      vcf2exp@meta <- meta
    } else {
      vcf2exp@meta <- c("##fileformat=VCFv4.2",
                        "##reference=xxx",
                        "##comment=Generated from format_curate_vcf function",
                        #"##filedate= 20065",
                        "##source='write.vcf of vcfR'",
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>")
    }

    ### fix element: add marker information
    #### define a default mrk.info data frame if not provided
    if(is.null(vcf.p2f) & is.null(mrk.info)){
      chr_pos <- rownames(vcf.num)
      mrk.info <- data.frame(CHROM=unlist(lapply(strsplit(chr_pos,"_"),function(x) x[1])),
                             POS=unlist(lapply(strsplit(chr_pos,"_"),function(x) x[2])),
                             ID=chr_pos,
                             REF="A",ALT="T",
                             QUAL=".",FILTER="PASS",INFO=".",
                             FORMAT="GT")
      mrk.info$CHROM <- gsub("Chr","chr",mrk.info$CHROM)

    } else {
      ## use the mrk.info from the vcf file, subset selected markers
      mrk.info$ID <- paste0(gsub("Chr","chr",mrk.info$CHROM),"_",mrk.info$POS)
      mrk.info <- mrk.info[match(rownames(vcf.num), mrk.info$ID),]
    }
    if(!all(c("INFO","FORMAT") %in% colnames(mrk.info))){
      mrk.info <- cbind(mrk.info,
                        #QUAL=".",FILTER="PASS",
                        INFO=".",FORMAT="GT")
      mrk.info <- mrk.info[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")]
    }
    mrk.info <- mrk.info[match(rowN, mrk.info$ID),]
    stopifnot(nrow(mrk.info) == nrow(vcf.num))
    mrk.info <- as.matrix(mrk.info)
    mrk.info[,"POS"] <- gsub(" ","",mrk.info[,"POS"])
    vcf2exp@fix <- mrk.info

    ### gt element, coded in 0/0, 0/1, 1/1
    vcf2exp@gt <- as.matrix(vcf.file)

    ### export
    vcfR::write.vcf(x=vcf2exp, file=p2f.export.vcf, mask=FALSE, APPEND=FALSE)

  } ## end if !is.null(p2f.export.vcf)


  ### Imputation of missing markers with kNNI
  if(length(imputation) >0 && imputation == "kNNI"){
    if(is.null(p2f.export.vcf)){
      stop("Need to export vcf file to perform imputation")
    }
    ## load rTASSEL and rJava packages
    if(!requireNamespace("rTASSEL", quietly = TRUE)){
      stop("Please install rTASSEL package to perform kNNI imputation")
    }
    if(!requireNamespace("rJava", quietly = TRUE)){
      stop("Please install rJava package to perform kNNI imputation")
    }

    tasGenoVCF <- rTASSEL::readGenotypeTableFromPath(path = p2f.export.vcf)
    ## use KNNI as imputation method, using default parameters
    tasGenoImp <- rTASSEL::imputeLDKNNi(tasObj=tasGenoVCF, highLDSSites = 30,
                                        knnTaxa = 10, maxDistance = 1e+07)
    ## use an internal function to convert this object to a data frame
    vcf.file <- getGenoTas_to_DF(tasGenoImp)

    ## format it to export it as a vcf file
    vcf2exp@gt <- as.matrix(vcf.file)
    # gp$geno <- t(vcf.file[rownames(gp$map),rownames(gp$geno)])
    # stopifnot(identical(colnames(gp$geno),rownames(gp$map)))

    ## detect and warn if missing values left

    if(verbose > 0) {
      ### proportion of missing value per column
      print("Summary of missing data percentage per marker: ")
      print(summary(colMeans(is.na(gp$geno)) *100))

      ### proportion of missing value per line
      print("Summary of missing data percentage per genotype: ")
      print(summary(rowMeans(is.na(gp$geno)) *100))

    }
    # # export VCF
    # synbreed::write.vcf(gp,file=p2f.export.vcf)
    ### export
    vcfR::write.vcf(x=vcf2exp, file=p2f.export.vcf, mask=FALSE, APPEND=FALSE)

  } else if(length(imputation) >0 && imputation == "Beagle"){
    if(is.null(p2f.export.vcf)){
      stop("Need to export vcf file to perform imputation")
    }
    if(!file.exists(p2f.beagle)){
      stop("Beagle jar file not found")
    }
    ## run Beagle imputation
    p2f.imp <- gsub(".vcf.gz","_imputed",p2f.export.vcf)
    system(paste0("java -Xmx2g -jar ",p2f.beagle," gt=",p2f.export.vcf,
                  " out=",p2f.imp,
                  " nthreads=4 ne=100000"))

    ## load imputed vcf
    vcf.imp <- vcfR::read.vcfR(paste0(p2f.imp,".vcf.gz"), verbose = FALSE)
    gt.imp <- vcfR::extract.gt(vcf.imp, element = "GT", as.numeric=FALSE, IDtoRowNames = FALSE)
    mrk.info.imp <- tibble::as_tibble(vcfR::getFIX(vcf.imp))
    vcf.file <- cbind(mrk.info.imp[,c("CHROM","POS")],gt.imp)
    vcf.file$POS <- as.numeric(vcf.file$POS)
    vcf.file <- dplyr::arrange(vcf.file, CHROM,POS)
    rownames(vcf.file) <- paste0(vcf.file$CHROM, "_",vcf.file$POS)
    cols2rem <- c("CHROM","POS")
    vcf.file[,cols2rem]<- NULL
    rowN <- rownames(vcf.file)
    ## recode numerically: replace 0/0=0, 0/1 or 1/0=1 and 1/1=2
    vcf.num <- apply(vcf.file, 2, function(x) {
      as.numeric(stringr::str_replace_all(gsub("[[:punct:]]", "", x),
                                          c("00"="0","01"="1","10"="1","11"="2")))})
    vcf.num <- as.data.frame(vcf.num)
    rownames(vcf.num) <- rowN
  } ## end if Beagle imput


  ## Remove markers with MAF below threshold after imputation
  if(!is.null(thresh.MAF)){
    thresh.MAF <- as.numeric(thresh.MAF)
    if(thresh.MAF < 0 | thresh.MAF > 0.5){
      stop("MAF should be between 0 and 0.5")
    }

    maf <- apply(vcf.num, 2, function(x) {
      t = table(x) ; t = t[names(t) %in% c(0,2)]
      if(length(t) < 2) maf = 0
      else maf = sort(t)[1]/length(x)
      return(maf)
    })

    mrk2rem <- names(maf)[maf < thresh.MAF]

    if(verbose > 0){
      print(paste0("Removed ",length(mrk2rem)," markers with MAF below ",thresh.MAF))
    }
    vcf.num <- vcf.num[!rownames(vcf.num) %in% mrk2rem,]
  }


  if(verbose >0){
    print(paste0("Final dimensions: ",ncol(vcf.num), " genotypes and ", nrow(vcf.num)," markers"))
  }
  return(t(vcf.num)) ## output in dimension geno x marker
}
