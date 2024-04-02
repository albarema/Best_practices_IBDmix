## Rscript calc_precision_200Mb.R -m {input.model} -l 100000,50000 --nean nea1_2,nea2_1 --popn 1,2 --age {params.age} --outf {output}

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(GenomicRanges))
quiet(library(tidyverse))
quiet(library(glue))
quiet(library("optparse"))


option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, help="infile name"),
  make_option(c("-m", "--model"), type="character", default=NULL, help="model name"),
  make_option(c("-d", "--lod"), type="character", default=3, help="LOD n"),
  make_option(c("-a", "--minlod"), type="character", default=3, help="minimum LOD to consider"),
  make_option(c("-l", "--segl"), type="character", default=50000, help="Segment length minimum"),
  make_option(c("-b", "--minsegl"), type="character", default=50000, help="Miimum Segment length to be considered"),
  make_option(c("-n", "--nean"), type="character", default=NULL, help="neanderthal ID"),
  make_option(c("-p", "--popn"), type="character", default=NULL, help="pop ID"),
  make_option(c("-c", "--chrn"), type="character", default=1, help="chr n"),
  make_option(c("-t", "--tracts"), type="character", default=NULL, help="tracts tsv"),
  make_option(c("-o", "--outf"), type="character", default=NULL, help="outfile name"),
  make_option(c("-e", "--acc"), type="character", default=NULL, help="missestimation name")
)
  
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tracts <- opt$tracts
infile <- opt$infile
mod <- opt$model
LOD <- opt$lod
chr <- opt$chrn
outfile <- opt$outf
segl <- as.character(opt$segl)
nea <- opt$nean
minlod <- opt$minlod
minsegl <- as.character(opt$minsegl)
pop2test <- opt$popn
missout <- opt$acc

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#mod = "modelA"
#LOD=3 # seq(2,5)
#segl= c('100000', '50000')
#chr <- 1 # chromosome 
#nea <- c('nea2_1', 'nea1_2') 
#popn = c(1,2)
# pop1_1 - pop1_86 pop1.ancient
#pop1_87 - pop1_186

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
read_tracks <- function(mod, tracts){
  # Read true segments
  ttracks <-read_tsv(tracts,show_col_types = FALSE) %>%
    dplyr:::rename(start = left, end = right)
  # merge individual ID
  true_segments <- ttracks %>%
    mutate(chrom = "chr1") %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  return(true_segments)
}

read_ibd <- function(infile, nea, chr, minsegl, minlod, pop2test){
  #ibd_estimated = read_tsv(paste(mod, "/ibdmix_temp_with_sites/ibd_summary_combined/", nea, "_", chr ,"_", LOD, "_", segl,".txt", sep=""),show_col_types = FALSE)
  ibd_estimated = read_tsv(infile, show_col_types = FALSE) %>% filter(pop == pop2test)  %>% distinct()
  maxlod = max(as.numeric(ibd_estimated$slod))
  maxseg = max(as.numeric(ibd_estimated$length))
  if(maxlod >= minlod & maxseg >= minsegl){
    ibd_estimated = ibd_estimated %>% filter(slod >= as.numeric(minlod), length  >= as.numeric(minsegl)) %>%
    mutate(chrom = "chr1") 
    if (dim(ibd_estimated)[1] == 0){ibd_estimated <- c()}
    if (dim(ibd_estimated)[1] >= 1){
      ibd_estimated=ibd_estimated %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      ibd_estimated$ances <- sapply(strsplit(ibd_estimated$ID, "_"),function(x) x[[1]])
    }
  }
  if(maxlod < minlod | maxseg < minsegl){ibd_estimated <- c()}
  return(ibd_estimated)
}


get_overlap <- function(query, subject,hits){
  overlaps<- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
  overlaps$proportion_subject <- width(overlaps) / width(subject[subjectHits(hits)])
  overlaps$proportion_query <- width(overlaps) / width(query[queryHits(hits)])
  return(overlaps)
}

get_confusionM_segcounts <- function(query_id, subject_id, hits,overlaps){
  # Get counts
  FP = length(query_id[which(!query_id %in% query_id[queryHits(hits)]),])
  FN = length(subject_id[which(!subject_id %in% subject_id[subjectHits(hits)]),])
  FN_info = subject_id[which(!subject_id %in% subject_id[subjectHits(hits)]),]
  TP=length(overlaps)
  recall_seg = TP/(TP+FN)
  precision_seg =  TP/(TP+FP)
  F1_seg=2*(precision_seg*recall_seg/(precision_seg+recall_seg))
  finalvec <- c(recall_seg, precision_seg, F1_seg) #nMCC_seg
  return(finalvec)
}

get_confucionM_window_subject_based <- function(overlaps, query_id, subject_id, hits){
  lwin=sum(width(subject_id))
  TP=sum(width(overlaps))/lwin
  FP=(sum(width(query_id)) - sum(width(overlaps)))/lwin
  FN = (sum(width(subject_id)) - sum(width(overlaps)))/lwin
  TN= ((4e8 - sum(width(subject_id))) - (sum(width(query_id)) - sum(width(overlaps))))/lwin # full minus the introgressed minus the FP 
  recall_wi = TP/(TP+FN)
  precision_wi =  TP/(TP+FP)
  F1_wi=2*(precision_wi*recall_wi/(precision_wi+recall_wi))
  #F1 = (2*TP)/(2*TP + FP+FN)
  MCC_wi=(as.numeric(TP)*as.numeric(TN)-as.numeric(FP)*as.numeric(FN))/sqrt(as.numeric(TP+FP)*as.numeric(TP+FN)*as.numeric(TN+FP)*as.numeric(TN+FN))
  nMCC_wi = (MCC_wi+1) /2
  FPR_wi=FP/(FP+TN)
  TPR_wi=TP/(TP+FN)
  acc_wi = (as.numeric(TP) + as.numeric(TN))/(as.numeric(TP) + as.numeric(FN) + as.numeric(TN) + as.numeric(FP))
  bal_acc_wi =0.5* ((as.numeric(TP)/(as.numeric(TP)+as.numeric(FN))) + (as.numeric(TN)/(as.numeric(TN)*as.numeric(FP))))
  finalvec <- c(recall_wi, precision_wi, F1_wi, nMCC_wi,FPR_wi,TPR_wi, acc_wi, bal_acc_wi, FP, FN, TN, TP)
  return(finalvec)
}

get_missestimation <- function(overlaps, query_id, subject_id, hits){
    overlap_introgression=sum(width(overlaps))/4e8
    detected_introgression = sum(width(query_id))/4e8
    true_introgression=sum(width(subject_id))/4e8
    # accuracy=paste(overlaps$proportion_subject, collapse=",")
    min_acc=min(overlaps$proportion_subject)
    max_acc=max(overlaps$proportion_subject)
    prop_true=mean(overlaps$proportion_subject)
    prop_detected=mean(overlaps$proportion_query)
    totvec <- c(mod, as.character(truelod), as.character(truesegl), as.character(minlod), as.character(minsegl), nea, pop2test, name[m], prop_detected, prop_true, detected_introgression,overlap_introgression,true_introgression)
          
    return(totvec)
}

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
true_segments <- read_tracks(mod, tracts)
ftab <- c()
misstab = c()

ibd_popX = read_ibd(infile, nea, chr, as.numeric(minsegl),as.numeric(minlod), pop2test)


if (length(ibd_popX) == 0){ 
  tbl_colnames= c('mod', 'truelod', 'truesegl', 'minlod', 'minsegl', 'nea', 'popX', 'name', 'recall_seg','precision_seg', 'F1_seg','recall_wi_truepc', 'precision_wi_truepc','F1_wi_truepc','nMCC_wi_truepc','FPR_wi_truepc','TPR_wi_truepc', 'acc_wi_truepc','bal_acc_wi_truepc','FP', 'FN', 'TN', 'TP')
  tib = as_tibble(matrix(nrow = 0, ncol = length(tbl_colnames)), .name_repair = ~ tbl_colnames)
  tbl2_colnames= c('mod', 'truelod', 'truesegl', 'minlod', 'minsegl', 'nea', 'popX', 'name', 'mean_prop_detected','mean_prop_true', 'detected_introgression','overlap_introgression','true_introgression')
  tib2 = as_tibble(matrix(nrow = 0, ncol = length(tbl2_colnames)), .name_repair = ~ tbl2_colnames)
}

if (length(ibd_popX) != 0){ 
  query <- ibd_popX

  ibd_popX <- as_tibble(ibd_popX)

  truelod <- min(ibd_popX$slod)
  truesegl <- min(ibd_popX$length)

  subject <- true_segments
  name <- unique(ibd_popX$ID)
  for (m in 1:length(name)){
          
    query_id <- query[query$ID == name[m]] %>% IRanges:::reduce()
    subject_id <- subject[subject$name == name[m]] %>% IRanges:::reduce()
          
    hits <- findOverlaps(query_id, subject_id)
    counts_queryinsubject <- countOverlaps(query_id, subject_id) # for two records in query there is only one in subject 
    overlaps <- get_overlap(query_id, subject_id, hits)
    conM_wi_sub <- get_confucionM_window_subject_based(overlaps, query_id,subject_id, hits)
    conM_seg <- get_confusionM_segcounts(query_id, subject_id, hits,overlaps)
    d <- c(mod, as.character(truelod), as.character(truesegl), as.character(minlod), as.character(minsegl), nea, pop2test, name[m])
          
    totvec = c(d, conM_seg, conM_wi_sub)
    ftab <- rbind(ftab, totvec)
    miss = get_missestimation(overlaps, query_id,subject_id, hits)
    misstab <- rbind(misstab, miss)
  }
  colnames(ftab)<- c('mod', 'truelod', 'truesegl', 'minlod', 'minsegl', 'nea', 'popX', 'name', 'recall_seg','precision_seg', 'F1_seg','recall_wi_truepc', 'precision_wi_truepc','F1_wi_truepc','nMCC_wi_truepc','FPR_wi_truepc','TPR_wi_truepc', 'acc_wi_truepc','bal_acc_wi_truepc','FP', 'FN', 'TN', 'TP')
  tib = as_tibble(ftab)
  colnames(misstab)<- c('mod', 'truelod', 'truesegl', 'minlod', 'minsegl', 'nea', 'popX', 'name', 'mean_prop_detected','mean_prop_true', 'detected_introgression','overlap_introgression','true_introgression')
  tib2 = as_tibble(misstab)
}

print("writing output...")
write_tsv(tib, outfile)
write_tsv(tib2, missout)
