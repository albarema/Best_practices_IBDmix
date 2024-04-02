quiet <- function(x) { suppressMessages(suppressWarnings(x)) }
quiet(library(tidyverse))

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

infile = args[1]
acfile <- args[2]
outfile <- args[3]
segl <- args[4]
mlod <- args[5]
model = args[6]

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

extract_digit_after_dot_s <- function(input_string) {
  # Use regular expression to match the digit after ".s"
  match_result <- regmatches(input_string, regexpr("\\.s(\\d+)", input_string))

  # Check if a match is found
  if (!is.na(match_result) && length(match_result) >= 1) {
    # Extract the matched digit
    digit_after_dot_s <-  gsub(".s", "", match_result)

    return(as.numeric(digit_after_dot_s))
  } else {
    # Return a default value if no match is found
    return(100)
  }
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

d = read_tsv(infile,show_col_types = FALSE) %>% filter(minlod == mlod, minsegl == segl)
d$samplesize <- lapply(d$popX, extract_digit_after_dot_s) %>% as.numeric()

tib = read_tsv(acfile,show_col_types = FALSE) %>% filter(minlod == mlod)
tib$samplesize <- lapply(tib$popX, extract_digit_after_dot_s) %>% as.numeric()

d$samplesize = factor(d$samplesize, levels=sort(as.numeric(unique(d$samplesize))))
tib$samplesize = factor(tib$samplesize, levels=sort(as.numeric(unique(tib$samplesize))))
d$minsegl <- as.character(d$minsegl)

get_summary <- function(statstab, smalltab, model, stats){
    smalltab = smalltab %>% mutate(overestimation=(detected_introgression - overlap_introgression)/true_introgression, underestimation = (true_introgression - overlap_introgression)/true_introgression)
    summ_o <- smalltab %>%
        group_by(samplesize) %>% 
        summarize(mean = mean(overestimation), median = median(overestimation), sd = sd(overestimation) ) %>% 
        mutate(stats = "overestimation")
    summ_u <- smalltab %>%
        group_by(samplesize) %>% 
        summarize(mean = mean(underestimation), median = median(underestimation), sd = sd(underestimation) )  %>% 
        mutate(stats = "underestimation")
    summ <- rbind(summ_o, summ_u)
    summ_tab = c()
    for (i in stats){
        summ2 <- statstab %>%
        group_by(samplesize) %>% 
        summarize(mean = mean(.data[[i]]), median = median(.data[[i]]), sd = sd(.data[[i]]))
        summ2$stats <- i
        summ_tab <- rbind(summ_tab, summ2)
    }
    tot_tab = rbind(summ,summ_tab)
    tot_tab$model = model  
    return(tot_tab)
}   

stats <- c("F1_wi_truepc", "nMCC_wi_truepc")
ftab = get_summary(d, tib, model, stats)


write_tsv(ftab, outfile)
