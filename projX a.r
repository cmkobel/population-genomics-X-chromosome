library(tidyverse)
library(plyr)
library(zoo)
setwd("~/Biologi/Pop Gen/12 project/data")
# Population Genetics on X-chromosome 
#The data consists of a vcf file of 150 male full X chromosomes, a bed file with callable regions, a gif gene annotation file, a metafile with information about the samples and a set of files for use with REHH.
# OK, let's go through the file formats.

# Investigate the following


## A. Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions,
#including at least the contrast 
# between Africa and Europe, 
# between Europe and East Asia, and 
# between East Asia and Africa. Identify the 10 strongest Fst outlier regions in each case. Identify their genomic position and the genes covered by thse Fat peaks. Discuss potential adaptive explanations.
# AM
# Please note that this requires the res_scan_<region> dataframes in order to work. So run ## B first.

# FST: Jeg er ikke sikker på at dette er den rigtige måde at gøre det på. Da jeg ikke skalerer med antallet af individer


cmean = function(x) {
    if (sum(sapply(x, function(z) sum(length(which(is.na(z))))))/length(x) < 0.5) {
        return(mean(x, na.rm=T))
    }
    else {
        return(NA)
    }
    
}


fst_rolling_window = function(pos, freq_A_1, freq_A_2, size) {
    fst_table = tibble(
        pos = pos,
        p_1 = freq_A_1,
        q_1 = 1 - p_1,
        p_2 = freq_A_2,
        q_2 = 1 - p_2
    ) %>%
        mutate(H_S = (2*p_1*q_1 + 2*p_2*q_2) / 2) %>%
        mutate(H_T = 2 * (p_1 + p_2)/2 * (q_1 + q_2)/2) %>%
        mutate(F_ST = 1 - H_S / (H_T))# %>%
    #na.omit() # remove fixed loci # men der er jo ingen grund til at fjerne dem, når cmean alligevel godt kan håndtere NAs
    
    # bør valideres, evt. med width 3, så man manuelt kan krydstjekke?
    return(cbind(
        fst_table,
        sliding_window=rollapply(fst_table$F_ST, size, cmean, fill = NA)
        ))
}



af_we_fst_10 = fst_rolling_window(res_scan_AF$POSITION, res_scan_AF$freq_A, res_scan_WE$freq_A, 10)

af_we_fst = fst_rolling_window(res_scan_AF$POSITION, res_scan_AF$freq_A, res_scan_WE$freq_A, 100)
we_ea_fst = fst_rolling_window(res_scan_AF$POSITION, res_scan_WE$freq_A, res_scan_EA$freq_A, 100)
ea_af_fst = fst_rolling_window(res_scan_AF$POSITION, res_scan_EA$freq_A, res_scan_AF$freq_A, 100)


ggplot(af_we_fst_10, aes(x=pos)) +
    geom_line(aes(y=sliding_window)) +
    xlab("chromosome X position") + ylab("F_ST (10 SNP sliding window)") +
    ggtitle("AF/WE: F_ST")

ggplot(we_ea_fst, aes(x=pos)) +
    geom_line(aes(y=sliding_window)) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    ggtitle("WE/EA: F_ST")

ggplot(ea_af_fst, aes(x=pos)) +
    geom_line(aes(y=sliding_window)) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    ggtitle("EA/AF: F_ST")

# merged
ggplot(af_we_fst, aes(x=pos)) +
    geom_line(aes(y=af_we_fst$sliding_window, color="AF/WE")) +
    geom_line(aes(y=we_ea_fst$sliding_window, color="WE/EA")) +
    geom_line(aes(y=ea_af_fst$sliding_window, color="EA/AF")) +
    
    xlab("chromosome X position") + ylab("filtered F_ST (100 SNP sliding window)") +
    ggtitle("F_ST between regions")


# An interesting question is - how long back are you looking when the bin size is different?
# When the window size varies, the we vary the generation time we're looking back?




