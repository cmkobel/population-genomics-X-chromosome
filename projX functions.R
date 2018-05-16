library(tidyverse)
library(zoo)
library(readr)



cmean = function(x) {
    return(mean(x, na.rm=T))
}

cmean_filtered = function(x, threshold) {
    if (sum(sapply(x, function(z) sum(length(which(is.na(z))))))/length(x) < 0.5) {
        return(mean(x, na.rm=T))
    }
    else {
        return(NA)
    }
    
}

cmedian = function(x) {
    return(median(x, na.rm=T))
}

cd2h = function(genotypes, snps) {
    return(
        data2haplohh(hap_file=genotypes, map_file=snps,
                     recode.allele=TRUE,
                     min_perc_geno.snp=100,
                     min_perc_geno.hap=80,
                     haplotype.in.columns=TRUE, chr.name=1)
    )
}
