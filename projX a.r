source("projX functions.r")
# Population Genetics on X-chromosome 
# The data consists of a vcf file of 150 male full X chromosomes, a bed file with callable regions, a gif gene annotation file, a metafile with information about the samples and a set of files for use with REHH.
# Investigate the following

##  A 
#-------

#   I: 
# ------
# Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions, including at least the contrast 
# between Africa and Europe, 
# between Europe and East Asia, and 
# between East Asia and Africa. 


for(i in (c("WE", "AF", "EA", "SA", "AM", "CAS", "O"))) {
    print(
        load(paste("data/res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}



fst_rolling_window = function(pos, freq_A_1, freq_A_2, size, mfunction) {
    fst_table = tibble(
        pos = pos,
        p_1 = freq_A_1,
        q_1 = 1 - p_1,
        p_2 = freq_A_2,
        q_2 = 1 - p_2
    ) %>%
        mutate(H_S = (2*p_1*q_1 + 2*p_2*q_2) / 2) %>%
        mutate(H_T = 2 * ((p_1 + p_2)/2) * ((q_1 + q_2)/2)) %>% 
        mutate(F_ST = 1 - H_S / (H_T))# %>%
    #na.omit() # nej, NA's skal jo fjernes så sent som muligt. Ellers ændres vinduestørrelsen (fejlagtigt) omvendt prop. med 
    
    # bør valideres, evt. med width 3, så man manuelt kan krydstjekke?
    return(cbind(
        fst_table,
        sliding_window=rollapply(fst_table$F_ST, size, mfunction, fill = NA)
        ))
}

# FST: Jeg er ikke sikker på at dette er den rigtige måde at gøre det på. Da jeg ikke skalerer med antallet af individer

fst_af_we = fst_rolling_window(res_scan_AF$POSITION, res_scan_AF$freq_A, res_scan_WE$freq_A, 1000, cmean)
fst_we_ea = fst_rolling_window(res_scan_AF$POSITION, res_scan_WE$freq_A, res_scan_EA$freq_A, 1000, cmean)
fst_ea_af = fst_rolling_window(res_scan_AF$POSITION, res_scan_EA$freq_A, res_scan_AF$freq_A, 1000, cmean)

fst_af_we_100 = fst_rolling_window(res_scan_AF$POSITION, res_scan_AF$freq_A, res_scan_WE$freq_A, 100, cmean)
fst_we_ea_100 = fst_rolling_window(res_scan_AF$POSITION, res_scan_WE$freq_A, res_scan_EA$freq_A, 100, cmean)
fst_ea_af_100 = fst_rolling_window(res_scan_AF$POSITION, res_scan_EA$freq_A, res_scan_AF$freq_A, 100, cmean)


# add hline
ggplot(fst_af_we, aes(x = pos)) +
    geom_point(aes(y = sliding_window), size = 0.05, alpha = 0.5) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    #geom_hline(aes(yintercept=1, linetype=cutoff), data=cutoff, show.legend=F) +
    ggtitle("Fst: AF/WE")
    ggplot2::ggsave(paste("plots/a_fst/fst_AF_WE.png"))

ggplot(fst_we_ea, aes(x = pos)) +
    geom_point(aes(y = sliding_window), size = 0.05, alpha = 0.5) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    ggtitle("Fst: WE/EA")
    ggplot2::ggsave(paste("plots/a_fst/fst_WE_EA.png"))

ggplot(fst_ea_af, aes(x = pos)) +
    geom_point(aes(y = sliding_window), size=0.05, alpha = 0.5) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    ggtitle("Fst: EA/AF")
    ggplot2::ggsave(paste("plots/a_fst/fst_EA_AF.png"))

# merged
ggplot(fst_af_we, aes(x=pos)) +
    geom_point(aes(y=fst_af_we$sliding_window, color="AF/WE"), size = 0.005) +
    geom_point(aes(y=fst_we_ea$sliding_window, color="WE/EA"), size = 0.005) +
    geom_point(aes(y=fst_ea_af$sliding_window, color="EA/AF"), size = 0.005) +
    xlab("chromosome X position") + ylab("filtered F_ST (1000 SNP sliding window)") +
    ggtitle("F_ST between regions")


#   II:
# --------
# Identify the 10 strongest Fst outlier regions in each case. (This will be done from the 100 snp windows)

sorted = sort(fst_af_we$sliding_window, index.return = T, decreasing = T)

# AF_WE, 7000 seems nice
n = 5000

plot(c(0, fst_af_we$pos[sorted$ix[1:n]]), c(0, sorted$x[1:n]))

# OK nu tager jeg de her index og sorterer dem
positions = sort(fst_af_we$pos[sorted$ix[1:n]], index.return = T)
plot(positions$x, positions$ix)

out = cbind(c(0, fst_af_we$pos[sorted$ix[1:n]]), c(0, sorted$x[1:n])) # and into excel
write.table(out[,1], file="thybajer1.tsv", sep=" ", row.names = F)
write.table(out[,2], file="thybajer2.tsv", sep=" ", row.names = F)

library(quantmod)
findPeaks(sin(1:23))
    


findPeaks(sin(1:10))

p <- findPeaks(sin(seq(1,10,.1)))
sin(seq(1,10,.1))[p]

plot(sin(seq(1,10,.1))[p])
plot(sin(seq(1,10,.1)),type='l')
points(p,sin(seq(1,10,.1))[p])


#   III:
# ---------
# Identify their genomic position and the genes covered by thse Fat peaks.


#   IV:
# -----
# Discuss potential adaptive explanations.



# An interesting question is - how long back are you looking when the bin size is different?
# When the window size varies, the we vary the generation time we're looking back?


# add the distribution


# density of snps



