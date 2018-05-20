# Author: Carl M. Kobel

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



# PLAN A
# Denne funktion bruges ikke hvis sorteringsmeoden er frugtbar
cut_fst = function(df_fst, n, return = FALSE) {
    sorted = sort(df_fst$sliding_window, index.return = T, decreasing = T)
    # AF_WE, 7000 seems nice
    result = cbind(c(365712, df_fst$pos[sorted$ix[1:n]], 155110877), c(0.0, sorted$x[1:n], 0.0))
    if (return) return(result)
    else print(plot(result, cex = 0.3))
}



# AF WE
#cut_fst(fst_af_we, 5000)
# Comment on window size: Because we want to finde parts of genes having high Fst, we need to look in windows 
af_we_thresh = 310; cut_fst(fst_af_we_100, af_we_thresh)
#write.table(cut_fst(fst_af_we_100, af_we_thresh, T), file = "10_peaks_fst_af_we.tsv")
af_we_peaks = tibble(start = c(19182660, 37878268, 46095748, 55959471, 66270950, 92414472, 104551933, 140734729, 141641747, 145276227),
                     end = start,
                     fst = c(0.2502342793, 0.2194004233, 0.2158843816, 0.2101256287, 0.2144259546, 0.2463066009, 0.2330120671, 0.2172732175, 0.2584422211, 0.214525942))
plot(af_we_peaks)


# WE EA <
we_ea_thresh = 700; cut_fst(fst_we_ea_100, we_ea_thresh) # plot
#write.table(cut_fst(fst_we_ea_100, we_ea_thresh, T), file = "10_peaks_fst_we_ea.tsv")
we_ea_peaks = tibble(start = c(25671846, 42579633, 71373407, 73295813, 74119733, 87742329, 108285879, 109355662, 116617318, 128684757),
                     end = start,
                      fst = c(0.2626619986, 0.2428685249, 0.2935684424, 0.3036629517, 0.2230540469, 0.2677614183, 0.2119459508, 0.2461968528, 0.2129359005, 0.2161413781))
plot(we_ea_peaks)


# EA AF
ea_af_thresh = 870; cut_fst(fst_ea_af_100, ea_af_thresh)
#write.table(cut_fst(fst_ea_af_100, ea_af_thresh, T), file = "10_peaks_fst_ea_af.tsv")
ea_af_peaks = tibble(start = c(19183926, 36636363, 55960598, 63969563, 66270950, 90576998, 112894056, 121069027, 124687575, 126784188),
                     end = start,
                     fst = c(0.2697129718, 0.3670665736, 0.2939780025, 0.2757454847, 0.390637325, 0.3747966217, 0.2606229121, 0.2809754149, 0.2643502852, 0.3067912248))
plot(ea_af_peaks[2:3])



# PLAN B
#testet:

# AF WE

# original kode
sorted = sort(fst_af_we_100$sliding_window, index.return = T, decreasing = T)
sortedn = tibble(x = sorted$x, ix = sorted$ix + 49, pos = fst_af_we_100$pos[sorted$ix + 49]) # adjust for window size
# nu kan det her indsættet i overlap med et range fra 1:n som giver et ønsket antal gener.
# Det kan godt være at koden kan skrives mere effektivt med argumentet partial (?sort), jeg gider bare ikke teste det.

old_af_we_all_peak_candidates = tibble(start = sortedn$pos,
                                   end = sortedn$pos,
                                   fst_rolwin_hecto = sortedn$x)

# parametriseret udgave
sort_fst = function(fst_region_100) {
    # fst_region_100 needs to have a column called sliding_window with fst values offset 49 rows, and a pos column with the respective genome positions for fst values
    sorted = sort(fst_region_100$sliding_window,
                  index.return = T,
                  decreasing = T)
    sortedn = tibble(x = sorted$x,
                     ix = sorted$ix + 49,
                     pos = fst_region_100$pos[sorted$ix + 49]) # adjust for window size
    # nu kan det her indsættet i overlap med et range fra 1:n som giver et ønsket antal gener.
    # Det kan godt være at koden kan skrives mere effektivt med argumentet partial (?sort), jeg gider bare ikke teste det.
    
    return_variable = tibble(start = sortedn$pos,
                          end = sortedn$pos,
                          fst_rolwin_hecto = sortedn$x)
    return(return_variable)
}
af_we_all_peak_candidates = sort_fst(fst_af_we_100)
we_ea_all_peak_candidates = sort_fst(fst_we_ea_100)
ea_af_all_peak_candidates = sort_fst(fst_ea_af_100)


#   III:
# ---------
# Identify their genomic position and the genes covered by thse Fat peaks.

fst_overlap = function(fsts, gene_annotation) {
    # Example of formatting of input, start and end is needed for both
    # fsts = data.table(start = c(12, 55, 19), end = c(12, 55, 19), fst = c(0.1, 0.2, 0.3))
    # gene_annotation = data.table(start = c(10,30,50,70), end = c(20,40,60,80), x_gen = c("gen1", "gen2", "gen3", "gen4"))
    setkey(fsts, start, end)
    return(
        foverlaps(gene_annotation, fsts, type="any", nomatch = 0)
    )
}


# import annotation
# require(GenomicRanges) #? not needed, right?
gtf <- as.data.table(rtracklayer::import("data/gencode.v17.annotation.gtf"))
gtf <- gtf[gtf$seqnames == 'chrX']   # to select only X chromosome


# PLAN A
# gammelt overlap fra manuelle fundne peaks.
overlap_af_we = fst_overlap(as.data.table(af_we_peaks), gtf)
overlap_af_we
View(overlap_af_we)


# PLAN B

# AF WE
lim = 12000
new_overlap_af_we = fst_overlap(as.data.table(af_we_all_peak_candidates[1:lim,]), gtf[1:lim,])
simple = new_overlap_af_we[,c(1, 2, 3, 15, 17, 18, 23)]
unique(new_overlap_af_we$gene_name)


# WE EA
lim = 12000
new_overlap_we_ea = fst_overlap(as.data.table(we_ea_all_peak_candidates[1:lim,]), gtf[1:lim,])
simple = new_overlap_we_ea[,c(1, 2, 3, 15, 17, 18, 23)]
unique(new_overlap_we_ea$gene_name)


# EA AF
lim = 12000
new_overlap_ea_af = fst_overlap(as.data.table(ea_af_all_peak_candidates[1:lim,]), gtf[1:lim,])
simple = new_overlap_ea_af[,c(1, 2, 3, 15, 17, 18, 23)]
unique(new_overlap_ea_af$gene_name)




#   IV:
# ------
# Discuss potential adaptive explanations.



# An interesting question is - how long back are you looking when the bin size is different?
# When the window size varies, the we vary the generation time we're looking back?


# add the distribution


# density of snps



