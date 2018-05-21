# Author: Carl M. Kobel
plot = F
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


# load in the scan_hh products
for(i in (c("WE", "AF", "EA", "SA", "AM", "CAS", "O"))) {
    print(
        load(paste("data/res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}


# Define the function that calculates Fst and makes the rolling window
fst_rolling_window = function(pos, freq_A_1, freq_A_2, size, mfunction) {
    # Fst part
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
    
    # Rolling window part
    return(cbind(
        fst_table,
        sliding_window=rollapply(fst_table$F_ST, size, mfunction, fill = NA)
        ))
}

# FST: Jeg er ikke sikker på at dette er den rigtige måde at gøre det på. Da jeg ikke skalerer med antallet af individer
# window 1000, only needed for plotting
if (plot) {
fst_af_we = fst_rolling_window(res_scan_AF$POSITION, res_scan_AF$freq_A, res_scan_WE$freq_A, 1000, cmean)
fst_we_ea = fst_rolling_window(res_scan_AF$POSITION, res_scan_WE$freq_A, res_scan_EA$freq_A, 1000, cmean)
fst_ea_af = fst_rolling_window(res_scan_AF$POSITION, res_scan_EA$freq_A, res_scan_AF$freq_A, 1000, cmean)
}


# window 100, for statistics
fst_af_we_100 = fst_rolling_window(res_scan_AF$POSITION, res_scan_AF$freq_A, res_scan_WE$freq_A, 100, cmean)
fst_we_ea_100 = fst_rolling_window(res_scan_AF$POSITION, res_scan_WE$freq_A, res_scan_EA$freq_A, 100, cmean)
fst_ea_af_100 = fst_rolling_window(res_scan_AF$POSITION, res_scan_EA$freq_A, res_scan_AF$freq_A, 100, cmean)


# # leg med densitet og percentiler
# quantile(fst_af_we_100$sliding_window, 0.95, na.rm = T)
# 
# # Ville det ikke være smart at smide densitets plot ind, måske nede i plan B2?
# plot(density(na.omit(fst_af_we_100$sliding_window)))
# abline(v=0.1607412)
# abline(v=0.1172794)



if (plot) {
# add hline
pdf("plots/a_fst/fst_AF_WE.pdf", width = 10, height = 6)
ggplot(fst_af_we, aes(x = pos)) + 
    geom_point(aes(y = sliding_window), size = 0.03) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    #geom_hline(aes(yintercept=1, linetype=cutoff), data=cutoff, show.legend=F) +
    ggtitle("Fst: AF/WE")
    #ggplot2::ggsave(paste("plots/a_fst/fst_AF_WE.png"))
dev.off()


pdf("plots/a_fst/fst_WE_EA.pdf", width = 10, height = 6)
ggplot(fst_we_ea, aes(x = pos)) +
    geom_point(aes(y = sliding_window), size = 0.03) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    ggtitle("Fst: WE/EA")
    #ggplot2::ggsave(paste("plots/a_fst/fst_WE_EA.png"))
dev.off()


pdf("plots/a_fst/fst_WE_AF.pdf", width = 10, height = 6)
ggplot(fst_ea_af, aes(x = pos)) +
    geom_point(aes(y = sliding_window), size=0.03) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +
    ggtitle("Fst: EA/AF")
    #ggplot2::ggsave(paste("plots/a_fst/fst_EA_AF.png"))
dev.off()
    
    
# merged
pdf("plots/a_fst/fst_all.pdf", width = 10, height = 6)
ggplot(fst_af_we, aes(x=pos)) +
    geom_point(aes(y=fst_af_we$sliding_window, color="AF/WE"), size = 0.02, alpha = 0.6) +
    geom_point(aes(y=fst_we_ea$sliding_window, color="WE/EA"), size = 0.02, alpha = 0.6) +
    geom_point(aes(y=fst_ea_af$sliding_window, color="EA/AF"), size = 0.02, alpha = 0.6) +
    xlab("chromosome X position") + ylab("filtered F_ST (1000 SNP sliding window)") +
    ggtitle("F_ST between regions")
    #ggplot2::ggsave(paste("plots/a_fst/fst_all.png"))
dev.off()
}

#   II:
# --------
# Identify the 10 strongest Fst outlier regions in each case. (This will be done from the 100 snp windows)



# sort_fst behøver ikke at sortere, men det er meget rart at have sliding window og positioner (star = end) sammen
not_sort_fst = function(fst_region_100) {
    return(
        na.omit(
            tibble(start = fst_region_100$pos,
                   end = start,
                   fst_rolwin_hecto = fst_region_100$sliding_window)
            )
        )
}

af_we_all_peak_candidates = not_sort_fst(fst_af_we_100)
we_ea_all_peak_candidates = not_sort_fst(fst_we_ea_100)
ea_af_all_peak_candidates = not_sort_fst(fst_ea_af_100)


#   III:
# ---------
# Identify their genomic position and the genes covered by thse Fat peaks.

fst_overlap = function(fsts, gene_annotation) {
    setkey(fsts, start, end)
    return(
        foverlaps(gene_annotation, fsts, type="any", nomatch = 0)
    )
}


# import annotation
# require(GenomicRanges) #? not needed, right?
gtf <- as.data.table(rtracklayer::import("data/gencode.v17.annotation.gtf"))
gtf <- gtf[gtf$seqnames == 'chrX']   # to select only X chromosome



# AF WE nyeste plan
# percentile = 0.998
# threshold = quantile(af_we_all_peak_candidates$fst_rolwin_hecto, percentile)
# new_overlap_af_we = fst_overlap(as.data.table(af_we_all_peak_candidates[af_we_all_peak_candidates$fst_rolwin_hecto >= threshold,]), gtf) # fst_af_we_100 burde kunne bruges lige så vel som reg_all_peak_cand.
# #View(new_overlap_af_we)
# unique(cbind(new_overlap_af_we$gene_name))
# 
# peak_result = new_overlap_af_we %>%
#     group_by(gene_name) %>%
#     summarise(position = start[which.max(fst_rolwin_hecto)],
#               fst_peak_nui = max(fst_rolwin_hecto),
#               transcript_type = transcript_type[which.max(fst_rolwin_hecto)]) # inserted in text


# parametriser:
get_peaks = function(fst_for_overlap, fst_column, percentile) {
    threshold = quantile(fst_column, percentile)
    new_overlap_region = fst_overlap(as.data.table(fst_for_overlap[fst_column >= threshold,]), gtf) # fst_af_we_100 burde kunne bruges lige så vel som reg_all_peak_cand.
    #View(new_overlap_region)
    unique(cbind(new_overlap_region$gene_name))
    
    peak_result = new_overlap_region %>%
        group_by(gene_name) %>%
        summarise(position = start[which.max(fst_rolwin_hecto)],
                  fst_peak = max(fst_rolwin_hecto),
                  transcript_type = transcript_type[which.max(fst_rolwin_hecto)]) # inserted in text
    
    return(peak_result)
    
    
    
}



get_peaks_af_we = get_peaks(af_we_all_peak_candidates, af_we_all_peak_candidates$fst_rolwin_hecto, 0.998)
get_peaks_af_we = get_peaks_af_we[order(get_peaks_af_we$fst_peak, decreasing = T),]
View(get_peaks_af_we)




#   IV:
# ------
# Discuss potential adaptive explanations.

# AF WE
# 
# GPR143 FALSE
# https://www.ncbi.nlm.nih.gov/gene/4935
# "binds to heterotrimeric G proteins and is targeted to melanosomes in pigment cells."

# WE EA


# EA AF


# An interesting question is - how long back are you looking when the bin size is different?
# When the window size varies, the we vary the generation time we're looking back?


# add the distribution


# density of snps



