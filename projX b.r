source("projX functions.r")

require(gridExtra) # for gridarrange in plotting
#   B: Extended haplotype scoring
#------------------------------------
# # See week 7 for insp.


#   I
# -----
# Perform an iHS (integrated haplotype score) scan of the whole X chromosome for at least three populations. 

# Install the package
library(rehh)

# Reading the data for each population:

# WE = WestEurasia
# AF = Africa
# EA = EastAsia
# SA = South Asia
# AM = America
# CAS = CentralAsiaSiberia
# O = Oceania





haplohh_AF = cd2h("data/genotypes_AF.hap", "data/snps_filtered.inp") # Africa         21 haplotypes and 411892 SNPs
haplohh_WE = cd2h("data/genotypes_WE.hap", "data/snps_filtered.inp") # WestEurasia    44 haplotypes and 411892 SNPs
haplohh_EA = cd2h("data/genotypes_EA.hap", "data/snps_filtered.inp") # EastAsia       24 haplotypes and 411892 SNPs

# Look at random positions
# andet argument er snp nummeret.
site_specific_ehh = calc_ehhs(haplohh_WE, 2)


# Compute ihh for all the snps in the halohh object considered.
# Computationally heavy
# res_scan_AF = scan_hh(haplohh_AF, threads = 4) # save this object to disk.
# res_scan_WE = scan_hh(haplohh_WE, threads = 4) # save this object to disk.
# res_scan_EA = scan_hh(haplohh_EA, threads = 4) # save this object to disk.

# save(res_scan_AF, file="res_scan_AF.rdata")
# save(res_scan_WE, file="res_scan_WE.rdata")
# save(res_scan_EA, file="res_scan_EA.rdata")

# ^ The above code products are loaded here (1000x faster than recomputing)
#for(i in (c("WE", "AF", "EA", "SA", "AM", "CAS", "O"))) {
for(i in (c("AF","WE", "EA"))) {
        
    print(
        load(paste("data/res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}

# Q2 Allele freq in pops.
pdf(file="../plots/b_ehh/freq_A_density_AF.pdf"); hist(res_scan_AF$freq_A); dev.off()
pdf(file="../plots/b_ehh/freq_A_density_WE.pdf"); hist(res_scan_WE$freq_A); dev.off()
pdf(file="../plots/b_ehh/freq_A_density_EA.pdf"); hist(res_scan_EA$freq_A); dev.off()

#Q3: How is the standardized iHH calculated? For what reason do they standardize iHS measure?

# Compute ihs (integrated haplotype score from integrated haplotype homozygosity)
wg_ihs_AF = ihh2ihs(res_scan_AF, freqbin = 0.16)
wg_ihs_WE = ihh2ihs(res_scan_WE, freqbin = 0.16)
wg_ihs_EA = ihh2ihs(res_scan_EA, freqbin = 0.16)
# Jeg ved ikke helt hvad freqbin gør, men nu har jeg øget den indtil den ikke brokker sig længere??
# Jeg går ud fra at freqbin bruges til oversættelsen til signifikans.


# Plotting the results.
# ihsplot(wg_ihs_AF, plot.pval = T, main="AF")
# ihsplot(wg_ihs_WE, plot.pval = T, main="WE", cex = 0.005, pch = 19)
# ihsplot(wg_ihs_EA, plot.pval = T, main="EA")

#Q4. Do you find outliers with significant iHS?
#Find the foverlap code


# set up plots in a fancy manner, better than ihsplot does 
plot_ihs = function(plot_df, title) {
    print(
        grid.arrange(
            ggplot(plot_df) + 
                geom_point(aes(x=POSITION, y=iHS), size=0.03) +
                xlab("") +
                ylab("iHS") + ggtitle(title),
            ggplot(plot_df) + 
                geom_point(aes(x=POSITION, y=`-log10(p-value)`), size=0.03) +
                xlab("chromosome X position") +
                ylab("-log10( p-value )"),
            layout_matrix = rbind(c(1),c(2))
        )
    )
}

#pdf("relation ihs pval.pdf"); plot(wg_ihs_WE$iHS$iHS, wg_ihs_WE$iHS$`-log10(p-value)`); dev.off()

pdf("../plots/b_ehh/ihs_WE.pdf"); plot_ihs(wg_ihs_WE$iHS, "WE"); dev.off()
pdf("../plots/b_ehh/ihs_AF.pdf"); plot_ihs(wg_ihs_AF$iHS, "AF"); dev.off()
pdf("../plots/b_ehh/ihs_EA.pdf"); plot_ihs(wg_ihs_EA$iHS, "EA"); dev.off()
pdf("../plots/b_ehh/ihs_SA.pdf"); plot_ihs(wg_ihs_SA$iHS, "SA"); dev.off()
pdf("../plots/b_ehh/ihs_AM.pdf"); plot_ihs(wg_ihs_AM$iHS, "AM"); dev.off()
pdf("../plots/b_ehh/ihs_CAS.pdf"); plot_ihs(wg_ihs_CAS$iHS, "CAS"); dev.off()
pdf("../plots/b_ehh/ihs_O.pdf"); plot_ihs(wg_ihs_O$iHS,"O"); dev.off()



#   II
# ------
# Identify the 10 most significant regions and associated with genes as in A.

# Denne funktion skal være generel og kan egentlig indsættes i projX functions.r

overlap = function(candidates, gene_annotation) {
    setkey(candidates, start, end)
    return(
        foverlaps(gene_annotation, candidates, type="any", nomatch = 0)
    )
}

gtf <- as.data.table(rtracklayer::import("data/gencode.v17.annotation.gtf"))
gtf <- gtf[gtf$seqnames == 'chrX']   # to select only X chromosome

# sorter en percentil
# skal parametriseres
# sorter med hele lortet denne gang.
# what = function(sorter_dette, sorter_efter_kolonne_nr, percentil) {
sorted = wg_ihs_AF$iHS[order(wg_ihs_AF$iHS$`-log10(p-value)`, decreasing = T),]
best_percentile = quantile(sorted$`-log10(p-value)`, 0.99976, na.rm = T)[[1]] # beregn percentil før filtrering
best_percentile
plot(density(
    na.omit(
        sorted$`-log10(p-value)`
    )
)); abline(v = best_percentile)
sorted_filtered = sorted[sorted$`-log10(p-value)` >= best_percentile,] # fjern alt under en tærskel
sorted_filtered_tibble = tibble(start = sorted_filtered$POSITION,
                                end = start,
                                ihs = sorted_filtered$iHS,
                                `-log10(p-value)` = sorted_filtered$`-log10(p-value)`)
overlap_out = overlap(na.omit(as.data.table(sorted_filtered_tibble)), gtf)
unique(overlap_out$gene_name)
# så er der ti !

to_overlap = tibble(start = wg_ihs_AF$iHS$POSITION,
                    end = start,
                    iHS = wg_ihs_AF$iHS$iHS,
                    ppval = wg_ihs_AF$iHS$`-log10(p-value)`) %>% 
             na.omit()


get_peaks = function(to_overlap, value_column, percentile) {
    # 1 c(start, end, value)
    # 2 c(value)
    # 3 percentile
    threshold = quantile(value_column, percentile)
    new_overlap_region = overlap(as.data.table(to_overlap[value_column >= threshold,]), gtf) # fst_af_we_100 burde kunne bruges lige så vel som reg_all_peak_cand.
    #return(new_overlap_region)
    #View(new_overlap_region)
    #return(unique(cbind(new_overlap_region$gene_name)))
    
    peak_result = new_overlap_region %>%
        group_by(gene_name) %>%
        summarise(position = start[which.max(iHS)],
                  iHS = iHS[which.max(iHS)],
                  ppval = ppval[which.max(iHS)],
                  gene_type = gene_type[which.max(iHS)],
                  comment = "")
                  
                  

    #               fst_peak = max(fst_rolwin_hecto),
    #               comment = "no comment") # inserted in text
    # 
    return(peak_result)
    
}

get_peaks_af = get_peaks(to_overlap, to_overlap$ppval, 0.99976)
get_peaks_af = get_peaks_af[order(get_peaks_af$ppval, decreasing = T),]
#View(get_peaks_af_we)


# Hvis jeg får noget mere forståelse og kan gå lidt i dybden i noget af outputtet, kunne denne sektion godt være færdig.
