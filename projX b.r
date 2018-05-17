setwd("~/Biologi/Pop Gen/12 project/data")
source("../projX functions.r")

require(gridExtra) # for gridarrange
## B. Perform an iHS (integrated haplotype score) scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A.
# see week 7 for insp.

# Install the package
#install.packages('rehh')
library(rehh)

# Reading the data for each population:

# WE = WestEurasia
# AF = Africa
# EA = EastAsia
# SA = South Asia
# AM = America
# CAS = CentralAsiaSiberia
# O = Oceania

haplohh_WE = cd2h("genotypes_WE.hap", "snps_filtered.inp") # WestEurasia    44 haplotypes and 411892 SNPs
haplohh_AF = cd2h("genotypes_AF.hap", "snps_filtered.inp") # Africa         21 haplotypes and 411892 SNPs
haplohh_EA = cd2h("genotypes_EA.hap", "snps_filtered.inp") # EastAsia       24 haplotypes and 411892 SNPs
haplohh_SA = cd2h("genotypes_SA.hap", "snps_filtered.inp") # SouthAsia      30 haplotypes and 411892 SNPs
haplohh_AM = cd2h("genotypes_AM.hap", "snps_filtered.inp") # America         6 haplotypes and 411892 SNPs
haplohh_CAS = cd2h("genotypes_CAS.hap", "snps_filtered.inp") # CntrlAsiaSiber8 haplotypes and 411892 SNPs
haplohh_O = cd2h("genotypes_O.hap", "snps_filtered.inp") # Oceania          12 haplotypes and 411892 SNPs


# Compute ihh for all the snps in the halohh object considered.
# Computationally heavy
# res_scan_WE = scan_hh(haplohh_WE, threads = 4) # save this object to disk.
# res_scan_AF = scan_hh(haplohh_AF, threads = 4) # save this object to disk.
# res_scan_EA = scan_hh(haplohh_EA, threads = 4) # save this object to disk.
# res_scan_SA = scan_hh(haplohh_SA, threads = 4) # save this object to disk.
# res_scan_AM = scan_hh(haplohh_AM, threads = 4) # save this object to disk.
# res_scan_CAS = scan_hh(haplohh_CAS, threads = 4) # save this object to disk.
# res_scan_O = scan_hh(haplohh_O, threads = 4) # save this object to disk.
# save(res_scan_WE, file="res_scan_WE.rdata")
# save(res_scan_AF, file="res_scan_AF.rdata")
# save(res_scan_EA, file="res_scan_EA.rdata")
# save(res_scan_SA, file="res_scan_SA.rdata")
# save(res_scan_AM, file="res_scan_AM.rdata")
# save(res_scan_CAS, file="res_scan_CAS.rdata")
# save(res_scan_O, file="res_scan_O.rdata")


# ^ The above code products are loaded here (1000x faster than recomputing)
for(i in (c("WE", "AF", "EA", "SA", "AM", "CAS", "O"))) {
    print(
        load(paste("res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}

# Q2 Allele freq in pops.
pdf(file="freq_A_density_WE.pdf"); hist(res_scan_WE$freq_A); dev.off()
pdf(file="freq_A_density_AF.pdf"); hist(res_scan_AF$freq_A); dev.off()
pdf(file="freq_A_density_EA.pdf"); hist(res_scan_EA$freq_A); dev.off()
pdf(file="freq_A_density_SA.pdf"); hist(res_scan_SA$freq_A); dev.off()
pdf(file="freq_A_density_AM.pdf"); hist(res_scan_AM$freq_A); dev.off()
pdf(file="freq_A_density_CAS.pdf"); hist(res_scan_CAS$freq_A); dev.off()
pdf(file="freq_A_density_O.pdf"); hist(res_scan_O$freq_A); dev.off()

#Q3: How is the standardized iHH calculated? For what reason do they standardize iHS measure?

# Compute ihs (integrated haplotype score from integrated haplotype homozygosity)
wg_ihs_WE = ihh2ihs(res_scan_WE, freqbin = 0.16)
wg_ihs_AF = ihh2ihs(res_scan_AF, freqbin = 0.16)
wg_ihs_EA = ihh2ihs(res_scan_EA, freqbin = 0.16)
wg_ihs_SA = ihh2ihs(res_scan_SA, freqbin = 0.16)
wg_ihs_AM = ihh2ihs(res_scan_AM, freqbin = 0.16)
wg_ihs_CAS = ihh2ihs(res_scan_CAS, freqbin = 0.16)
wg_ihs_O = ihh2ihs(res_scan_O, freqbin = 0.16)

# Jeg ved ikke helt hvad freqbin gør, men nu har jeg øget den indtil den ikke brokker sig længere??


# Plotting the results.
# ihsplot(wg_ihs_WE, plot.pval = T, main="WE", cex = 0.005, pch = 19)
# ihsplot(wg_ihs_AF, plot.pval = T, main="AF")
# ihsplot(wg_ihs_EA, plot.pval = T, main="EA")
# ihsplot(wg_ihs_SA, plot.pval = T, main="SA")
# ihsplot(wg_ihs_AM, plot.pval = T, main="AM")
# ihsplot(wg_ihs_CAS, plot.pval = T, main="CAS")
# ihsplot(wg_ihs_O, plot.pval = T, main="O")

#Q4. Do you find outliers with significant iHS?
#Find the foverlap code


# set up plots in a fancy manner, better than ihsplot does 
plot_ihs = function(plot_df) {
    print(
        grid.arrange(
            ggplot(plot_df) + 
                geom_point(aes(x=POSITION, y=iHS),
                           size=0.05) +
                xlab("chromosome X position") +
                ylab("iHS") + ggtitle("AF"),
            
            ggplot(plot_df) + 
                geom_point(aes(x=POSITION, y=`-log10(p-value)`),
                           size=0.05) +
                xlab("chromosome X position") +
                ylab("p-value"), layout_matrix = rbind(c(1),c(2)))
    )
}

pdf("ihs_WE.pdf"); plot_ihs(wg_ihs_WE$iHS); dev.off()
pdf("ihs_AF.pdf"); plot_ihs(wg_ihs_AF$iHS); dev.off()
pdf("ihs_EA.pdf"); plot_ihs(wg_ihs_EA$iHS); dev.off()
pdf("ihs_SA.pdf"); plot_ihs(wg_ihs_SA$iHS); dev.off()
pdf("ihs_AM.pdf"); plot_ihs(wg_ihs_AM$iHS); dev.off()
pdf("ihs_CAS.pdf"); plot_ihs(wg_ihs_CAS$iHS); dev.off()
pdf("ihs_O.pdf"); plot_ihs(wg_ihs_O$iHS); dev.off()





# XP-EHH (sammenligning)
wg_xpehh_AM_CAS = ies2xpehh(res_scan_AM, res_scan_CAS, popname1 = "AM", popname2 = "CAS", method = "bilateral")
# Skal det her egentlig stå på den anden side af opgave C?

xpehhplot(wg_xpehh_AM_CAS, plot.pval = T)

# Hvis jeg får noget mere forståelse og kan gå lidt i dybden i noget af outputtet, kunne denne sektion godt være færdig.
