## B. Perform an iHS (integrated haplotype score) scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A.
# see week 7 for insp.

# Install the package
#install.packages('rehh')
# library(rehh)

# cd2h = function(genotypes, snps) {
#     return(
#         data2haplohh(hap_file=genotypes ,map_file=snps,
#                      recode.allele=TRUE,
#                      min_perc_geno.snp=100,
#                      min_perc_geno.hap=80,
#                      haplotype.in.columns=TRUE, chr.name=1)
#     )
# }



# Reading the data for each population:

# WE = WestEurasia
# AF = Africa
# EA = EastAsia
# SA = South Asia
# AM = America
# CAS = CentralAsiaSiberia
# O = Oceania

# haplohh_WE = cd2h("genotypes_WE.hap", "snps_filtered.inp") # WestEurasia    44 haplotypes and 411892 SNPs
# haplohh_AF = cd2h("genotypes_AF.hap", "snps_filtered.inp") # Africa         21 haplotypes and 411892 SNPs
# haplohh_EA = cd2h("genotypes_EA.hap", "snps_filtered.inp") # EastAsia       24 haplotypes and 411892 SNPs
# haplohh_SA = cd2h("genotypes_SA.hap", "snps_filtered.inp") # SouthAsia      30 haplotypes and 411892 SNPs
# haplohh_AM = cd2h("genotypes_AM.hap", "snps_filtered.inp") # America         6 haplotypes and 411892 SNPs
# haplohh_CAS = cd2h("genotypes_CAS.hap", "snps_filtered.inp") # CntrlAsiaSiber8 haplotypes and 411892 SNPs
# haplohh_O = cd2h("genotypes_O.hap", "snps_filtered.inp") # Oceania          12 haplotypes and 411892 SNPs



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


## hejse
## load again
for(i in (c("WE", "AF", "EA", "SA", "AM", "CAS", "O"))) {
    print(
        load(paste("res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}


# Compute ihs (integrated haplotype score from integrated haplotype homozygosity)
wg_ihs_AM = ihh2ihs(res_scan_AM, freqbin = 0.16)
wg_ihs_CAS = ihh2ihs(res_scan_CAS, freqbin = 0.16)
# Jeg ved ikke helt hvad freqbin gør, men nu har jeg øget den indtil den ikke brokker sig længere.


# Plotting the results.
ihsplot(wg_ihs_AM, plot.pval = T)
ihsplot(wg_ihs_CAS, plot.pval = T)



# XP-EHH (sammenligning)
wg_xpehh_AM_CAS = ies2xpehh(res_scan_AM, res_scan_CAS, popname1 = "AM", popname2 = "CAS", method = "bilateral")
# Skal det her egentlig stå på den anden side af opgave C?

xpehhplot(wg_xpehh_AM_CAS, plot.pval = T)

# Hvis jeg får noget mere forståelse og kan gå lidt i dybden i noget af outputtet, kunne denne sektion godt være færdig.
