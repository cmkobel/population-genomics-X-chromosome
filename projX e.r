source("projX functions.r")

## E. Perform any additional analysis of your own choice, such as (diversity along the C X chromosome)

#Plot ancestral freq window
for(i in (c("WE", "AF", "EA", "SA", "AM", "CAS", "O"))) {
    print(
        load(paste("data/res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}



# I don't know what I'm doing here.
roll_window_mean = rollapply(res_scan_AF$freq_A, 1000, mean, fill=NA)
ggplot(res_scan_AF) + 
    geom_line(aes(x=POSITION, y=roll_window_mean, color="mean")) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)")
    


#   Linkage Disequilibrium
# ----------------------------   


# Make it more digital. 0s for ancestral alleles, 1s for derived allele, simply.
ancestral2zero = function(df_ancestral, df_snps) {
    cbind(df_ancestral, df_snps) %>%
    apply(1, function(snps) {
        zeroed = gsub(snps[1], 0, snps[-1]) # turn ancestral snps into zeroes
        as.integer(gsub("[A,C,T,G]", 1, zeroed)) # turn derived snps into ones
        # I thin there is a smarter way to do the two above lines. There is no reason why you wouldn't be able to do a case
    }) %>% 
    t() %>%
    return()
}


linkage_disequilibrium = function(binary_snp_data, title, size, plot = F) {
    # Parametriser
    # Dependencies
    #   snp_metadata is the metadata file
    #   binary_snps is the output of ancestral2zero with a specific genotype file from the populations.
    
    n_snps = dim(binary_snp_data)[1]
    n_snps_per_window = size
    n_folds = round(n_snps/n_snps_per_window)
    e_results = tibble(pos=NA, LD=NA)[-1,]
    if(n_snps > n_folds) {
    print(paste("fold size: ", n_snps/n_folds))
        for (fold in 1:n_folds) { # find a way to skip when the SD == 0
            start = floor((fold-1)*(n_snps/n_folds)+1)
            end = floor((fold)*(n_snps/n_folds))
            #print(paste(fold, ":", start, end))
            
            correlation = cor(t(binary_snp_data[start:end,]))**2 # transpose to get between-individuals corr.
            # Is there a nice way to catch NAs?
            diag(correlation) = NA
            
            #print(correlation)
            res_mean =  mean(correlation, na.rm = T)
            print(res_mean)
            
            e_results = rbind(e_results, tibble(pos = snp_metadata$pos[start + abs(start-end)/2], LD = res_mean))
        }
    }

    if (plot) { # plot or not
        ggplot(e_results) + 
            #geom_line(aes(x=pos, y=LD), size = 0.3) +
            geom_point(aes(x=pos, y=LD), size = 0.3) +
            xlab("chromosome X position") + ylab(paste("mean LD (", n_folds, " sequential windows with ", n_snps_per_window, " SNPs in each)", sep="")) +
            ggtitle(paste("Linkage disequilibrium calculated with r^2:", title))
        ggplot2::ggsave(paste("plots/e_LD/LD_rsq_", title, "_", n_snps_per_window,"_point.pdf", sep=""), width=10, height=6)
    }
    else return(e_results)
}


# In order to do LD we need the snp files.
snp_metadata <- read.delim("data/snps_filtered.inp", header=FALSE); names(snp_metadata) = c("name", "chromosome", "pos", "ancestral", "derived")

genotypes_WE <- read.delim("data/genotypes_WE.hap", header=FALSE)
genotypes_AF <- read.delim("data/genotypes_AF.hap", header=FALSE)
genotypes_EA <- read.delim("data/genotypes_EA.hap", header=FALSE)
genotypes_SA <- read.delim("data/genotypes_SA.hap", header=FALSE)
genotypes_AM <- read.delim("data/genotypes_AM.hap", header=FALSE)
genotypes_CAS <- read.delim("data/genotypes_CAS.hap", header=FALSE)
genotypes_O <- read.delim("data/genotypes_O.hap", header=FALSE)

binary_snps_WE = ancestral2zero(snp_metadata$ancestral, genotypes_WE)
binary_snps_AF = ancestral2zero(snp_metadata$ancestral, genotypes_AF)
binary_snps_EA = ancestral2zero(snp_metadata$ancestral, genotypes_EA)
binary_snps_SA = ancestral2zero(snp_metadata$ancestral, genotypes_SA)
binary_snps_AM = ancestral2zero(snp_metadata$ancestral, genotypes_AM)
binary_snps_CAS = ancestral2zero(snp_metadata$ancestral, genotypes_CAS)
binary_snps_O = ancestral2zero(snp_metadata$ancestral, genotypes_O)

# ignore following warnings, as they simply mean that the SD == 0 thus the corr. can't be calculated.
# Make LD plots
linkage_disequilibrium(binary_snps_WE, "WE", 100, T)
linkage_disequilibrium(binary_snps_AF, "AF", 100, T)
linkage_disequilibrium(binary_snps_EA, "EA", 100, T)
linkage_disequilibrium(binary_snps_SA, "SA", 100, T)
linkage_disequilibrium(binary_snps_AM, "AM", 100, T)
linkage_disequilibrium(binary_snps_CAS, "CAS", 100, T)
linkage_disequilibrium(binary_snps_O, "O", 100, T)

res_WE = linkage_disequilibrium(binary_snps_WE, "", 100)
res_AF = linkage_disequilibrium(binary_snps_AF, "", 100)
res_EA = linkage_disequilibrium(binary_snps_EA, "", 100)
res_SA = linkage_disequilibrium(binary_snps_SA, "", 100)
res_AM = linkage_disequilibrium(binary_snps_AM, "", 100)
res_CAS = linkage_disequilibrium(binary_snps_CAS, "", 100)
res_O = linkage_disequilibrium(binary_snps_O, "", 100)



plot(density(res_WE$LD, na.rm = T))
plot(density(res_AF$LD, na.rm = T))
plot(density(res_EA$LD, na.rm = T))
plot(density(res_SA$LD, na.rm = T))
plot(density(res_AM$LD, na.rm = T))
plot(density(res_CAS$LD, na.rm = T))
plot(density(res_O$LD, na.rm = T))



# Add LD density
