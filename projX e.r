source("projX functions.r")


## E. Perform any additional analysis of your own choice, such as (diversity along the C X chromosome)

#Plot ancestral freq window
for(i in (c("WE", "AF", "EA", "SA", "AM", "CAS", "O"))) {
    print(
        load(paste("res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}

roll_window_mean = rollapply(res_scan_AF$freq_A, 1000, mean, fill=NA)

ggplot(res_scan_AF) + 
    geom_line(aes(x=POSITION, y=roll_window_mean, color="mean")) +
    xlab("chromosome X position") + ylab("F_ST (1000 SNP sliding window)") +



    
    


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


linkage_disequilibrium = function(binary_snp_data, title, size) {
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
            
            correlation = cor(t(binary_snp_data[start:end,])) # transpose to get between-individuals corr.
            # Is there a nice way to catch NAs?
            diag(correlation) = NA
            
            #print(correlation)
            res_mean =  mean(correlation, na.rm = T)
            print(res_mean)
            
            e_results = rbind(e_results, tibble(pos = snp_metadata$pos[start + abs(start-end)/2], LD = res_mean))
        }
    }
    
    ggplot(e_results) + 
        #geom_line(aes(x=pos, y=LD), size = 0.3) +
        geom_point(aes(x=pos, y=LD), size = 0.3) +
        xlab("chromosome X position") + ylab(paste("mean LD (", n_folds, " sequential windows with ", n_snps_per_window, " SNPs in each)", sep="")) +
        ggtitle(paste("Linkage disequilibrium:", title))
    ggplot2::ggsave(paste("plots/e_LD/LD_", title, "_", n_snps_per_window,"_point.pdf", sep=""), width=10, height=6)
    
    #return(e_results)
}


# In order to do LD we need the snp files.
snp_metadata <- read.delim("~/Biologi/Pop Gen/12 project/data/snps_filtered.inp", header=FALSE); names(snp_metadata) = c("name", "chromosome", "pos", "ancestral", "derived")

genotypes_WE <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_WE.hap", header=FALSE)
genotypes_AF <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_AF.hap", header=FALSE)
genotypes_EA <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_EA.hap", header=FALSE)
genotypes_SA <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_SA.hap", header=FALSE)
genotypes_AM <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_AM.hap", header=FALSE)
genotypes_CAS <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_CAS.hap", header=FALSE)
genotypes_O <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_O.hap", header=FALSE)


binary_snps_WE = ancestral2zero(snp_metadata$ancestral, genotypes_WE)
binary_snps_AF = ancestral2zero(snp_metadata$ancestral, genotypes_AF)
binary_snps_EA = ancestral2zero(snp_metadata$ancestral, genotypes_EA)
binary_snps_SA = ancestral2zero(snp_metadata$ancestral, genotypes_SA)
binary_snps_AM = ancestral2zero(snp_metadata$ancestral, genotypes_AM)
binary_snps_CAS = ancestral2zero(snp_metadata$ancestral, genotypes_CAS)
binary_snps_O = ancestral2zero(snp_metadata$ancestral, genotypes_O)

# ignore the following warnings, as they simply mean that the SD == 0 thus the corr. can't be calculated.
linkage_disequilibrium(binary_snps_WE, "WE", 100)
linkage_disequilibrium(binary_snps_AF, "AF", 100)
linkage_disequilibrium(binary_snps_EA, "EA", 100)
linkage_disequilibrium(binary_snps_SA, "SA", 100)
linkage_disequilibrium(binary_snps_AM, "AM", 100)
linkage_disequilibrium(binary_snps_CAS, "CAS", 100)
linkage_disequilibrium(binary_snps_O, "O", 100)



# Add LD density
    


