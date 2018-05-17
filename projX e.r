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


# Jeg vil gerne beregne LD.




# in order to do LD we need the snp files.
genotypes_AF <- read.delim("~/Biologi/Pop Gen/12 project/data/genotypes_AF.hap", header=FALSE)
snp_metadata <- read.delim("~/Biologi/Pop Gen/12 project/data/snps_filtered.inp", header=FALSE); names(snp_metadata) = c("name", "chromosome", "pos", "ancestral", "derived")
#snp = cbind(snp_metadata, genotypes_AF)





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

binary_snps = ancestral2zero(snp_metadata$ancestral, genotypes_AF)

# Parametriser
# Dependencies
#   snp_metadata is the metadata file
#   binary_snps is the output of ancestral2zero with a specific genotype file from the populations.

n_snps = dim(binary_snps)[1]
n_snps_per_window = 1000
n_folds = round(n_snps/n_snps_per_window)
e_results = tibble(pos=NA, LD=NA)[-1,]
if(n_snps > n_folds) {
print(paste("fold size: ", n_snps/n_folds))
    for (fold in 1:n_folds) { # finda way to skip when the SD = 0
        
        start = floor((fold-1)*(n_snps/n_folds)+1)
        end = floor((fold)*(n_snps/n_folds))
        #print(paste(fold, ":", start, end))
        
        correlation = cor(t(binary_snps[start:end,])) # transpose to get between-individuals corr.
        # Is there a nice way to catch NAs
        diag(correlation) = NA
        
        #print(correlation)
        res_mean =  mean(correlation, na.rm = T)
        print(res_mean)
        
        #e_results[fold] = res_mean
        
        e_results = rbind(e_results, tibble(pos = snp_metadata$pos[start + abs(start-end)/2], LD = res_mean))
        # Optionally, add the chromosomal position?
    }
}

ggplot(e_results) + 
    geom_line(aes(x=pos, y=LD), size = 0.3) +
    xlab("chromosome X position") + ylab(paste("mean LD for", n_folds, "sequential windows with", n_snps_per_window, "SNPs in each")) +
    ggtitle("Linkage disequilibrium: AF")


