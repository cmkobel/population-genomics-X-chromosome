source("projX functions.r")

#   C
# ------


#   I:
# -----
# Perform an XP-EHH scan of the whole X chromosome for at least three populations. 


for(i in (c("AF","WE", "EA"))) {
    
    print(
        load(paste("data/res_scan_", i, ".rdata", sep=""), verbose=T)
    )
}





#   XP-EHH 
# ----------
# (pairwise population tests) 

wg_xpehh_AF_WE = ies2xpehh(res_scan_AF, res_scan_WE, popname1 = "AF", popname2 = "WE", method = "bilateral")
wg_xpehh_WE_EA = ies2xpehh(res_scan_WE, res_scan_EA, popname1 = "WE", popname2 = "EA", method = "bilateral")
wg_xpehh_EA_AF = ies2xpehh(res_scan_EA, res_scan_AF, popname1 = "EA", popname2 = "AF", method = "bilateral")


# Skal det her egentlig stå på den anden side af opgave C?

plot_xpehh = function(plot_df, title) {
    print(
        grid.arrange(
            ggplot(plot_df) + 
                geom_point(aes(x=POSITION, y=plot_df[3]), size=0.03) +
                xlab("") +
                ylab("XPEHH") + ggtitle(title),
            ggplot(plot_df) + 
                geom_point(aes(x=POSITION, y=`-log10(p-value) [bilateral]`), size=0.03) +
                xlab("chromosome X position") +
                ylab("-log10( p-value )"),
            layout_matrix = rbind(c(1),c(2))
        )
    )
}


# Africa and Europe, 
# Europe and East Asia
# East Asia and Africa 
pdf("plots/c_xpehh/xpehh_AF_WE.pdf"); plot_xpehh(wg_xpehh_AF_WE, "AF_WE"); dev.off()
pdf("plots/c_xpehh/xpehh_WE_EA.pdf"); plot_xpehh(wg_xpehh_WE_EA, "WE_EA"); dev.off()
pdf("plots/c_xpehh/xpehh_EA_AF.pdf"); plot_xpehh(wg_xpehh_EA_AF, "EA_AF"); dev.off()

head(wg_xpehh_AF_WE[3])


#   II:
# -------
# Identify the 10 most significant regions and associated with genes as in A.


overlap = function(candidates, gene_annotation) {
    setkey(candidates, start, end)
    return(
        foverlaps(gene_annotation, candidates, type="any", nomatch = 0)
    )
}


gtf <- as.data.table(rtracklayer::import("data/gencode.v17.annotation.gtf"))
gtf <- gtf[gtf$seqnames == 'chrX']   # to select only X chromosome


af_we_to_overlap = tibble(start = wg_xpehh_AF_WE$POSITION,
                          end = start,
                          xpehh = wg_xpehh_AF_WE[,3],
                          ppval = wg_xpehh_AF_WE$`-log10(p-value) [bilateral]`) %>% 
    na.omit()

we_ea_to_overlap = tibble(start = wg_xpehh_WE_EA$POSITION,
                          end = start,
                          xpehh = wg_xpehh_WE_EA[,3],
                          ppval = wg_xpehh_WE_EA$`-log10(p-value) [bilateral]`) %>% 
    na.omit()

ea_af_to_overlap = tibble(start = wg_xpehh_EA_AF$POSITION,
                          end = start,
                          xpehh = wg_xpehh_EA_AF[,3],
                          ppval = wg_xpehh_EA_AF$`-log10(p-value) [bilateral]`) %>% 
    na.omit()




get_peaks = function(to_overlap, value_column, percentile) { # value column skal nok være et tal den her gang
    # 1 c(start, end, value)
    # 2 c(value)
    # 3 percentile
    threshold = quantile(unlist(to_overlap[,4]), percentile)
    new_overlap_region = overlap(as.data.table(to_overlap[to_overlap[,4] >= threshold,]), gtf) # fst_af_we_100 burde kunne bruges lige så vel som reg_all_peak_cand.
    #return(new_overlap_region)
    #View(new_overlap_region)
    #return(unique(cbind(new_overlap_region$gene_name)))
    
    peak_result = new_overlap_region %>%
        group_by(gene_name) %>%
        summarise(position = start[which.max(ppval)],
                  xpehh = xpehh[which.max(ppval)],
                  ppval = ppval[which.max(ppval)],
                  transcript_type = gene_type[which.max(ppval)],
                  comment = "")
    
    
    
    return(peak_result)
    
}

#library(openxlsx)

get_peaks_af_we = get_peaks(af_we_to_overlap, 4, 0.99927)
get_peaks_af_we = get_peaks_af_we[order(get_peaks_af_we$ppval, decreasing = T),]
write.xlsx(get_peaks_af_we, "frames_out/xpehh_af_we_peaks_dataframe.xlsx")
#View(get_peaks_af)

get_peaks_we_ea = get_peaks(we_ea_to_overlap, 4, 0.9994)
get_peaks_we_ea = get_peaks_we_ea[order(get_peaks_we_ea$ppval, decreasing = T),]
write.xlsx(get_peaks_we_ea, "frames_out/xpehh_we_ea_peaks_dataframe.xlsx")
#View(get_peaks_af)

get_peaks_ea_af = get_peaks(ea_af_to_overlap, 4, 0.99829)
get_peaks_ea_af = get_peaks_ea_af[order(get_peaks_ea_af$ppval, decreasing = T),]
write.xlsx(get_peaks_ea_af, "frames_out/xpehh_ea_af_peaks_dataframe.xlsx")
#View(get_peaks_af)


