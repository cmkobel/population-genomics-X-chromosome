source("projX functions.r")

#   C
# ------
# Perform an XP-EHH scan of the whole X chromosome for at least three populations. 
# Identify the 10 most significant regions and associated with genes as in A.







#   XP-EHH 
# ----------
# (pairwise population tests) 

wg_xpehh_AF_WE = ies2xpehh(res_scan_AM, res_scan_CAS, popname1 = "AF", popname2 = "WE", method = "bilateral")
wg_xpehh_WE_EA = ies2xpehh(res_scan_AM, res_scan_CAS, popname1 = "WE", popname2 = "EA", method = "bilateral")
wg_xpehh_EA_AF = ies2xpehh(res_scan_AM, res_scan_CAS, popname1 = "EA", popname2 = "AF", method = "bilateral")


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
pdf("../plots/b_ehh/xpehh_AF_WE.pdf"); plot_xpehh(wg_xpehh_AF_WE, "AF_WE"); dev.off()
pdf("../plots/b_ehh/xpehh_WE_EA.pdf"); plot_xpehh(wg_xpehh_WE_EA, "WE_EA"); dev.off()
pdf("../plots/b_ehh/xpehh_EA_AF.pdf"); plot_xpehh(wg_xpehh_EA_AF, "EA_AF"); dev.off()

head(wg_xpehh_AF_WE[3])

