library("lattice")
library("magrittr")
library("tidyr")
library("dplyr")
library("readr")
library("colorRamps")
library("RColorBrewer")

# plot conflicts out in a contour plot

# read in data and create subset
if (!exists("the.data.long"))
{
    the.data <- read.table("battleground_data_vary_df_nf.csv"
            ,sep=";"
            ,header=T)

    the.data <- the.data %>% subset(l %in% c(1.0,0.5))

    # change into long form
    the.data.long <- the.data %>% pivot_wider(
            names_from = c("expression")
            ,values_from = c("af","am"))
}

#at.range <- c(-100,seq(-0.5,0.5,0.02),100)

#at.range <- c(seq(0,0.5,0.005),100)

at.range <- c(0.05, 0.2, 0.5)

#col.reg <- colorRampPalette(brewer.pal(8,"Dark2"))(length(table.data))
#rwb <- colorRampPalette(colors = c("red", "white", "blue"))

col.reg <- c("white","grey80","grey70","grey50","grey40")

the.strip <- function(strip.levels,...) { strip.default(strip.levels=T,...) }

ylim <- c(1,6.9)

pdf("levelplot_battleground_AF_mi_vs_pi.pdf")
print(
        contourplot(
                abs(af_PI - af_MI) ~ df * nf | genetic_system * l
                ,data=the.data.long
                ,at=at.range
                ,strip=the.strip
                ,ylim=ylim
                ,label=F
                ,col.regions=col.reg
                )
        )
dev.off()

pdf("levelplot_battleground_AM_mi_vs_pi.pdf")
print(
        contourplot(
                abs(am_PI - am_MI) ~ df * nf | genetic_system * l
                ,data=the.data.long
                ,at=at.range
                ,strip=the.strip
                ,ylim=ylim
                ,label=F
                ,col.regions=col.reg
                )
        )
dev.off()
                
pdf("levelplot_battleground_AF_mother_offspring.pdf")
print(
        levelplot(
                abs(af_ - af_M) ~ df * nf | genetic_system * l
                ,data=the.data.long
                ,at=at.range
                ,strip=the.strip
                ,ylim=ylim
                ,col.regions=col.reg
                )
        )
dev.off()
