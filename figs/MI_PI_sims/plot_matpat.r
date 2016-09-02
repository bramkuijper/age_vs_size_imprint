the.data <- read.table("summary_local_mating_mi_pi.csv",sep=";",header=T)

pdf("matpat.pdf")
print(xyplot(meanam + meanaf ~ dm | nfp * nmp * imprint * system,
                data=the.data,
                auto.key=T,
                strip=function(...,strip.levels) { 
                    strip.default(...,strip.levels=T) }
                ))
dev.off()

the.data <- read.table("summary_nonlocal.csv",sep=";",header=T)
pdf("mi_pi_nonlocal.pdf")
print(xyplot(meanam + meanaf ~ dm | l * imprint * system,
                data=the.data,
                auto.key=T,
                strip=function(...,strip.levels) { 
                    strip.default(...,strip.levels=T) }
                ))
dev.off()
