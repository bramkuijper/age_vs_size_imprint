dat <- read.table("summary_matpat_new.csv",sep=";",header=T)

pdf("fig_mp.pdf")
print(xyplot(meanaf + meanam ~ dm | nmp * nfp * system * imprint,
                auto.key=T,
                strip=function(...,strip.levels) { strip.default(..., strip.levels=T) },
                data=dat
                ))
dev.off()

dat <- read.table("summary_matpat_nonlocal.csv",sep=";",header=T)

pdf("fig_mp_nonlocal.pdf")
print(xyplot(meanaf + meanam ~ dm | nmp * nfp * system * imprint * l,
                auto.key=T,
                strip=function(...,strip.levels) { strip.default(..., strip.levels=T) },
                data=dat
                ))
dev.off()
