# plot a contour plot that shows the level of polymorphism
# present in a numerical analysis that assumes N mitochondria
# per cell and compare this with a haploid population genetics
# model


library("lattice")
library("grid")
library("rootSolve")
library("VGAM")
source("/home/bram/R/src/bramlib.r")

# heights and widths of the panels that display the graphs
heights <- c(0.2, 1,0.2,1,0.3)
widths <- c(0.3,1,0.3,1,0.8)

# parameters for the graphics
line.lwd <- 0.3
lwd <- 0.3
tick.cex <- 0.4
label.cex <- 0.7
legend.cex <- 0.5
key.cex <- 0.7
plot.tck <- -0.4
plot.lwd <- 0.5
main.contour.level = 0.05

lim <- 1.0 

#file.name <- "leakage_contour_bottleneck"

E <- exp(1)

Power <- function(x, y) {
    return(x^y);
}

thetam <- 2
thetaf <- -2
sigmam <- 1.5
sigmaf <- 1.5


source("mathematica_expressions.txt")

# custom panel function to be called by xyplot
panel.opt <- function(x,y,
        cols, # list of line colors
        ltys, # list of line types
        lwds, # list of line widths
        params,
        gradient_expressions, # names of the gradient expressions to be plotted
                                # corresponding to the expression names in mathematica_expressions.r
        af, # two booleans saying we need to plot af and am or one of them
        deltas, # list of deltas, 
                #constants to be subtracted or added from a line so that lines don't overlap
        ...)
{
    dm_list <- seq(0,1.0,0.05)

    for (i in 1:length(gradient_expressions))
    {
        gradient_expression_i <- gradient_expressions[[i]]
        af.grad <- gradient_expression_i$af
        am.grad <- gradient_expression_i$am
        dWdz <- function(x, params)
        {
            params <- c(list(af=x[1],am=x[2]),params)

            aftplus1 <- with(params,eval(af.grad))
            amtplus1 <- with(params,eval(am.grad))

            return(c(aftplus1,amtplus1))
        }

        vals.af <- c()
        vals.am <- c()

        # multivariate root solving
        for (dm.i in dm_list)
        {
            vals <- multiroot(f = dWdz, start=c(0.5,0.5), params=c(params,list(dm=dm.i,df=1-dm.i)))

            vals.af <- c(vals.af,vals$root[[1]])
            vals.am <- c(vals.am,vals$root[[2]])
        }

        y.vals <- vals.af

        if (!af)
        {
            y.vals <- vals.am
        }
        # plot af
        panel.xyplot(x=dm_list,y=y.vals+deltas[[i]],
                        col=cols[[i]],
                        lty=ltys[[i]],
                        lwd=lwds[[i]],
                        type="l"
                        )
    }
}


block <- function(row, col, xlab="", ylab="",
                    print.xlabels=T,
                    print.ylabels=T,
                    ind.label="A",
                    label="",
                    cols,
                    lwds,
                    ltys,
                    ylim=c(-0.1,4.5),
                    af,
                    params,
                    gradient_expressions,
                    deltas=list(c(0,0,0,0)),
                    sub.label="")
{
    xp <- xyplot(seq(0,1,0.01) ~ seq(0,1.0,0.01),
                    xlim=c(-.05,1.05),
                    ylim=ylim,
                    panel=panel.opt,
                    cols=cols,
                    gradient_expressions=gradient_expressions,
                    lwds=lwds,
                    af=af,
                    ltys=ltys,
                    params=params,
                    deltas=deltas
                    )
            
    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        do.call("panel.opt",trellis.panelArgs(xp,1))

        #        grid.rect(gp=gpar(lwd=lwd,fill="transparent"))
        grid.lines(x=c(0,1),y=c(0,0),gp=gpar(lwd=lwd))
        grid.lines(x=c(0,0),y=c(0,1),gp=gpar(lwd=lwd))

        grid.text(x=0.5,y=1.1,just="centre",label=label,gp=gpar(cex=legend.cex))
        grid.text(x=0.83,y=0.27,just="centre",label=sub.label,gp=gpar(cex=label.cex))
        grid.text(x=0.05,y=0.95,just="left",label=ind.label,gp=gpar(cex=label.cex))

    upViewport()
    
    pushViewport(viewport(layout.pos.row=row+1,
                            layout.pos.col=col,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

        single.axis(range=xp$x.limits,
                        side="top",
                        labels=print.xlabels,
                        cex=tick.cex,
                        lwd=line.lwd,
                        labelcex=label.cex,
                        tck=plot.tck,
                        distance=0.5,
                        y.text.off=0.5,
                        nsub=5,
                        text=ifelse(row==length(heights),xlab,"")
                        )
    upViewport()
   
   # x-axis 
    pushViewport(viewport(layout.pos.row=row,
                            layout.pos.col=col-1,
                            xscale=xp$x.limits,
                            yscale=xp$y.limits
                            ))

            expr <- ylab
        single.axis(range=xp$y.limits,
                        side="right",
                        tck=plot.tck,
                        labels=print.ylabels,
                        lwd=line.lwd,
                        cex=tick.cex,
                        x.text.off=0.3,
                        labelcex=label.cex,
                        distance=0.5,
                        nsub=5,
                        text=ifelse(col==2,expr,""))
    upViewport()
}

file.name <- paste("philopatry",sep="")

all.expressions <- 
    list(
    dip=list(
                    list(af=afdip,am=amdip),
                    #                    list(af=afdipM,am=amdipM),
                    #                    list(af=afdipP,am=amdipP),
                    list(af=afdipMI,am=amdipMI),
                    list(af=afdipPI,am=amdipPI)
            ),
    hap=list(
                    list(af=afhapdip,am=amhapdip),
                    #                    list(af=afhapdipM,am=amhapdipM),
                    #                    list(af=afhapdipP,am=amhapdipP),
                    list(af=afhapdipMI,am=amhapdipMI),
                    list(af=afhapdipPI,am=amhapdipPI)
    )
)

print(names(all.expressions))


for (expression_name in names(all.expressions))
{
init.plot(filename=paste(file.name,expression_name,sep=""),width=550,background="transparent",height=500,font="helvetica",type="pdf")

lo <- grid.layout(
                    ncol=length(widths),
                    nrow=length(heights),
                    heights=heights,
                    widths=widths)

expression_set_i <- all.expressions[[expression_name]]

pushViewport(viewport(layout=lo))

    lwds_general=c(rep(0.5,times=4),.75,.75)
    # use the block function to generate one graph
    block(row=2,col=2,
            print.xlabels=T,
            ind.label=expression(paste("A, local male mating only, ",italic(ℓ),"=1")),
            ylim=c(-0.1,6),        
            params=list(l=1.0, nf=2, nm=2, μ=1.0, k=1.0/3),
            gradient_expressions=expression_set_i,
            af=T,
            cols=c("red","blue","darkgreen","orange","purple"),
            ltys=list(1,"22","11","22",1,1),
            lwds=c(rep(0.75,times=5)),
            deltas=c(0,0,0,0,0,0)
            )
    
    block(row=4,col=2,
            print.xlabels=T,
            ind.label=expression(paste("C, nonlocal male mating, ",italic(ℓ),"=0.5")),
            ylim=c(-0.1,6),        
            params=list(l=0.5, nf=2, nm=2, μ=1.0, k=1.0/3),
            gradient_expressions=expression_set_i,
            af=T,
            cols=c("red","blue","darkgreen","orange","purple"),
            ltys=list(1,"22","11","22",1,1),
            lwds=c(rep(0.75,times=5)),
            deltas=c(0,0,0,0,0,0)
            )

    block(row=2,col=4,
            print.xlabels=T,
            ind.label=expression(paste("B, local male mating only, ",italic(ℓ),"=1")),
            ylim=c(-0.1,3),        
            params=list(l=1.0, nf=2, nm=2, μ=1.0, k=1.0/3),
            af=F,
            gradient_expressions=expression_set_i,
            cols=c("red","blue","darkgreen","orange","purple"),
            ltys=list(1,"22","11","22",1,1),
            lwds=c(rep(0.75,times=5)),
            deltas=c(0,-0.01,0.01,-0.02,-0.02,0)
            )

    block(row=4,col=4,
            print.xlabels=T,
            ind.label=expression(paste("D, nonlocal male mating, ",italic(ℓ),"=0.5")),
            ylim=c(-0.1,3),        
            params=list(l=0.5, nf=2, nm=2, μ=1.0, k=1.0/3),
            af=F,
            gradient_expressions=expression_set_i,
            cols=c("red","blue","darkgreen","orange","purple"),
            ltys=list(1,"22","11","22",1,1),
            lwds=c(rep(0.75,times=5)),
            deltas=c(0,-0.01,0.01,-0.02,-0.02,0)
            )
    # axes labels
    grid.text(x=0.43,y=0.06,
                label=expression(paste("Male dispersal probability, ",italic(d)["m"]," = ",1-{},"  ",italic(d)["f"])),
                gp=gpar(cex=label.cex))

    grid.text(x=0.04,y=0.5,
                rot=90,
                label=expression(paste("Female age at maturity, ",italic(a)["f"])),
                gp=gpar(cex=label.cex))

    grid.text(x=0.43,y=0.5,
                rot=90,
                label=expression(paste("Male age at maturity, ",italic(a)["m"])),
                gp=gpar(cex=label.cex))

    # legend
    pushViewport(viewport(layout.pos.row=2,
                            layout.pos.col=5))

        y.pos <- 0.98
        grid.text(x=0.1,
                    y=y.pos,
                    just="left",
                    label="",
                    gp=gpar(cex=key.cex))

        y.pos <- 0.9
        grid.lines(x=c(0.1,0.2),
                    y=c(y.pos,y.pos),
                    gp=gpar(col="red",
                            lwd=0.75,
                            lineend="square",
                            lty=1))
        grid.text(x=0.25,
                    y=y.pos,
                    just="left",
                    label="Offspring expression",
                    gp=gpar(cex=key.cex))


        y.pos <- 0.82
        grid.lines(x=c(0.1,0.2),
                    y=c(y.pos,y.pos),
                    gp=gpar(col="blue",
                            lwd=0.75,
                            lineend="square",
                            lty="22"))
        grid.text(x=0.25,
                    y=y.pos,
                    just="left",
                    label="Expression MI allele",
                    gp=gpar(cex=key.cex))

        
        y.pos <- 0.74
        grid.lines(x=c(0.1,0.2),
                    y=c(y.pos,y.pos),
                    gp=gpar(col="darkgreen",
                            lwd=0.75,
                            lineend="square",
                            lty="11"))
        grid.text(x=0.25,
                    y=y.pos,
                    just="left",
                    label="Expression PI allele",
                    gp=gpar(cex=key.cex))


            #        y.pos <- 0.66
            #        grid.lines(x=c(0.1,0.2),
            #                    y=c(y.pos,y.pos),
            #                    gp=gpar(col="orange",
            #                            lwd=0.75,
            #                            lineend="square",
            #                            lty="22"))
            #        grid.text(x=0.25,
            #                    y=y.pos,
            #                    just="left",
            #                    label="Expression MI-allele",
            #                    gp=gpar(cex=key.cex))
            #
            #
            #        y.pos <- 0.58
            #        grid.lines(x=c(0.1,0.2),
            #                    y=c(y.pos,y.pos),
            #                    gp=gpar(col="purple",
            #                            lwd=0.75,
            #                            lineend="square",
            #                            lty=1))
            #        grid.text(x=0.25,
            #                    y=y.pos,
            #                    just="left",
            #                    label="Expression PI-allele",
            #                    gp=gpar(cex=key.cex))
            #

    upViewport()
upViewport()

dev.off()

}


