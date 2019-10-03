library("rootSolve")

source("mathematica_expressions_r.txt")


get_sols <- function(nf, nm, l, k, μ, df, dm) {

    # expression patterns
    exp.pat <- c(
            "" # autosomal, offspring
            ,"M" # maternal
            ,"P" # paternal
            ,"PI" # padumnal
            ,"MI" # madumnal
            )

    # genetic system
    gen.sys <- c("dip","hapdip")

    # allocate data frame to store results
    results <- expand.grid(expression=exp.pat
            ,genetic_system=gen.sys
            ,af=NA
            ,am=NA
            ,nm=nm
            ,nf=nf
            ,l=l
            ,k=k
            ,mu=μ
            ,df=df
            ,dm=dm
            )

    params <- list(
            nf=nf
            ,nm=nm
            ,l=l
            ,k=k
            ,μ=μ
            ,df=df
            ,dm=dm
            ,E=exp(1.0))

    Power <- function(a,b) { return(a^b) }

    for (row_i in 1:nrow(results))
    {
        # specify selection gradient function
        dWdz <- function(x, params) {

            params <- c(list(af=x[1], am=x[2]), params)

            female_grad <- with(params, eval(
                            expr=eval(parse(text=
                                    paste("af"
                                            ,results[row_i,"genetic_system"]
                                            ,results[row_i,"expression"]
                                            ,sep="")
                                    ))
                            ))
            
            male_grad <- with(params, eval(
                            expr=eval(parse(text=
                                    paste("am"
                                            ,results[row_i,"genetic_system"]
                                            ,results[row_i,"expression"]
                                            ,sep="")
                                    ))
                            ))

            return(c(female_grad,male_grad))
        }

        vals <- multiroot(f = dWdz
                ,start=c(1.5,1.5)
                ,params=c(params))

        results[row_i,"af"] <- vals$root[[1]]
        results[row_i,"am"] <- vals$root[[2]]
    }

    return(results)
}

dfvals <- seq(0.01,0.99,0.01)

nfvals <- seq(1,7.05,0.05)

lvals <- c(0.1,0.5, 0.9, 1.0)

the.df <- NULL

ctr <- 1

for (l.i in lvals)
{
    for (df.i in dfvals)
    {
        for (nf.i in nfvals)
        {
            print(ctr)
            temp.df <- get_sols(nf=nf.i, nm=8-nf.i, l=l.i, k=1.0/3, μ=1.0, df=df.i, dm=1.0 - df.i)

            the.df <- rbind(the.df, temp.df)

            ctr <- ctr + 1
        }
    }
}

write.table(x=the.df,file="battleground_data_vary_df_nf.csv",sep=";",row.names=F)

#print("df done.")
#
#the.df <- NULL
#
#nfvals <- seq(1.0, 8.0, 0.05)
#
#lvals <- c(0.1,0.5, 0.9, 0.99)
#
#for (l.i in lvals)
#{
#    for (nf.i in nfvals)
#    {
#        print(nf.i)
#        print(max(nfvals) + 1 - nf.i)
#        temp.df <- get_sols(nf=nf.i, nm=max(nfvals) + 1 - nf.i, l=l.i, k=1.0/3, μ=1.0, df=0.5, dm=0.5)
#        
#        the.df <- rbind(the.df, temp.df)
#    }
#}
#
#write.table(x=the.df,file="battleground_data_vary_nf.csv",sep=";",row.names=F)
