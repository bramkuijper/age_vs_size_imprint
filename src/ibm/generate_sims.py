#!/usr/bin/env python

import numpy as np
import math, sys

#nf_nm = [[2,2],[3,2],[2,3]]
nf_nm = [[2,2]]
dm = list(np.arange(0,1.05,0.05))
q = 1
k = 1.0/3

type = [ 0, 1, 2, 3, 4 ]
diploid = [ 0, 1 ] 

nonlocalmating = [ 0.5 ] 

exe = "./xage_vs_size_imprint"

mu = 0.01
sdmu = 0.02

reps = 5

ctr = 0

for nfnm_i in nf_nm:
    nf_i = nfnm_i[0]
    nm_i = nfnm_i[1]
    for dm_i in dm:
        df_i = 1-dm_i
        for nonlocal_i in nonlocalmating:
            for type_i in type:
                for diploid_i in diploid:
                    for rep in range(0, reps):
                        print("echo " + str(ctr))
                        ctr += 1

                        exe_i = exe
                        exe_i += "_nf" + str(nf_i) + "_nm" + str(nm_i) + " "


                        print(exe_i + 
                                str(mu) + " " + str(sdmu) + " " + 
                                str(df_i) + " " + str(dm_i) + " " +
                                str(k) + " " + str(q) + " " + str(nonlocal_i) + " " +
                                str(diploid_i) + " " + str(type_i))

