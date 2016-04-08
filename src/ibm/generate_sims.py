#!/usr/bin/env python

import numpy as np
import math, sys

nf_nm = [[2,2],[3,2],[2,3]]
df = [ 0.1, 0.5, 0.9 ]
dm = list(np.arange(0,1.05,0.05))
q = 1
k = 1.0/3

type = [ 3, 4 ]
diploid = [ 0, 1 ] 

exe = "./xage_vs_size_imprint"
#exe = "/home/uccoaku/age_vs_size_imprint/src/ibm/xage_vs_size_imprint"

mu = 0.01
sdmu = 0.02

reps = 1

ctr = 0

for nfnm_i in nf_nm:
    nf_i = nfnm_i[0]
    nm_i = nfnm_i[1]
    for df_i in df:
        for dm_i in dm:
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
                                str(k) + " " + str(q) + " " +
                                str(diploid_i) + " " + str(type_i))

