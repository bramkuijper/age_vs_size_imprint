import os, re

filename = "mathematica_expressions.r"

f = open(filename)
fl = f.readlines()
f.close()


regexps = [(μ,mu)]

for line in fl:
    if 
