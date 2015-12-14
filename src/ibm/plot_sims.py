#!/usr/bin/env python

import sys, re, os.path, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# read in the csv file
dictdat = csv.DictReader(open(filename,"r"), delimiter=";")

# process the parameters at the end of the file
def process_params(dictionary, rowctr):

    fo = open(filename,"r")
    fl = fo.readlines()

    params = {};

    for line in fl[rowctr:]:
        if line.strip() != "":
            splitted = line.strip().split(";")
            params[splitted[0]] = splitted[1]

    return params



# data can now only accessed through looping row by row
# whereas we want lists of each column
# this function does that
def get_csvdict_by_column(the_raw_dict):

    # initialize a empty dict to contain the data
    by_column_dat = {}

    rowctr = 0

    # loop through the rows of the csv file and
    # put data in the dictionary
    for row in dictdat:

        if rowctr == 0:
            for key in row.keys():
                if key != "":
                    by_column_dat[key] = []

        for key, val in row.iteritems():

            if val == "system":
                params = process_params(dictdat,rowctr+1)

                return (params,by_column_dat)
                
            elif key != "" and val != None and val != "":
                by_column_dat[key].append(float(val))

        rowctr += 1

    return({},by_column_dat)
    

params,histdat = get_csvdict_by_column(dictdat)

# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 2

# add first subplot
plt.subplot(num_rows,1,1)
plt.plot(histdat["generation"],histdat["meanaf"],'r',
        histdat["generation"],histdat["meanam"],'b',linewidth=1)
#plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'age at maturity $a_{\mathrm{f}}$, $a_{\mathrm{m}}$')
plt.legend((r'$a_{\mathrm{f}}$',r'$a_{\mathrm{m}}$'))
plt.ylim(0,plt.ylim()[1])

plt.subplot(num_rows,1,2)
plt.plot(histdat["generation"],histdat["varaf"],'c',
        histdat["generation"],histdat["varam"],'m',linewidth=1)
#plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
plt.ylabel(r'variance $\sigma_{a_{\mathrm{f}}}^{2}$, $\sigma_{a_{\mathrm{m}}}^{2}$')
plt.legend((r'$\sigma_{a_{\mathrm{f}}}^{2}$',r'$\sigma_{a_{\mathrm{m}}}^{2}$'))

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
