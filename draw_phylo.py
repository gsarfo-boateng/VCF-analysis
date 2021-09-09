#!/usr/bin/env python3

from Bio import Phylo
import sys
import os
import matplotlib
matplotlib.use('Agg')
import pylab

import networkx

file_=sys.argv[1]

tree=Phylo.read(file_,'newick')
Phylo.draw(tree,do_show=False)
#plt.savefig('Phylo.jpg')

workdir=os.getcwd()

pylab.savefig('Phylo.jpg',dpi=1500)
pylab.savefig('Phylo.svg',dpi=1500)


print("Phylogenetic tree generated and saved as Phylo.jpg and Phylo.svg")
print("Phylogenetic tree generated and saved in directory:%s"%workdir)
