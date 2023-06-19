"""
phylo.py

Author: Jackson Sheppard
Last Edit: 06/19/23

Script to read aligned MAPT/MAP2/MAP4 Protein sequences for reproduction of
Bayesian consensus phylogenetic tree in Sundermann et al. BMC Genomics (2016)
17:264 Fig. 2a.
"""

import nexusformat.nexus as nex
import os

fname = os.path.join("TreeBASE", "1_1458651913_MAP-102x1953_ExaBayesConsensusEMR.nexorg")
a = nex.nxload(fname)
print(a.tree)