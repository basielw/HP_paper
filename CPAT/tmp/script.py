#!/usr/bin/python


import string
from optparse import OptionParser
import warnings
import string
import collections
import sets
import signal
from numpy import mean,median,std,nansum
from string import maketrans
import subprocess

#import 3rd party modules
#import pysam
#from bx.bbi.bigwig_file import BigWigFile
import numpy as np

#import my own modules
from cpmodule 	import fickett
from cpmodule  import orf
from cpmodule  import fasta
from cpmodule  import annoGene
from cpmodule  import FrameKmer
from cpmodule  import ireader

sequence = raw_input("Paste sequence here: ")
tmp = orf.ORFFinder(sequence)

print
print "When using the original:"
print tmp.longest_orf('+')
print
print "When using the modified -LONGEST:"
print tmp.all_orfs("+",['ATG'],['TAG','TAA','TGA'],25,"LONGEST")
print
print "When using the modified -ALL:"
print tmp.all_orfs("+",['ATG'],['TAG','TAA','TGA'],25,"ALL")
print