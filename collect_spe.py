from glob import glob
import os
import sys

import re

path = r'/home/inbaloz/Projects/NT_on_surfaces/interlayer_potential_Itai'

if len(sys.argv) > 1:
	path += sys.argv[1]



g = glob(path+"/SPE*dat")
g.sort(key=lambda x: int(x.split(".")[-2]))
print g
l = [open(i,'rb').read() for i in g]

#ris = [re.findall("E_KC= (.*) RI= (.*)", i.splitlines()[1])[0] for i in l]
ris = [re.findall("ToT_Energy= (.*)", i.splitlines()[0])[0] for i in l]
print "\n".join(ris)


