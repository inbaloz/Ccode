from glob import glob
import os
import re

path = r'/home/inbaloz/Projects/NT_on_surfaces/KC_RI_Itai/'
g = glob(path+"SPE*dat")
g.sort(key=lambda x: int(x.split(".")[-2]))
print g
l = [open(i,'rb').read() for i in g]

ris = [re.findall("E_KC= (.*) RI= (.*)", i.splitlines()[0])[0] for i in l]
print "\n".join("\t".join(x) for x in ris)


