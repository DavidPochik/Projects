import numpy as np
import math
import sys
pythonpath=str(sys.argv[2])
sys.path.insert(0, pythonpath)
import athena_read

filename = 'GM_value.txt'
M_PNS    = float(sys.argv[1])   # Msun
G        = 6.6743e-8            # cgs
GM       = G * M_PNS * 1.989e33 # cgs
c        = [GM]
with open(filename,"w") as file:
	for x in zip(c):
		file.write(str(*x))
