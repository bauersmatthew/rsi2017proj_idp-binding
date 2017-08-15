#!/usr/bin/python
import subprocess
import sys

prots = []
with open(sys.argv[1]) as fin:
    for line in fin:
        line = line.rstrip()
        if line:
            prots.append(line.split())

for prot in prots:
    subprocess.call('wget https://files.rcsb.org/download/{}.pdb'.format(prot[5]), shell=True)
    subprocess.call('mv {}.pdb {}/{}.pdb'.format(prot[5], sys.argv[2], prot[0]), shell=True)


