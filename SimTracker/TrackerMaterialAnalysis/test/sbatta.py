#!/usr/bin/env python

from ROOT import *
import re

ifile = open("dataRadLenEta.log", "r")
ofile = TFile("dataRadLenEta.root", "RECREATE")
histo = TProfile("RadLen", "RadLen", 250, -5., 5., 0., 10.)
for line in ifile:
  m = re.match('.*eta: (?P<ETA>-?\d+\.\d+).*(?P<RL>\d+\.\d+)', line)
  if m:
    histo.Fill(float(m.groupdict()['ETA']), float(m.groupdict()['RL']))
histo.Write()
ofile.Close()
