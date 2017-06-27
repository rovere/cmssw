#! /usr/bin/env python

from ROOT import *
import sys, re

files = sys.argv[1:]


keepAlive = []
colors = [kRed, kOrange+1, kBlue+1, kGreen]

pixels = {}
for f in files:
  r = TFile(f)
  m = re.match('.*R000(\d+).*', f)
  if m:
    keepAlive.append(r)
    for layer in range(1,5):
      if m.group(1) not in pixels.keys():
        pixels[m.group(1)] = []
      pixels[m.group(1)].append(r.FindObjectAny("PixelClusterSizeY_%d" % layer))


for run in pixels.keys():
  l = TLegend(0.1,0.7,0.48,0.9)
  c = TCanvas("Smile_%s" % run, "Smile_%s" % run, 1024, 1024)
  for layer in range(0,4):
    p = pixels[run][layer].ProfileX("%s_%d" % (run, layer))
    keepAlive.append(p)
    p.SetLineColor(colors[layer])
    p.SetMarkerStyle(kFullCircle)
    p.SetMarkerColor(colors[layer])
    l.AddEntry(p, "PixelBarrelLayer%d" % layer)
    p.Draw("SAME")
    l.Draw()
  c.SaveAs("Run%d.png" % int(run))
print pixels
