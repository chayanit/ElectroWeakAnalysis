#!/bin/tcsh

root.exe -l <<EOF

TChain ch("invHiggsInfo/InvHiggsInfo")
ch.Add("invHiggsInfo_*.root")
TFile *file = TFile::Open("NEW.root","RECREATE")
file->mkdir("invHiggsInfo")->cd()
ch.Merge(file,10001)
EOF
