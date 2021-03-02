time root -l > genLstandNoBKG.log 2>&1 <<EOF
.x ~/Scaricati/load.C
.L ~/Scaricati/na60plus/na60fastsim/files/macros/runLambdaCandidatesmio.C+
GenerateD0SignalCandidates(1000000,160.,"~/Scaricati/na60plus/na60fastsim/files/setup-5um-maps_Eff1.txt","/home/carolinda/Scaricati/PtY_spectra.root","/home/carolinda/Scaricati/na60plus/na60fastsim/files/decaytables/USERTABLAM.DEC")
.q
EOF
