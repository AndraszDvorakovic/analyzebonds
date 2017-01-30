######## INPUT FOR PLOTTING SCRIPT #######not for distanalyzer

Ntraj=3                   #number of trajectories
mainfolder='ION1TRAJ_0ST'   #folder with trajectories, same as in inanalyze.in

######### END OF INPUT  ######
analyzedfolder=$mainfolder'analyzed/'
i=1
j=1

cat >> ingnux.in << EOF
set terminal png size 800,600
set xlabel "Time (fs)"
set ylabel "Bond lenght (A)" 
EOF

while [ $i -le $Ntraj ]

do
inputfile1='./'$analyzedfolder$i'_bonds1to6.txt'
inputfile2='./'$analyzedfolder$i'_bonds7to12.txt'
inputfile3='./'$analyzedfolder$i'_bonds13to18.txt'
outputfile=$i'_dist.png'

cat >> ingnux.in << EOF
set output '$outputfile'
plot for [col=2:7] '$inputfile1' using 1:col with lines title columnheader, for [col=2:7] '$inputfile2' using 1:col with lines title columnheader,for [col=2:7] '$inputfile3' using 1:col with lines title columnheader 
EOF
let i++

done

gnuplot ingnux.in

cd $analyzedfolder
mkdir GRAPHS
cd ..
while [ $j -le $Ntraj ]

do
graph=$j'_dist.png'
cp $graph $analyzedfolder'GRAPHS'
rm $graph

let j++

done

rm ingnux.in
