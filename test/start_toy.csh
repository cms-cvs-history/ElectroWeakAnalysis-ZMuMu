#########################
#                       # 
# author: Pasquale Noli #
# INFN Naples           #
# script to run ToyMC   #
#                       #
#########################

#!/bin/csh
rm fitResalt*
mkdir outputToy
set i=1
while ($i <= $1)
     	toyMonteCarlo -p analysis_Z_133pb_trackIso_3.root -n 1 -s $1 -T 0.998364  -S 0.989626  -I 0.980196  -H 0.915496 -y 50381.3
	mergeTFileServiceHistograms -o analysis_$i.root -i zmm_1.root bkg_1.root
	zChi2Fit analysis_$i.root >& log_fit_$i.txt
	rm -f analysis_$i.root
	mv *ps outputToy
	mv log_fit_$i.txt outputToy
	set i=`expr $i + 1`
	echo  $i
end
	root -l -q create_tree_for_toyMC.C
	root -l -q pulls.C
	mv pulls.eps outputToy
echo "Pulls are saved into pulls.eps"
    
gv outputToy/pulls.eps 
