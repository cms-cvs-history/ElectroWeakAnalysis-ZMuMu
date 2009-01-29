#########################
#                       # 
# author: Pasquale Noli #
# INFN Naples           #
# script to run ToyMC   #
#                       #
#########################

#!/bin/csh
rm fitResalt*
set i=1
while ($i <= $1)
     	$CMSSW_BASE/test/slc4_ia32_gcc345/generator  -n 1 -s $1
	mergeTFileServiceHistograms -o analysis_$i.root -i zmm_0.root bkg_0.root
	set i=`expr $i + 1`
	echo  $i
end
	root -l -q create_tree_for_toyMC.C
echo "done."
