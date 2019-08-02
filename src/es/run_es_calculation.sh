VERSION_STAMP=190412
for DATA_PREFIX in tabula_muris mousebrain_all campbell_lvl1 campbell_lvl2
do 
	# echo $DATA_PREFIX $VERSION_STAMP &> es_calculation.log.$DATA_PREFIX.txt & 
	time Rscript es_calculation.R $DATA_PREFIX $VERSION_STAMP &> es_calculation.log.$DATA_PREFIX.txt &
done

wait;
echo "Done bash script"

### ONE-LINER [Does not work]
# VERSION_STAMP=190412
# for datprefix in tabula_muris mousebrain campbell_lvl1 campbell_lvl2; do time Rscript es_calculation.R $datprefix $VERSION_STAMP &> es_calculation.log.$datprefix.txt &; done
