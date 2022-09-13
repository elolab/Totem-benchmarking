WRKDIR="/wrk/local2/B22007_nonlinear_scShaper/B20006_Trajectory"
RESULTDIR=${WRKDIR}"/comparison/Totem_123456_slingshot"
mkdir $RESULTDIR
cd $WRKDIR

BENCHMARKING_FILES=$(find data/linear_tree -type f -name "*rds")

R_SCRIPT="run_dynverse_benchmarking_Totem_slingshot_distance_with_args.R"

for BENCHMARKING_FILE in $BENCHMARKING_FILES
do
	echo $RESULTDIR
	echo $BENCHMARKING_FILE
	sbatch -J R --mem 50000 --partition long  --wrap="module add R/4.1.0 ; Rscript --slave $R_SCRIPT $BENCHMARKING_FILE $RESULTDIR"
	#exit 0

done
