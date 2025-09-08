jupyter nbconvert --to script *.ipynb
echo "prep features"; Rscript 0_prep_features.r;
echo "prep cohorts"; Rscript 1_prep_cohorts.r;
echo "run tests"; Rscript 2_run.r;
echo "prep results"; Rscript 3_ready.r;