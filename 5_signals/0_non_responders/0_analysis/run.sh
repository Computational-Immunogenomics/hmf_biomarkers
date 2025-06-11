jupyter nbconvert --to script *.ipynb
echo "Prep features"; Rscript 0_prep_features.r;
echo "Prep cohorts"; Rscript 1_prep_cohorts.r;
echo "Run marginal tests"; Rscript 2_run_marginal.r;
echo "Add combination tests"; Rscript 3_run_interaction.r;
echo "Run simulation .."; Rscript 5_simulation.r;
