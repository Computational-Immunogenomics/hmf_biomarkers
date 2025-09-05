jupyter nbconvert --to script *.ipynb
echo "Plot cohorts"; Rscript 0_cohorts.r
echo "Plot volcano"; Rscript 1_volcano.r
echo "Plot summaries"; Rscript 2_summaries.r
echo "Plot highlights"; Rscript 3_highlights.r
echo "Plot PFS"; Rscript 4_highlights_pfs.r
echo "Plot interaction"; Rscript 5_interaction.r
echo "Plot power"; Rscript 6_power.r
