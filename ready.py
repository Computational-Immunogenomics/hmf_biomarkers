import os
aqui = os.getcwd()

steps = ['0_clinical/', '1_database/', '2_biomarkers/', '3_signals/0_analysis/', '3_signals/1_figures/']

for i in steps:
    print("Go! " + i)
    os.system('python prep.py ' + i)
