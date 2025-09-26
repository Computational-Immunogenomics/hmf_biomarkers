import os
aqui = os.getcwd()

steps = ['0_clinical/', '1_database/', '2_biomarkers/']

for i in steps:
    print("Go! " + i)
    os.system('python prep.py ' + i)
