import os
aqui = os.getcwd()

steps = ['/0_curate_clinical', '/1_create_database', '/2_prep_dna', '/3_prep_rna']

for i in steps:
    print("Go! " + i)
    os.chdir(aqui + i)
    os.system('./run.sh')
