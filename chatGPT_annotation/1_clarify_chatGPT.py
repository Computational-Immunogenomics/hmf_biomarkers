

import pandas as pd
import pickle

# Open the pickle file in binary read mode
with open('go_annotation_raw.pkl', 'rb') as file:
    data = pickle.load(file)