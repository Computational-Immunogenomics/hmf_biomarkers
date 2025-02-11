
from openai import OpenAI
import os 
import pandas as pd
import pickle

### Use ChatGPT AI key ### 
client = OpenAI(
  api_key="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"
)

### ChatGPT whisperers ### 

def format_question( genes_string: str):
	question = "Can you provide a annotated table for the following genes (" + genes_string + ") with columns that describe wheter they are cancer driver genes, oncogene or tumor suppressor, cancer pathway, protein type, sub-cellular location, and functional role in the cell?"
	return question 

def ask_for_annotation( genes_string: str ):
	question = format_question(genes_string)	
	answer = client.chat.completions.create(model="gpt-4", store=True, messages=[{"role": "user", "content": question}])
	return answer

#### Run a marathon #### 
annotation = pd.read_csv("/Users/josephusset/Desktop/work/hartwig_driver_annotation_chatGPT.csv", on_bad_lines='warn', sep = ";")
#annotation = pd.read_csv("all_genes.csv", on_bad_lines='warn', sep = ";")
genes = [i for i in annotation['Gene']][1:20]

answers = []
group_size = 20
for i in range(0, len(genes), group_size):
  print(i)
  genes_string = ", ".join(genes[i:min(i + group_size, len(genes))])
  tell_me = ask_for_annotation(genes_string)
  answers.append(tell_me)

with open("go_annotation_raw.pkl", "wb") as file:
    pickle.dump(answers, file)
