from openai import OpenAI
import pandas as pd
import pickle
import sys

### Enter the API key to use perplexity client ### 
if len(sys.argv) > 1:
    api_key_go = sys.argv[1]
    print("Good to go!")
else:
    print("Please provide a perplexity API key.")

perplexity_client = OpenAI(api_key= api_key_go, base_url="https://api.perplexity.ai")

### ChatGPT whisperers ### 
def format_question( cohortGo: str):
	therapy = cohortGo.split(" ## ")[1]
	cohort = cohortGo.split(" ## ")[0]
	question = """
	"Can you generate a table of 15 known biomarkers linked to response/non-response for (""" + therapy + """) anti-cancer therapies in (""" + cohort + """) cohorts?  
    Please include the following columns for each biomarker:
	1. **Biomarker**: Name the biomarker.
	2. **Gene**: Gene name of affected biomarker.
	3. **Driver Event**:  (Amplification, Deletions, Mutation, Disruption). 
	4. **Association**: Indicate if the biomarker is linked to response or non-response.
	5. **How is biomarker measured?**: Describe how the biomarker is measured (e.g., RNA-seq, IHC, PCR, Somatic Mutations, Driver event).  
	6. **Evidence level**: Chose between Low, Medium, and High. 
	7. **Protein Type**: Specify the type of protein encoded (e.g., enzyme, transcription factor).
	8. **Functional Role**: Briefly describe the biomarker’s role in affecting treatment response.
	9. **Key Citation**: Provide a main literature reference.
	10. **Databases**: List databases used for annotation (e.g., TCGA, GEO, COSMIC).
	Please ensure that each biomarkers data is presented in a separate row of the table, with the appropriate details in each column."
	"""
	return question

def ask_perplexity(cohortGo: str = "Anti-AR ## Prostate"):
    question = format_question(cohortGo)
    message = [{"role": "user", "content": question}]
    return perplexity_client.chat.completions.create(model="sonar-pro",messages = message, temperature=0.2, top_p = .3)

### Run it 
annotation = pd.read_csv("/mnt/petasan_immunocomp/datasets/hartwig/biomarkers/share/top_mechanisms.csv", on_bad_lines='warn', sep = ",")
treatment_cohorts = [i for i in annotation['cohortGo'] if i != "Pan-Cancer"]
print(treatment_cohorts)

answers = {}
for i in treatment_cohorts:
    if i != "Pan-Cancer":
	    print(i)
	    answers[i] = ask_perplexity(i)

with open('0_perplexity_prior_knowledge.pkl', 'wb') as f:
    pickle.dump(answers, f)
