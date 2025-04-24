import pandas as pd
import pickle
import re

with open('0_perplexity_prior_knowledge.pkl', 'rb') as file:
    data = pickle.load(file)

extract = {}
for cohort in data.keys():
	extract[cohort] = [i for i in data[cohort].to_dict()['choices'][0]['message']['content'].split("\n")[4:] if "|" in i]

citations = {i: data[i].to_dict()['citations'] for i in data.keys()}

df['mapped_urls'] = df.apply(lambda row: map_keys_to_urls(row['key_str_col'], row['dict_col']), axis=1)


def go(cohortGo: str = "Anti-AR ## Prostate", column_names: list = ["Biomarkers", "Gene", "Event", "Association", "Measurement", "Evidence", "Protein", "Function", "Citation", "Database"]):
	base = [[i.strip() for i in line.split("|")[1:-1]] for line in extract[cohortGo]]
	df = pd.DataFrame(base, columns=column_names)
	df.insert(0, 'CohortGo', cohortGo)
	df.insert(1, 'Cohort', cohortGo.split(" ## ")[0])
	df.insert(2, 'Treatment_Mechanism', cohortGo.split(" ## ")[1])
	return df

share = pd.concat([go(i) for i in extract.keys()], ignore_index=True)

### Add the citations
citations = {i: data[i].to_dict()['citations'] for i in data.keys()}

url_dicts = {
    cohort: {f"[{i+1}]": url for i, url in enumerate(urls)}
    for cohort, urls in citations.items()
}

def map_keys_to_urls( key_string, cohortGo):
	mapping_dict = url_dicts[cohortGo]
	keys = re.findall(r'\[(\d+)\]', key_string)
	urls = [mapping_dict.get(f'[{key}]') for key in keys if f'[{key}]' in mapping_dict]
	return ', '.join(urls)

share['Citations'] = share.apply(
    lambda row: map_keys_to_urls(row['Citation'], row['CohortGo']) if pd.notna(row['Citation']) and pd.notna(row['CohortGo']) else None,
    axis=1
)

share.to_csv("1_process.csv", index = False)
