#%%

import numpy as np 
import pandas as pd 
import tqdm
data = pd.read_csv('../analysis/data_analysisresults/secretion_microbiomestudies.csv')
studies = pd.read_csv('../analysis/data_curated_microbiome/relabundance_colData.csv')
valid_samples = pd.read_csv('./samples_forinteractive.csv')
data = data[data['Unnamed: 0'].isin(valid_samples['sample'].values)]

families = ['Bacteroidaceae',
            'Bifidobacteriaceae',
            'Coriobacteriaceae',
            'Enterobacteriaceae',
            'Lachnospiraceae',
            'Oscillospiraceae',
            'Prevotellaceae',
            'Tannerellaceae']

fermentation_products = ['glucose', 'maltose', 'acetate', 'butyrate', 'formate', 'lactate',
                         'propionate', 'succinate']

#%%
# %% Composition data
comp_df = pd.DataFrame([])
ferm_df = pd.DataFrame([])
for i in tqdm.tqdm(range(len(data))):
    d = data.iloc[i]
    g = d['Unnamed: 0']
    study_info = studies[studies['Unnamed: 0']==g]
    _comp_df = pd.DataFrame([])     
    for i, f in enumerate(families): 
        __comp_df = pd.DataFrame({'study_name': study_info['study_name'].values[0],
               'subject_id': study_info['subject_id'].values[0],
               'age_category': study_info['age_category'].values[0],
               'disease_state': study_info['disease'].values[0],
               'bacterial_family': f,
               'bm_fraction':d[f'familylevel_{f}']/100}, index=[0])
        _comp_df = pd.concat([_comp_df, __comp_df])
    tot = _comp_df['bm_fraction'].unique().sum()
    __comp_df = pd.DataFrame({'study_name': study_info['study_name'].values[0],
                         'subject_id': study_info['subject_id'].values[0],
                         'age_category': study_info['age_category'].values[0],
                         'disease_state': study_info['disease'].values[0],
                         'bacterial_family': 'not characterized', 
                         'bm_fraction':1-tot}, index=[0])
    comp_df = pd.concat([comp_df, _comp_df, __comp_df], sort=False)
comp_df.to_csv('./processed_data/interactive_composition_data.csv', index=False)

#%%
ferm_df = pd.DataFrame([])
for i in tqdm.tqdm(range(len(data))):
    d = data.iloc[i]
    g = d['Unnamed: 0']
    study_info = studies[studies['Unnamed: 0']==g] 
    _data_dict = {'study_name': study_info['study_name'].values[0],
                  'subject_id': study_info['subject_id'].values[0],
                  'age_category': study_info['age_category'].values[0],
                  'disease_state': study_info['disease'].values[0], 
                  'bm_fraction': d[[f'familylevel_{f}' for f in families]].sum()}
    # break
    for i, fp in enumerate(fermentation_products):
        _data_dict[fp] = d[f'familylevel_{fp}']
    _df = pd.DataFrame(_data_dict, index=[0])
    ferm_df = pd.concat([ferm_df, _df], sort=False)
ferm_df['total_uptake'] = np.abs(ferm_df['glucose'] + 2 * ferm_df['maltose'])
ferm_df.to_csv('./processed_data/interactive_fermentation_products.csv', index=False)

