#%%
import numpy as np 
import pandas as pd 
import tqdm
data = pd.read_csv('../analysis/data_analysisresults/secretion_microbiomestudies.csv')
studies = pd.read_csv('../analysis/data_curated_microbiome/relabundance_colData.csv')
families = ['Bacteroidaceae',
            'Bifidobacteriaceae',
            'Coriobacteriaceae',
            'Enterobacteriaceae',
            'Lachnospiraceae',
            'Oscillospiraceae',
            'Prevotellaceae',
            'Tannerellaceae']

fermentation_products = ['acetate', 'butyrate', 'formate', 'lactate',
                         'propionate', 'succinate']

#%%
comp_df = pd.DataFrame([])
for g, d in tqdm.tqdm(data.groupby('Unnamed: 0')):
    study_info = studies[studies['Unnamed: 0']==g]
    _comp_df = pd.DataFrame([]) 
    
    for i, f in enumerate(families):
        __comp_df = pd.DataFrame({'study_name': study_info['study_name'].values[0],
               'subject_id': study_info['subject_id'].values[0],
               'age_category': study_info['age_category'].values[0],
               'disease_state': study_info['disease'].values[0],
               'bacterial_family': f,
               'bm_fraction':d[f'familylevel_{f}'].values[0]/100}, index=[0])
        _comp_df = pd.concat([_comp_df, __comp_df])
    tot = _comp_df['bm_fraction'].sum()
    __comp_df = pd.DataFrame({'study_name': study_info['study_name'].values[0],
                         'subject_id': study_info['subject_id'].values[0],
                         'age_category': study_info['age_category'].values[0],
                         'disease': study_info['disease'].values[0],
                         'bacterial_family': 'not characterized',
                         'bm_fraction':1-tot}, index=[0])
    comp_df = pd.concat([comp_df, _comp_df, __comp_df], sort=False)
#%%
comp_df.to_csv('./processed_data/family_composition.csv', index=False)

#%%
comp_df = comp_df[comp_df['disease_state']=='healthy']
comp_df.to_csv('./processed_data/family_composition_healthy.csv', index=False)

#%%
prod_df = pd.DataFrame()
for g, d in tqdm.tqdm(data.groupby('Unnamed: 0')):
    study_info = studies[studies['Unnamed: 0']==g]
    for i, fp in enumerate(fermentation_products):
        for j, fam in enumerate(families):
            _prod_df = pd.DataFrame({
                'study_name': study_info['study_name'].values[0],
                'subject_id': study_info['subject_id'].values[0],
                'age_category': study_info['age_category'].values[0],
                'disease_state': study_info['disease'].values[0],
                'ferm_product': fp,
                'bacterial_family':fam,
                'excretion': d[f'familylevel_{fp}'].values[0]
        }, index=[0])
        prod_df = pd.concat([prod_df, _prod_df])
#%%
prod_df.to_csv('./processed_data/family_excretion.csv', index=False)
prod_df = prod_df[prod_df['disease_state']=='healthy']
prod_df.to_csv('./processed_data/family_excretion_healthy.csv', index=False)
#%%
