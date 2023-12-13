#%%

import numpy as np
import pandas as pd
import bokeh.io 
import bokeh.plotting 
import bokeh.models
import bokeh.layouts
import os
import utils

cor, pal = utils.bokeh_style()
FNAME = './output/interactive_excretion.html'
bokeh.io.output_file(FNAME)
data = pd.read_csv('./processed_data/interactive_fermentation_products.csv')
data.loc[data['study_name'].str.contains('Tett'), 'study_name'] = 'TettAJ_2019'

# Generate mappers for study names and diseases 
study_data = pd.read_csv('./selected_studies_forinteractive_mapper.csv')
study_dict = {k:v for k, v in zip(study_data['study_name'].values, study_data['key'].values)}
age_dict = {'newborn': 'Newborn',
            'child': 'Child',
            'schoolage': 'Schoolage',
            'adult': 'Adult',
            'senior': 'Senior'}  

disease_dict = {
    'healthy': 'Healthy',
    'IBD': 'Inflammatory Bowel Disease',
    'T1D': 'Type I Diabetes',
    'T2D': 'Type II Diabetes',
    'T2D;respiratoryinf': 'Type II Diabetes (+ Resp. Infection)',
    'respiratoryinf': 'Respiratory Infection',
    'adenoma': 'Adenoma',
    'CRC': 'Colorectal Cancer',
    'ACVD': 'Atherosclerotic Cardiovascular Disease',
    'IGT': 'Impaired Glucose Tolerance',
    'IGT;respiratoryinf': 'Impaired Glucose Tolerace (+ Resp. Infection)',
    'STH': 'Soil-Transmitted Helminthiasis',
    'IBD;perianal_fistula': 'Inflammatory Bowel Disease (+ Perianal Fistula)',
    'IBD': 'Inflammatory Bowel Disease',
    'carcinoma_surgery_history': 'History of Carcinoma Surgery',
    'few_polyps': 'Colorectal Polyps)'
}

data['study_name'] = [study_dict[k] for k in data['study_name'].values]
data['age_category'] = [age_dict[k] for k in data['age_category'].values]
data['disease_state'] = [disease_dict[k] for k in data['disease_state'].values]

data.dropna(inplace=True)

study_state_dict = {}
for g, d in data.groupby('study_name'):
    age_dict_ = {}
    for _g, _d in  d.groupby('age_category'):
        age_dict_[_g] = list(_d['disease_state'].unique())
    study_state_dict[g] = age_dict_


study_menu = list(data['study_name'].unique())
study_menu.append('All Studies')

INIT_VALUE = 50
INIT_BINS = 25
GLUCOSE_MASS = 0.18016 #in g per mmol
fps = ['glucose', 'maltose', 'acetate', 
       'formate', 'butyrate', 'propionate', 'lactate',
       'succinate']
for f in fps:
    data[f] = data[f].fillna(0)

data['drymass'] = (INIT_VALUE/GLUCOSE_MASS) / data['total_uptake']

data['total_excretion'] = data[fps[2:]].sum(axis=1)
data['total_excretion_mass'] = data['total_excretion'] * data['drymass']
fps = fps[2:]
data.sort_values(by='drymass', inplace=True)

source = bokeh.models.ColumnDataSource(data)
biomass_hist, biomass_bins = np.histogram(data['drymass'], bins=INIT_BINS)
total_fp_hist, total_fp_bins = np.histogram(data['total_excretion_mass'], bins=INIT_BINS)
composition_hist, composition_bins = np.histogram(data['bm_fraction'], bins=INIT_BINS)

fp_cds_dict = {}
for f in fps:
    hist, bins = np.histogram(data[f] * data['drymass'], bins=INIT_BINS)
    cds = bokeh.models.ColumnDataSource({'top':hist, 'bottom':np.zeros_like(hist),
                                         'left':bins[:-1], 'right':bins[1:]})
    fp_cds_dict[f] = cds

biomass_dist_data = bokeh.models.ColumnDataSource({'top':biomass_hist, 'bottom':np.zeros_like(biomass_hist), 
                                                   'left':biomass_bins[:-1], 'right':biomass_bins[1:]})
total_fp_dist_data = bokeh.models.ColumnDataSource({'top':total_fp_hist, 'bottom':np.zeros_like(total_fp_hist), 
                                                   'left':total_fp_bins[:-1], 'right':total_fp_bins[1:]})
total_composition_dist_data = bokeh.models.ColumnDataSource({'top':composition_hist, 'bottom':np.zeros_like(composition_hist), 
                                                   'left':composition_bins[:-1], 'right':composition_bins[1:]})

# ##############################################################################
# INPUT WIDGET DEFINITION
# ##############################################################################

input_slider = bokeh.models.Slider(start=10, end=100, value=INIT_VALUE, step=0.1,
                            title='Lower-Intestine Carbohydrate Load [g / day]',
                            sizing_mode='stretch_width', bar_color=cor['primary_black'])


study_selector = bokeh.models.Select(value='All Studies', options=study_menu, title='Study')
age_selector = bokeh.models.MultiChoice(value=["Adult"], options=list(data['age_category'].unique()), title='Age category', sizing_mode='stretch_width')
disease_selector = bokeh.models.MultiChoice(value=["Healthy"], options=list(data['disease_state'].unique()), title='Health status', sizing_mode='stretch_width')

# ##############################################################################
# AXIS DEFINITION
# ##############################################################################
biomass_ax = bokeh.plotting.figure(x_axis_label='total bacterial drymass [g]', 
                                   y_axis_label='number of individuals',
                                   width=300, height=260,
                                   x_range=(0, 80),
                                   y_range=(0, 3500))

total_fp_ax = bokeh.plotting.figure(x_axis_label='total daily fermentation product [mmol]', 
                                   y_axis_label='number of individuals',
                                   width=300, height=260,
                                   x_range=(0, 3000),
                                   y_range=(0, 6000))

total_composition_ax = bokeh.plotting.figure(x_axis_label='characterized fraction [%]', 
                                   y_axis_label='number of individuals',
                                   width=450, height=200,
                                   x_range=(0, 100))

fp_ax = {}
for i, f in enumerate(fps):
    if i == 0:
        x_range = (0, 1000)
        y_range = (0, 5000)
    else:
        x_range = fp_ax['acetate'].x_range
        y_range = (0, 3500)
    fp_ax[f] = bokeh.plotting.figure(x_axis_label='total daily amount [mmol / day]',
                                     y_axis_label='number of individuals',
                                     title=f'secreted {f}', width=250, height=180,
                                     x_range=x_range,
                                     y_range=y_range)


################################################################################
# CANVAS POPULATION
################################################################################
fps_cor = {'acetate':cor['dark_blue'], 'formate':cor['gold'], 
           'propionate':cor['light_blue'], 'butyrate':cor['dark_green'],
           'lactate':cor['light_green'], 'succinate':cor['primary_blue']}
nox = ['acetate', 'formate', 'propionate', 'lactate']
for f in fps:
    fp_ax[f].quad(top='top', bottom='bottom', left='left', right='right',
                fill_color=fps_cor[f], line_color=fps_cor[f],
                source=fp_cds_dict[f])
    if f in nox:
        fp_ax[f].xaxis.axis_label = ''
        fp_ax[f].xaxis.axis_line_width = 0
        fp_ax[f].xaxis.minor_tick_in = 0
        fp_ax[f].xaxis.minor_tick_out = 0
        fp_ax[f].xaxis.major_label_text_font_size = "0pt"

biomass_ax.quad(top='top', bottom='bottom', left='left', right='right',
                fill_color=cor['light_green'], line_color=cor['light_green'],
                source=biomass_dist_data)

total_fp_ax.quad(top='top', bottom='bottom', left='left', right='right',
                fill_color=cor['light_purple'], line_color=cor['light_purple'],
                source=total_fp_dist_data)

total_composition_ax.quad(top='top', bottom='bottom', left='left', right='right',
                fill_color=cor['light_black'], line_color=cor['light_black'],
                source=total_composition_dist_data)

################################################################################               
# CALLBACK DEFINITION
################################################################################               
args = {'input_slider': input_slider,
        'source': source,
        'biomass_bin_source': biomass_dist_data,
        'total_fp_bin_source': total_fp_dist_data,
        'fp_cds':fp_cds_dict,
        'study_selector': study_selector,
        'age_selector':age_selector,
        'disease_selector':disease_selector,
        'study_state_dict':study_state_dict,
        'age_dict':age_dict,
        'disease_dict':disease_dict}

cb = utils.load_js('interactive_excretion.js', args=args)
input_slider.js_on_change('value', cb)
study_selector.js_on_change('value', cb)
age_selector.js_on_change('value', cb)
disease_selector.js_on_change('value', cb)

################################################################################
# LAYOUT SPECIFICATION
################################################################################
selector_col = bokeh.layouts.column(study_selector, age_selector, disease_selector)
selector_row = bokeh.layouts.row(selector_col, total_composition_ax)
totals_grid = bokeh.layouts.gridplot([[biomass_ax], [bokeh.layouts.Spacer(height=20)], [total_fp_ax]])
fp_grid = bokeh.layouts.gridplot([[fp_ax['acetate'], fp_ax['propionate']],
                                  [fp_ax['formate'], fp_ax['lactate']],
                                  [fp_ax['butyrate'], fp_ax['succinate']]])
bottom_row = bokeh.layouts.row(totals_grid, bokeh.layouts.Spacer(width=50), fp_grid)
layout = bokeh.layouts.column(selector_row, input_slider, bottom_row)

################################################################################
# EXPORT 
################################################################################
def save_with_d3(layout, fname, jsdelivr=['npm/d3@7', 
                                          'npm/mathjs@12.2.0/lib/browser/math.min.js']):
    bokeh.io.save(layout)
    with open(fname, 'r') as f:
        lines = f.readlines()
    newfile = ""
    added = False
    for l in lines:
        if ('<script' in l) & (added == False):
            for js in jsdelivr:
                script = f'<script src="https://cdn.jsdelivr.net/{js}"></script>\n'
                newfile += script
            added = True  
        newfile += l
    os.remove(fname)
    with open(fname, 'w') as f:
        f.write(newfile)
    print(f"Saved figure as '{FNAME}'.")

save_with_d3(layout, FNAME)
