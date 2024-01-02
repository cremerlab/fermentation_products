# %%

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
data = data[data['disease_state'] != 'respiratoryinf']
data.loc[data['study_name'].str.contains('Tett'), 'study_name'] = 'TettAJ_2019'

# Generate mappers for study names and diseases
study_data = pd.read_csv('./selected_studies_forinteractive_mapper.csv')
study_dict = {k: v for k, v in zip(
    study_data['study_name'].values, study_data['key'].values)}
fps_enthalpies = {'acetate': 0.21, 'butyrate': 0.52, 'formate': 0.0,
                  'lactate': 0.33, 'propionate': 0.37, 'succinate': 0.36}
# tot = np.sum(tot, axis=0)
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
    'T2D;respiratoryinf': 'Type II Diabetes',
    'respiratoryinf': 'Respiratory Infection',
    'adenoma': 'Adenoma',
    'CRC': 'Colorectal Cancer',
    'ACVD': 'Atherosclerotic Cardiovascular Disease',
    'IGT': 'Impaired Glucose Tolerance',
    'IGT;respiratoryinf': 'Impaired Glucose Tolerace',
    'STH': 'Soil-Transmitted Helminthiasis',
    'IBD;perianal_fistula': 'Inflammatory Bowel Disease',
    'carcinoma_surgery_history': 'History of Carcinoma Surgery',
    'few_polyps': 'Colorectal Polyps'
}

data['study_name'] = [study_dict[k] for k in data['study_name'].values]
data['age_category'] = [age_dict[k] for k in data['age_category'].values]
data['disease_state'] = [disease_dict[k] for k in data['disease_state'].values]
data.dropna(inplace=True)
study_state_dict = {}
for g, d in data.groupby('study_name'):
    age_dict_ = {}
    for _g, _d in d.groupby('age_category'):
        age_dict_[_g] = list(_d['disease_state'].unique())
    study_state_dict[g] = age_dict_


study_menu = list(data['study_name'].unique())
study_menu.append('All Studies')

INIT_VALUE = 50
INIT_BINS = 25
GLUCOSE_MASS = 0.18016  # in g per mmol
fps = ['glucose', 'maltose', 'acetate',
       'formate', 'butyrate', 'propionate', 'lactate',
       'succinate']
for f in fps:
    data[f] = data[f].fillna(0)

data['drymass'] = (INIT_VALUE/GLUCOSE_MASS) / data['total_uptake']
tot = [data[k] * v * data['drymass'] for k, v in fps_enthalpies.items()]
data['total_energy'] = np.sum(tot, axis=0)
data['total_excretion'] = data[fps[2:]].sum(axis=1)
data['total_excretion_mass'] = data['total_excretion'] * data['drymass']
fps = fps[2:]
data.sort_values(by='drymass', inplace=True)

source = bokeh.models.ColumnDataSource(data)
_data = data[(data['age_category'] == 'Adult') &
             (data['disease_state'] == 'Healthy')]
biomass_hist, biomass_bins = np.histogram(_data['drymass'], bins=INIT_BINS)
total_fp_hist, total_fp_bins = np.histogram(
    _data['total_excretion_mass'], bins=INIT_BINS)
composition_hist, composition_bins = np.histogram(
    _data['bm_fraction'], bins=INIT_BINS)
energy_hist, energy_bins = np.histogram(_data['total_energy'], bins=INIT_BINS)


table_dict = {'Quantity': ['Number of individuals', 'Median Characterized Fraction', 'Microbiota Accessible Carbohydrates', 'Average Bacterial Drymass Growth',
                           'Average Acetate', 'Average Propionate', 'Average Formate', 'Average Lactate',
                           'Average Butyrate', 'Average Succinate', 'Average Total Fermentation Products', 'Average Daily Energy From Fermentation Products'],
              'Units': ['#', '%', 'g / person / day', 'g / person / day', 'mmol / person/ day',
                        'mmol / person/ day', 'mmol / person / day',
                        'mmol / person / day', 'mmol / person / day',
                        'mmol / person / day', 'mmol / person / day',
                        'kcal / person / day'],
              'Value': [len(_data),
                        round(_data['bm_fraction'].median(), 2),
                        INIT_VALUE,
                        round(_data['drymass'].mean(), 2),
                        round(_data['acetate'].mean(), 2),
                        round(_data['propionate'].mean(), 2),
                        round(_data['formate'].mean(), 2),
                        round(_data['lactate'].mean(), 2),
                        round(_data['butyrate'].mean(), 2),
                        round(_data['succinate'].mean(), 2),
                        round(_data['total_excretion_mass'].mean(), 2),
                        round(_data['total_energy'].mean(), 2)
                        ]}
table_cds = bokeh.models.ColumnDataSource(table_dict)

fp_cds_dict = {}
for f in fps:
    hist, bins = np.histogram(_data[f] * _data['drymass'], bins=INIT_BINS)
    cds = bokeh.models.ColumnDataSource({'top': hist, 'bottom': np.zeros_like(hist),
                                         'left': bins[:-1], 'right': bins[1:]})
    fp_cds_dict[f] = cds

biomass_dist_data = bokeh.models.ColumnDataSource({'top': biomass_hist, 'bottom': np.zeros_like(biomass_hist),
                                                   'left': biomass_bins[:-1], 'right': biomass_bins[1:]})
total_fp_dist_data = bokeh.models.ColumnDataSource({'top': total_fp_hist, 'bottom': np.zeros_like(total_fp_hist),
                                                   'left': total_fp_bins[:-1], 'right': total_fp_bins[1:]})
total_composition_dist_data = bokeh.models.ColumnDataSource({'top': composition_hist, 'bottom': np.zeros_like(composition_hist),
                                                             'left': composition_bins[:-1], 'right': composition_bins[1:]})
total_energy_dist_data = bokeh.models.ColumnDataSource({'top': energy_hist, 'bottom': np.zeros_like(energy_hist),
                                                        'left': energy_bins[:-1], 'right': energy_bins[1:]})

# ##############################################################################
# INPUT WIDGET DEFINITION
# ##############################################################################

input_slider = bokeh.models.Slider(start=10, end=100, value=INIT_VALUE, step=0.1,
                                   title='Total Microbiota Available Carbohydrate Load [g / day]',
                                   sizing_mode='stretch_width', bar_color=cor['primary_black'])
starch_slider = bokeh.models.Slider(
    start=10, end=1000, value=200, title='consumed starch [g / day]', bar_color=cor['primary_black'], width=400)
fiber_slider = bokeh.models.Slider(
    start=10, end=60, value=25, title='consumed fiber [g / day]', bar_color=cor['primary_black'], width=400)
diet_selector = bokeh.models.RadioButtonGroup(labels=["British Reference Diet", "NHANES US Diet", "Hadza Diet"], active=None,
                                              sizing_mode='stretch_width')

study_selector = bokeh.models.Select(
    value='All Studies', options=study_menu, title='Study')
age_selector = bokeh.models.MultiChoice(value=["Adult"], options=list(
    data['age_category'].unique()), title='Age category', width=400)
disease_selector = bokeh.models.MultiChoice(value=["Healthy"], options=list(
    data['disease_state'].unique()), title='Health status', width=400)

mac_input_panes = bokeh.models.TabPanel(
    child=input_slider, title="Adjust Total Microbiota Available Carbohydrates")
starch_fiber_composition = bokeh.models.TabPanel(child=bokeh.layouts.row(
    starch_slider, fiber_slider, sizing_mode='stretch_width'), title="Adjust Consumed Starch and Fiber")
diet_radio = bokeh.models.TabPanel(
    child=diet_selector, title='Select Reference Diet')
tab = bokeh.models.Tabs(
    tabs=[mac_input_panes, starch_fiber_composition, diet_radio], sizing_mode='stretch_width')

# ##############################################################################
# AXIS DEFINITION
# ##############################################################################
biomass_ax = bokeh.plotting.figure(x_axis_label='total bacterial drymass growth\n[g / person / day]',
                                   y_axis_label='number of individuals',
                                   width=300, height=180,)
#    x_range=(0, 80))

total_fp_ax = bokeh.plotting.figure(x_axis_label='total daily fermentation product\n[mmol / person / day]',
                                    y_axis_label='number of individuals',
                                    width=300, height=180,)
#    x_range=(0, 3000))

total_composition_ax = bokeh.plotting.figure(x_axis_label='characterized fraction of biomass [%]',
                                             y_axis_label='number of individuals',
                                             width=450, sizing_mode='stretch_height')
total_energy_ax = bokeh.plotting.figure(x_axis_label='fermentation product energy\ncontribution [kcal / person / day]',
                                        y_axis_label='number of individuals',
                                        width=300,  height=180)


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
                                     title=f'secreted {f}', width=250, height=180)
    #  x_range=x_range)


################################################################################
# CANVAS POPULATION
################################################################################
fps_cor = {'acetate': cor['light_blue'], 'formate': cor['light_gold'],
           'propionate': cor['light_green'], 'butyrate': cor['light_purple'],
           'lactate': cor['pale_blue'], 'succinate': cor['light_black']}
nox = ['acetate', 'formate', 'propionate', 'lactate']
for f in fps:
    fp_ax[f].quad(top='top', bottom='bottom', left='left', right='right',
                  fill_color=fps_cor[f], line_color='white',
                  source=fp_cds_dict[f])
biomass_ax.quad(top='top', bottom='bottom', left='left', right='right',
                fill_color=cor['black'], line_color='white',
                source=biomass_dist_data)

total_fp_ax.quad(top='top', bottom='bottom', left='left', right='right',
                 fill_color=cor['black'], line_color='white',
                 source=total_fp_dist_data)

total_composition_ax.quad(top='top', bottom='bottom', left='left', right='right',
                          fill_color=cor['black'], line_color='white',
                          source=total_composition_dist_data)

total_energy_ax.quad(top='top', bottom='bottom', left='left', right='right',
                     fill_color=cor['primary_red'], line_color='white',
                     source=total_energy_dist_data)
#
################################################################################
# TABLE DEFINITION
################################################################################
columns = [
    bokeh.models.TableColumn(field='Quantity', title='Quantity'),
    bokeh.models.TableColumn(field='Value', title='Value'),
    bokeh.models.TableColumn(field='Units', title='Unit')
]

data_table = bokeh.models.DataTable(
    source=table_cds, columns=columns, width=850, height=350)


################################################################################
# CALLBACK DEFINITION
################################################################################
args = {'input_slider': input_slider,
        'source': source,
        'biomass_bin_source': biomass_dist_data,
        'total_fp_bin_source': total_fp_dist_data,
        'characterized_bin_source': total_composition_dist_data,
        'fp_cds': fp_cds_dict,
        'study_selector': study_selector,
        'age_selector': age_selector,
        'disease_selector': disease_selector,
        'study_state_dict': study_state_dict,
        'age_dict': age_dict,
        'disease_dict': disease_dict,
        'table_cds': table_cds,
        'total_energy_bin_source': total_energy_dist_data,
        'tab_selector': tab,
        'starch_slider': starch_slider,
        'fiber_slider': fiber_slider,
        'diet_selector': diet_selector}

cb = utils.load_js('interactive_excretion.js', args=args)
input_slider.js_on_change('value', cb)
study_selector.js_on_change('value', cb)
age_selector.js_on_change('value', cb)
disease_selector.js_on_change('value', cb)
tab.js_on_change('active', cb)
starch_slider.js_on_change('value', cb)
fiber_slider.js_on_change('value', cb)
diet_selector.js_on_event('button_click', cb)

################################################################################
# LAYOUT SPECIFICATION
################################################################################
selector_col = bokeh.layouts.column(
    study_selector, age_selector, disease_selector)
selector_row = bokeh.layouts.row(selector_col, total_composition_ax)
totals_grid = bokeh.layouts.gridplot(
    [[biomass_ax], [total_fp_ax], [total_energy_ax]])
fp_grid = bokeh.layouts.gridplot([[fp_ax['acetate'], fp_ax['propionate']],
                                  [fp_ax['formate'], fp_ax['lactate']],
                                  [fp_ax['butyrate'], fp_ax['succinate']]])
bottom_row = bokeh.layouts.row(totals_grid, bokeh.layouts.Spacer(width=15), bokeh.models.Div(
    text='<div style="border-left:1px solid #afafaf;height:550px"></div>', sizing_mode='stretch_height', width=10), fp_grid)
row_div1 = bokeh.models.Div(
    text='<div style="border-bottom:1px solid #afafaf;width:860px"></div>')
row_div2 = bokeh.models.Div(
    text='<div style="border-bottom:1px solid #afafaf;width:860px"></div>')
layout = bokeh.layouts.column(selector_row, row_div1, tab, bokeh.models.Div(
    text='<div style="border-bottom:1px solid #afafaf;width:860px"></div>'), bottom_row, row_div2, data_table)

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
