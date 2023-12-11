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
INIT_VALUE = 50
GLUCOSE_MASS = 0.18016 #in g per mmol

data['drymass'] = (INIT_VALUE/GLUCOSE_MASS) / data['total_uptake']
data.sort_values(by='drymass', inplace=True)

input = bokeh.models.Slider(start=0.1, end=100, value=INIT_VALUE, step=0.1, 
                            title='Lower-Intestine Carbohydrate Load [g / day]')

biomass_ax = bokeh.plotting.figure(x_axis_label='total bacterial drymass [g]')
x = np.sort(data['drymass'].values)
y = np.arange(len(data)) / len(data)
biomass_ax.step(x=x, y=y, line_color=cor['primary_blue'])

args = {'input_slider': input,
        'X_VAL': x}
cb = utils.load_js('interactive_excretion.js', args=args)

input.js_on_change('value', cb)
layout = bokeh.layouts.column(input, biomass_ax)

def save_with_d3(layout, fname, jsdelivr='npm/d3@7'):
    bokeh.io.save(layout)
    _jsdelivr = f"https://cdn.jsdelivr.net/{jsdelivr}"
    script = f'<script type="text/javascript" src="{_jsdelivr}"></script>\n'
    with open(fname, 'r') as f:
        lines = f.readlines()
    newfile = ""
    added = False
    for l in lines:
        if ('<script' in l) & (added == False):
            newfile += script
            added = True  
        newfile += l
    os.remove(fname)
    with open(fname, 'w') as f:
        f.write(newfile)
    print(f"Saved figure as '{FNAME}'.")

save_with_d3(layout, FNAME)
