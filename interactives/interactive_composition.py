#%%
import numpy as np 
import pandas as pd 
import bokeh.plotting 
import bokeh.io 
import bokeh.models
import bokeh.palettes
from utils import bokeh_style, load_js
cor, pal = bokeh_style() 

bokeh.io.output_file('./output/interactive_composition.html')
data = pd.read_csv('./processed_data/family_composition_healthy.csv')
_data = data[-80*9:]
fam = _data['bacterial_family'].unique()
subj = _data['subject_id'].unique()
color = []
data_dict = {g:d['bm_fraction'].values for g, d in _data.groupby('bacterial_family')}
data_dict['subj'] = np.arange(len(subj))
data_dict

p = bokeh.plotting.figure(width=1000, x_range=subj, tools="hover", tooltips="$name @subject_id: @$name",
                          y_axis_label="fractional composition",
                          x_axis_label="subject ID")
                          
p.vbar_stack(fam, x='subj', width=0.9, source=data_dict, 
             color=bokeh.palettes.Category10_9)

p.y_range.start = 0
p.x_range.range_padding = 0.1
p.xgrid.grid_line_color = None
p.axis.minor_tick_line_color = None
p.outline_line_color = None
# p.legend.location = "top_left"
bokeh.io.show(p)
bokeh.io.save(p)