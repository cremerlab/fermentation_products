let GLUCOSE_MASS = 0.18016;
let data = source.data;
let biomass_bin_data = biomass_bin_source.data;
let total_fp_bin_data = total_fp_bin_source.data;
let total_energy_bin_data = total_energy_bin_source.data;
let characterized_bin_data = characterized_bin_source.data;
let study = study_selector.value;
let age = age_selector.value;
let disease = disease_selector.value;
const binningFn = d3.bin();
const charBinningFn = d3.bin().domain([0, 100]);
if (tab_selector.active == 0) {
    var mCarbInput = input_slider.value;
}
else if (tab_selector.active == 1) {
    var mCarbInput = 0.13 * starch_slider.value + 0.5 * fiber_slider.value;
}
else {

    if (diet_selector.active == 0) {
        var mCarbInput = 35.7;
    }
    else if (diet_selector.active == 1) {
        var mCarbInput = 24.6;
    }

    else {
        var mCarbInput = 57.6;
    }
}
// Populate the age and disease selector
let age_menu = []
if (study == "All Studies") {
    age_menu = ['Newborn', 'Child', 'Schoolage', 'Adult', 'Senior']
}
else {
    age_menu = Object.keys(study_state_dict[study]);
}
if (study == "All Studies") {
    var studies = Object.keys(study_state_dict)
} else {
    var studies = [study]
}

let disease_menu = []

for (var i = 0; i < studies.length; i++) {
    let study_ = studies[i];
    for (var j = 0; j < age.length; j++) {
        if (Object.keys(study_state_dict[study_]).includes(age[j])) {
            var diseases_ = study_state_dict[study_][age[j]]
            for (var k = 0; k < diseases_.length; k++) {
                if (disease_menu.includes(diseases_[k]) == false) {
                    disease_menu.push(diseases_[k])
                }
            }
        }

    }
}

// Update the selector options
age_selector.options = age_menu;
disease_selector.options = disease_menu;


// Apply the study and disease selection
let drymass = []
let total_fp = []
let characterized = []
let total_energy = []
let fps_keys = ['acetate', 'propionate', 'formate', 'lactate', 'butyrate', 'succinate']
let fps = { 'acetate': [], 'butyrate': [], 'formate': [], 'lactate': [], 'propionate': [], 'succinate': [] }
let fps_enthalpies = { 'acetate': 0.21, 'butyrate': 0.52, 'formate': 0.0, 'lactate': 0.33, 'propionate': 0.37, 'succinate': 0.36 }
let table_vals = []

for (var i = 1; i < data['drymass'].length; i++) {

    // Determine if the selected study is valid
    if (study == data['study_name'][i] || study == 'All Studies') {
        var valid_study = true
    }
    else {
        var valid_study = false
    }

    // Determine if the selected age is valid
    if (age.includes(data['age_category'][i])) {
        var valid_age = true
    }
    else {
        var valid_age = false
    }

    // Determine if the disease state is valid
    if (disease.includes(data['disease_state'][i])) {
        var valid_disease = true
    }
    else {
        var valid_disease = false
    }

    if (valid_study && valid_age && valid_disease) {
        var energy_ = 0
        let drymass_ = mCarbInput / (GLUCOSE_MASS * data['total_uptake'][i])
        drymass.push(drymass_)
        characterized.push(data['bm_fraction'][i])
        total_fp.push(drymass_ * data['total_excretion'][i])
        for (var j = 0; j < fps_keys.length; j++) {
            let val_ = drymass_ * data[fps_keys[j]][i]
            fps[fps_keys[j]].push(val_);
            energy_ += fps_enthalpies[fps_keys[j]] * val_
        }
        total_energy.push(energy_)
    }
}

table_vals.push(drymass.length)
if (drymass.length > 2) {
    table_vals.push(parseFloat(math.median(...characterized).toPrecision(3)))
}
else {
    table_vals.push('Not enough individuals to compute median.')
}
table_vals.push(parseFloat(mCarbInput.toPrecision(3)))

if (drymass.length > 2) {
    table_vals.push(parseFloat(math.mean(...drymass).toPrecision(3)))
}
else {
    table_vals.push('Not enough individuals to compute mean.')
}
for (var i = 0; i < fps_keys.length; i++) {
    if (drymass.length > 2) {
        table_vals.push(parseFloat(math.mean(...fps[fps_keys[i]]).toPrecision(3)))

    }
    else {
        table_vals.push('Not enough individuals to compute mean.')
    }
}

if (drymass.length > 2) {
    table_vals.push(parseFloat(math.mean(...total_fp).toPrecision(3)))

}
else {
    table_vals.push('Not enough individuals to compute mean.')
}

if (drymass.length > 2) {
    table_vals.push(parseFloat(math.mean(...total_energy).toPrecision(3)))

}
else {
    table_vals.push('Not enough individuals to compute mean.')
}
table_cds.data['Value'] = table_vals
table_cds.change.emit()
for (var i = 0; i < fps_keys.length; i++) {
    var data_ = fp_cds[fps_keys[i]].data;
    let binned = binningFn(fps[fps_keys[i]])
    let bottom = [];
    let top = [];
    let left = [];
    let right = [];
    for (var j = 0; j < binned.length; j++) {
        bottom.push(0);
        top.push(binned[j].length);
        left.push(binned[j].x0);
        right.push(binned[j].x1);
    }

    data_['bottom'] = bottom;
    data_['top'] = top;
    data_['left'] = left;
    data_['right'] = right;
    fp_cds[fps_keys[i]].change.emit();
}

let drymass_binned = binningFn(drymass);
let total_fp_binned = binningFn(total_fp);
let total_energy_binned = binningFn(total_energy);
let characterized_binned = charBinningFn(characterized);
let total_energy_bottom = [];
let total_energy_top = [];
let total_energy_left = [];
let total_energy_right = [];
let characterized_bottom = [];
let characterized_top = [];
let characterized_left = [];
let characterized_right = [];
let drymass_bottom = [];
let drymass_top = [];
let drymass_left = [];
let drymass_right = [];
let total_fp_bottom = [];
let total_fp_top = [];
let total_fp_left = [];
let total_fp_right = [];
for (var i = 1; i < drymass_binned.length; i++) {
    drymass_bottom.push(0);
    drymass_top.push(drymass_binned[i].length);
    drymass_left.push(drymass_binned[i].x0);
    drymass_right.push(drymass_binned[i].x1);
}

for (var i = 1; i < total_fp_binned.length; i++) {
    total_fp_bottom.push(0);
    total_fp_top.push(total_fp_binned[i].length);
    total_fp_left.push(total_fp_binned[i].x0);
    total_fp_right.push(total_fp_binned[i].x1);
}

for (var i = 1; i < characterized_binned.length; i++) {
    characterized_bottom.push(0);
    characterized_top.push(characterized_binned[i].length);
    characterized_left.push(characterized_binned[i].x0);
    characterized_right.push(characterized_binned[i].x1);
}

for (var i = 1; i < total_energy_binned.length; i++) {
    total_energy_bottom.push(0);
    total_energy_top.push(total_energy_binned[i].length);
    total_energy_left.push(total_energy_binned[i].x0);
    total_energy_right.push(total_energy_binned[i].x1);
}



biomass_bin_data['top'] = drymass_top;
biomass_bin_data['bottom'] = drymass_bottom;
biomass_bin_data['left'] = drymass_left;
biomass_bin_data['right'] = drymass_right;

total_fp_bin_data['top'] = total_fp_top;
total_fp_bin_data['bottom'] = total_fp_bottom;
total_fp_bin_data['left'] = total_fp_left;
total_fp_bin_data['right'] = total_fp_right;

characterized_bin_data['top'] = characterized_top;
characterized_bin_data['bottom'] = characterized_bottom;
characterized_bin_data['left'] = characterized_left;
characterized_bin_data['right'] = characterized_right;

total_energy_bin_data['top'] = total_energy_top;
total_energy_bin_data['bottom'] = total_energy_bottom;
total_energy_bin_data['left'] = total_energy_left;
total_energy_bin_data['right'] = total_energy_right;

biomass_bin_source.change.emit();
total_fp_bin_source.change.emit();
characterized_bin_source.change.emit();
total_energy_bin_source.change.emit();