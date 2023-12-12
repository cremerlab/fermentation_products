let mCarbInput = input_slider.value;
let GLUCOSE_MASS = 0.18016;
let data = source.data;
let biomass_bin_data = biomass_bin_source.data;
let total_fp_bin_data = total_fp_bin_source.data;
let drymass = []
let total_fp = []
const binningFn = d3.bin();
let fps = ['acetate', 'butyrate', 'formate', 'lactate', 'propionate', 'succinate']
for (var i = 1; i < data['drymass'].length; i++) {
    var drymass_ = mCarbInput / (GLUCOSE_MASS * data['total_uptake'][i])
    drymass.push(drymass_)
    total_fp.push(drymass_ * data['total_excretion'][i])
}


console.log(data[fps[1]])
for (var i = 0; i < fps.length; i++) {
    var data_ = fp_cds[fps[i]].data;
    let excretion_ = []
    for (var j = 0; j < drymass.length; j++) {
        excretion_.push(drymass[j] * data[fps[i]][j]);
    }
    let binned = binningFn(excretion_)
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
    ata_['right'] = right;
    fp_cds[fps[i]].change.emit();
}

let drymass_binned = binningFn(drymass);
let total_fp_binned = binningFn(total_fp);
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

biomass_bin_data['top'] = drymass_top;
biomass_bin_data['bottom'] = drymass_bottom;
biomass_bin_data['left'] = drymass_left;
biomass_bin_data['right'] = drymass_right;

total_fp_bin_data['top'] = total_fp_top;
total_fp_bin_data['bottom'] = total_fp_bottom;
total_fp_bin_data['left'] = total_fp_left;
total_fp_bin_data['right'] = total_fp_right;


biomass_bin_source.change.emit();
total_fp_bin_source.change.emit();
