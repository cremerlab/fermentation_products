let m_carb_input = input_slider.value;
let GLUCOSE_MASS = 0.18016;
let data = source.data;
let binning_fn = d3.bin().thresholds(50)
