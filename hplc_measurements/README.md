
Our experiments include data on growth and raw chromatograms to determine sugar and fermentation product concentrations in culture samples.

# Raw experimental data

The growth data is provided in the csv file growthcurves.csv. For each culture, this file utilizes five columns to provide time (columns 1 and 2), optical density (columns 3 and 4), and the filename of the HPLC sample taken (column 5) during a growth curve. The species name is provided in the column header.

The HPLC raw data is provided in the folder chromatogram_raw_data, with one text file per chromatogram, named as listed in the table growthcurves.csv. 

# HPLC analysis

To analyze the chromatograms and extract concentration values, we utilized the Python pipeline hplc-py (Chure et al, https://doi.org/10.21105/joss.06270)  we had developed for this purpose. This pipeline allows for the efficient fitting of Gaussian peaks to chromatograms, dealing with overlaps of peaks. We recommend following the documentation provided on https://github.com/cremerlab/hplc-py to utilize this pipeline. The scripts we provide here are an extension to run this pipelines for all chromatography samples from one growth curve and analyze observed trends. We  provide an example analysis of this approach in the folder hplc_analysis_example. This includes also the integration of calibration curves. Finally, we also provide in the folder data_calibrationthe chromatograms and their analysis used for calibration. 
