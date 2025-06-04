
Our experiments include data on growth and raw chromatograms to determine sugar and fermentation product concentrations. 

# Raw experimental data

The growth data is provided in the csv file growthcurves.csv. For each culture, this file utilizes five columns to provide time (columns 1 and 2), optical density (columns 3 and 4), and the filename of the HPLC sample taken (column 5). The species name is provided in the column header.

The HPLC raw data is provided in the folder chromatogram_raw_data, with one text file per chromatogram, named as listed in the table growthcurves.csv. 

# HPLC analysis

To analyze the chromatograms and extract concentration values, we utilized the Python pipeline hplc-py (Chure et al, https://doi.org/10.21105/joss.06270)  we had developed for this purpose, allowing for the efficient fitting of Gaussian peaks to chromatograms to deal with overlapping peaks. We recommend following the documentation provided on https://github.com/cremerlab/hplc-py to utilize this pipeline. We also provide an example analysis of our data in the folder hplc_analysis_example, including the integration of growth rate data and calibration curves. Finally, for completeness, we also provide the calibration curves we used in our analysis in the data_calibration folder.

# Illustration of HPLC results

To facilitate the comparative analysis of different species, we provide plots illustrating the trends in optical density, sugar, and fermentation product abundances for each analyzed culture over time.

We provide one pdf per media condition probed. In folder plots_hplc_results:

1. HPLC_trends_YCA.pdf
2. HPLC_trends_BHIS.pdf
3. HPLC_trends_epsilon.pdf
4. HPLC_trends_gamma.pdf