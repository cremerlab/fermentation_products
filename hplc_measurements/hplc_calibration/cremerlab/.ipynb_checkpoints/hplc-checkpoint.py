from time import sleep
import pandas as pd 
import numpy as np
from io import StringIO
import shutil
import scipy.signal
import scipy.optimize
import scipy.special
import tqdm
import os 
import json
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import glob 

def scrape_metadata(file, delimiter=','):
    """
    Scrapes the sample information from the output of a Shimadzu HPLC ASCII 
    file and returns a dictionary of the metadata entries. 
    Parameters
    ----------
    file : str 
        The contents of the file as a string with newline characters (`\\n`) present. 
    delimiter: str 
        The delimeter used in the file. If  tab-delimited, use `\\t`.
    Returns
    -------
    metadata_dict : dictionary
        A dictionary of the metadata listed under the 'Sample Information', excluding 
        ISTD Amounts.
    Notes
    -----
    This function assumes that the file contains metadata fields `[Sample Information]`
    followed by `[Original Files]`. If either `[Sample Information]` or `[Original Files]`
    fields are missing or have other names, a ValueError exception is thrown. 
    Raises
    ------
    TypeError
        Raised if file is not of type `str`
    ValueError:
        Raised if `[Sample Information]`, `[Original Files]`, `\\n`, or delimiter
        is not present in filebounpar
        
    """
    # Make sure the 'file' provided is a string and has necessary fields. 
    if type(file) != str:
        raise TypeError(f'Argument `file` must be a string. type {type(file)} was provided.')
    if '\n' not in file:
        raise ValueError(f'Newline characters (`\\n`) are not in file, but must be for proper parsing.')
    if ('[Sample Information]' not in file) | ('[Original Files]' not in file):
        raise ValueError('`[Sample Information]` or `[Original Files]` field is missing.')
    if delimiter not in file:
        raise ValueError(f'Delimiter {f} not present in the file.')

    # Get the sample information and split by newline
    metadata = file.split('[Sample Information]\n')[1].split('[Original Files]\n')[0]
    
    # Set up the dictionary and loop through lines
    metadata_dict = {}
    for line in metadata.split('\n'):

        # Split by the delimiter.
        entry = line.split(delimiter)
        
        # Add the entry to the dictionary if longer than 1 and does not 
        # contain ISTD Amount.
        if (len(entry) > 1) & ('ISTD' not in entry[0]):
            metadata_dict[entry[0]] = entry[1]
    return metadata_dict

def scrape_chromatogram(file, detector='B', delimiter=',', metadata=True):
    """
    Scrapes the chromatogram for a given detector from a Shimadzu HPLC
    ASCII output.
    Parameters
    ----------
    file : str 
        The contents of the file as a string with newline characters (`\\n`) present. 
    detector: str, 'A' or 'B'
        The name of the detector in the file. Default is `B`. Note that only that only 
        'A' or 'B' is an acceptable input and not the complete detector name such as 
        'LC Chromatogram(Detector B-Ch1)'.
    delimiter: str 
        The delimeter used in the file. If  tab-delimited, use `\\t`.
    metadata : bool
        If `true`, a dictionary with the metadata about the detector is returned. 
        Default is True. 
    Returns
    -------
    chrom : pandas DataFrame
        A tidy pandas DataFrame with two columns -- `time_minutes` and `intensity_mV`
    metadata_dict : dictionary
        A dictionary of the metadata associated with the desired detector channel.
        if `metadata` is not `True`, this is not returned. 
    Notes
    -----
    This function assumes that the detector name follows the convention 
    `[LC Chromatogram(Detector A/B-ChX]` where `A/B` is the detector label 
    and `X` is the channel number. 
    Raises
    ------
    TypeError :
        Raised if file, detctor, or delimiter is not of type `str`.
    ValueError:
        Raised if `[LC Chromatogram(Detector`, `\\n`, or delimiter is not present in file.
        Also raised if `detector.upper()` is not `A` or `B`
    """
    # Do type checks. 
    for arg in (file, detector, delimiter):
        if type(arg) is not str:
            raise TypeError(f'Type of `{arg}` must be `str`. Type {type(arg)} was provided.')
    for st in ['[LC Chromatogram(Detector', '\n', delimiter]:
        if st not in file:
            raise ValueError(f'Pattern {st} not present in file.')
    if (detector.upper() != 'A') & (detector.upper() != 'B'):
        raise ValueError(f'Detector must be `A` or `B`. String of `{detector.upper()}` was provided.')

    # Parse and split the file to get the chromatogram
    chrom = file.split(
                    f'[LC Chromatogram(Detector {detector.upper()}-Ch'
                    )[1]
    if '[' in chrom:
        chrom = chrom.split('[')[0]
    
    chrom = '\n'.join(chrom.split('\n')[1:-2])
    # Read into as a dataframe.
    df = pd.read_csv(StringIO(chrom), delimiter=delimiter, skiprows=6)

    # Rename the columns
    df.columns = ['time_min', 'intensity_mV']
    df['detector'] = detector

    # Dropnas
    df.dropna(inplace=True)
    # Determine if the metadata should be scraped and returned
    if metadata:
        metadata_dict = {}
        for line in chrom.split('\n')[:6]:
            entry = line.split(delimiter)
            if len(entry) > 0:
                metadata_dict[entry[0]] = entry[1]
        out = [df, metadata_dict]
    else:
        out = df
    return out

def scrape_peak_table(file, detector='B', delimiter=','):
    """
    Scrapes the Shimadzu-generated peak table output from the ASCII output.

    Parameters
    ----------
    file : str 
        The contents of the file as a string with newline characters (`\\n`) present. 
    detector: str, 'A' or 'B'
        The name of the detector in the file. Default is `B`. Note that only that only 
        'A' or 'B' is an acceptable input and not the complete detector name such as 
        'Peak Table(Detector B)'. 
    delimiter: str 
        The delimeter used in the file. If  tab-delimited, use `\\t`.

    Returns
    -------
    peaks : pandas DataFrame
        A tidy pandas DataFrame with the identified peaks, times, areas, and 
        heights.

    Notes
    -----
    This function assumes that the detector name follows the convention 
    `[Peak Table(Detector A/B]` where `A/B` is the detector label 
    and `X` is the channel number. 

    Raises
    ------
    TypeError :
        Raised if file, detctor, or delimiter is not of type `str`.
    ValueError:
        Raised if `[LC Chromatogram(Detector`, `\\n`, or delimiter is not present in file.
        Also raised if `detector.upper()` is not `A` or `B`
    """
    # Do type checks. 
    for arg in (file, detector, delimiter):
        if type(arg) is not str:
            raise TypeError(f'Type of `{arg}` must be `str`. Type {type(arg)} was provided.')
    for st in ['[Peak Table(Detector', '\n', delimiter]:
        if st not in file:
            raise ValueError(f'Pattern {st} not present in file.')
    if (detector.upper() != 'A') & (detector.upper() != 'B'):
        raise ValueError(f'Detector must be `A` or `B`. String of `{detector.upper()}` was provided.')

    # Parse and split the file to get the chromatogram
    peaks = file.split(
                    f'[Peak Table(Detector {detector.upper()}'
                    )[1]
    if '[' in peaks:
        peaks = peaks.split('[')[0]
    
    peaks = '\n'.join(peaks.split('\n')[1:])
    
    # Read into as a dataframe.
    df = pd.read_csv(StringIO(peaks), delimiter=delimiter, skiprows=1)

    # Rename the columns
    df = df[['Peak#', 'R.Time', 'I.Time', 'F.Time', 'Area', 'Height']]
    df.rename(columns={'Peak#':'peak_idx', 'R.Time':'retention_time', 
                       'I.Time': 'arrival_time', 'F.Time':'departure_time',
                       'Area':'area', 'Height':'height'}, inplace=True)
    df['detector'] = detector

    # Dropnas
    df.dropna(inplace=True)

    # Determine if the metadata should be scraped and returned
    return df 

def convert(file_path, detector='B', delimiter=',', peak_table=False, 
            output_dir=None,  save_prefix=None, save_suffix=None, 
            verbose=True, overwrite_output_dir=False,hplc_machine=None, use_samplename_for_output=False):
    """
    Reads the ASCII output from a Shimadzu HPLC and returns a DataFrame of the
    chromatogram. Converted files can also be saved to disk with human-readable 
    metadata. 
    Parameters
    ----------
    file_path : str or list of str
        The file path of the ASCII file to read. Multiple files can be provided 
        as a list of file paths. 
    detector : str or list of str, ['A' or 'B']
        The detector channel chromatogram to return. Must be either 'A' or 'B' 
        or both, if provided as a list. Default is just 'B'.
    delimiter: str 
        The delimiter used in the file. If  tab-delimited, use `\t`. Default is 
        a comma, `,`.
    peak_table: bool
        If True, the peak table is also parsed and is saved in the directory 
        with the extra suffix `_peaks`
    output_dir : str or None
        The output directory (provided as a global path, starting with '/') if 
        the dataframes are to be saved to disk. If `None` and `save = True`,
        the dataframe will be saved in the directory of the `file` in
        an `output` folder.
    save_prefix : str
        A prefix to append to the file name if saved to disk. If None, 
        saved filed will just be `SAMPLE_NAME.csv`.       
    save_suffix : str
        A suffix to append to the file name if saved to disk. If None, 
        saved filed will just be `SAMPLE_NAME.csv`.       
    verbose : bool 
        If True a progress bar will print if there's more than 5 files to 
        process.
    overwrite_output_dir : bool
        If True, empty `output_dir` and populate exclusively with the converted files
    """
    if type(file_path) is not list:
            file_path = [file_path]
            
            
    if hplc_machine in [None,"Shimadzu_ProminenceLC2030C"]:
        # Determine the size of the  
        

        # Determine if there should be a verbose output
        if (verbose is True) & (len(file_path) >= 5):
            iterator = enumerate(tqdm.tqdm(file_path, desc='Converting ASCII output.'))
        else:
            iterator = enumerate(file_path)

        #TODO: Make sure this works on a windows system
        if output_dir is None:
            output_path = '/'.join(file_path[0].split('/')[:-1])
            if os.path.isdir(f'{output_path}/converted'):
                shutil.rmtree(f'{output_path}/converted')
            if os.path.isdir(f'{output_path}/converted') == False:
                os.mkdir(f'{output_path}/converted')
            output_path += '/converted'
        else:
            output_path = output_dir
            if os.path.isdir(output_path)==True:
                if overwrite_output_dir:
                    if os.path.isdir(output_path):
                        shutil.rmtree(output_path)
                else:
                    raise FileExistsError(output_path)
            else:
                # recursively makes all intermediate dirs, if missing
                os.makedirs(output_path)

        for _, f in iterator:
          try:
            with open(f, 'r') as file:
                raw_file = file.read()
                # Get the file metadata 
                file_metadata = scrape_metadata(raw_file, delimiter=delimiter)

                # Get chromatogram
                if type(detector) == str:
                    detector = [detector]
                chroms = []
                peaks = []
                chrom_metadata = {}
                chrom_file_metadata = []
                for d in detector:

                    # Parse the chromatogram
                    _chrom, _chrom_metadata = scrape_chromatogram(raw_file, 
                                                                detector=d, 
                                                                delimiter=delimiter)    
                    if peak_table:
                        _peaks = scrape_peak_table(raw_file,
                                                   detector=d,
                                                   delimiter=delimiter)
                        peaks.append(_peaks)
                    chroms.append(_chrom)
                    chrom_metadata[f'Detector {d}'] = _chrom_metadata

                    # Generate the metadata for file saving
                    _metadata = f"""#
# Detector {d.upper()}
# ---------------------
# Acquisition Interval: {_chrom_metadata['Interval(msec)']} ms
# Intensity Units: {_chrom_metadata['Intensity Units']}
# Intensity Multiplier: {_chrom_metadata['Intensity Multiplier']}
#
"""
                    chrom_file_metadata.append(_metadata)
                if len(detector) == 1: 
                    chrom = chroms[0]
                    if peak_table:
                        peaks = peaks[0]
                else:
                    chrom = pd.concat(chrom, sort=False)
                    if peak_table:
                        peaks = pd.concat(peaks, sort=False)
                if type(detector) == str:
                    detector = [detector]

                # Assemble the comment file 
                header = f"""#
# {file_metadata['Sample Name']}
# ------------------------------
# Acquired: {file_metadata['Acquired']}
# Converted: {datetime.now().strftime('%m/%d/%Y %I:%m:%S %p')}
# Vial: {file_metadata['Vial#']}
# Injection Volume: {file_metadata['Injection Volume']} uL
# Sample ID: {file_metadata['Sample ID']}
"""
                for m in chrom_file_metadata:
                    header += m
                # Process prefix and suffix to the save file
                name = ''
                if save_prefix is not None:
                    name += save_prefix + '_'
                if use_samplename_for_output:
                        name += file_metadata['Sample Name']
                else:
                        name += os.path.basename(f)[:-4]
                        
                if save_suffix is not None:
                    name += '_' + save_suffix
                if peak_table:
                    peak_name = name + '_peak_table.csv'
                    peak_save_name = f'{output_path}/{peak_name}'
                name += '_chromatogram.csv'
                save_name = f'{output_path}/{name}'
                print("save name")
                print(save_name)
                

                if os.path.isfile(save_name):
                    exists = Truep
                    iter = 0
                    while exists:
                        new_name = f'{save_name.split(".csv")[0]}_{iter}.csv'
                        if peak_table:
                            new_peak_name = f'{peak_save_name.split(".csv")[0]}_{iter}.csv'
                        if os.path.isfile(new_name):
                            iter += 1
                        else:
                            exists = False
                    save_name = new_name
                    peak_save_name = new_peak_name
                with open(save_name, 'a') as save_file:
                    save_file.write(header)
                    chrom.to_csv(save_file, index=False) 
                    
                if peak_table:
                    with open(peak_save_name, 'a') as save_file:
                        save_file.write(repr(header))
                        peaks.to_csv(save_file, index=False) 
          except:
            print("file "+f+" could not be converted")
        print(f'Converted file(s) saved to `{output_path}`')


    elif hplc_machine=="Zurich_ThermoUltimate_3000":
        foldername=os.path.dirname(file_path[0])
        os.makedirs(os.path.join(foldername,'converted'), exist_ok=True)
        for i in  range (len(file_path)):

            #Create dataframe from the output file
            raw_data = pd.read_csv(file_path[i], sep="\t", names=['time_min','step','intensity_mV'])
            raw_data = raw_data.drop(raw_data.index[range(40)])
            for x in range (2400):
                raw_data.iloc[x,2] = raw_data.iloc[x,2].replace("- ", "-")

            raw_data = raw_data.astype(float)

            #Create a new folder with csv files with only the chromatogram values
            raw_data.to_csv(os.path.join(os.path.join(foldername,'converted'),os.path.basename(file_path[i]).replace('.txt','')+'_chromatogram.csv'), index = False)
    else:
        error_machine_not_found


class Chromatogram(object):
    """
    Base class for dealing with HPLC chromatograms
    """
    def __init__(self, file=None, time_window=None, baselinecorrection=False,
                    cols={'time':'time_min', 'intensity':'intensity_mV'},
                    csv_comment='#'):
        """
        Instantiates a chromatogram object on which peak detection and quantification
        is performed.
        Parameters
        ----------
        file: str or pandas DataFrame, optional
            The path to the csv file of the chromatogram to analyze or 
            the pandas DataFrame of the chromatogram. If None, a pandas DataFrame 
            of the chromatogram must be passed.
        dataframe : pandas DataFrame, optional
            a Pandas Dataframe of the chromatogram to analyze. If None, 
            a path to the csv file must be passed
        time_window: list [start, end], optional
            The retention time window of the chromatogram to consider for analysis.
            If None, the entire time range of the chromatogram will be considered.
        cols: dict, keys of 'time', and 'intensity', optional
            A dictionary of the retention time and intensity measurements 
            of the chromatogram. Default is 
            `{'time':'time_min', 'intensity':'intensity_mV'}`.
        csv_comment: str, optional
            Comment delimiter in the csv file if chromatogram is being read 
            from disk.
        """

        # Peform type checks and throw exceptions where necessary. 
        if file is None:
            raise RuntimeError(f'File path or dataframe must be provided')
        if (type(file) is not str) & (type(file) is not pd.core.frame.DataFrame):
            raise RuntimeError(f'Argument must be either a filepath or pandas DataFrame. Argument is of type {type(file)}')
        if (time_window is not None):
            if type(time_window) != list:
                raise TypeError(f'`time_window` must be of type `list`. Type {type(time_window)} was proivided')
            if len(time_window) != 2:
                raise ValueError(f'`time_window` must be of length 2 (corresponding to start and end points). Provided list is of length {len(time_window)}.')

        # Assign class variables 
        self.time_col = cols['time']
        self.int_col = cols['intensity']

        # Load the chromatogram and necessary components to self. 
        if type(file) is str:
            dataframe = pd.read_csv(file, comment='#')
        else:
            dataframe = file 
        self.df = dataframe

        # Prune to time window
        if time_window is not None:
            self.crop(time_window)
        else: 
            self.df = dataframe

         # Correct for a negative baseline 
        df = self.df
        min_int = df[self.int_col].min() 
        if baselinecorrection:
            intensity = df[self.int_col] - min_int

        # Blank out vars that are used elsewhere
        self.window_df = None
        self.window_props = None
        self.peaks = None
        self.peak_df = None

    def crop(self, time_window=None, return_df=False):
        """
        Restricts the time dimension of the DataFrame
        Parameters
        ----------
        time_window : list [start, end], optional
            The retention time window of the chromatogram to consider for analysis.
            If None, the entire time range of the chromatogram will be considered.
        return_df : bool
            If `True`, the cropped DataFrame is 
        Returns
        -------
        cropped_df : pandas DataFrame
            If `return_df = True`, then the cropped dataframe is returned.
        """
        if type(time_window) != list:
                raise TypeError(f'`time_window` must be of type `list`. Type {type(time_window)} was proivided')
        if len(time_window) != 2:
                raise ValueError(f'`time_window` must be of length 2 (corresponding to start and end points). Provided list is of length {len(time_window)}.')
        self.df = self.df[(self.df[self.time_col] >= time_window[0]) & 
                          (self.df[self.time_col] <= time_window[1])]
        if return_df:
            return self.df

    def backgroundsubstraction(self,num_iterations=10, return_df=False):
        """
        Substracts background for entire chromatogram using the algorithm by Miroslav Morhac et al 

        Parameters
        ----------
        num_iterations : int
            The number of iterations to run. For each iteration one additional pixel is included and one should chose the number of iterations as the typical width of a peak
        return_df : bool
            If `True`, then chromatograms (before and after background correction) are returned
        Returns
        -------
        corrected_df : pandas DataFrame
            If `return_df = True`, then the original and the corrected chromatogram are returned.
        """
        df = self.df
        try:
            intensity = self.df[self.int_col+"_nobackgroundcorrection"].values
        except:
             intensity = self.df[self.int_col].values
            
        
        intensity_old=intensity.copy()
        intensity = intensity*np.heaviside(intensity,0)
        #transform to log scale
        intensity_transf=np.log(np.log(np.sqrt(intensity+1)+1)+1)
        #start itteration
        for il in range(0,num_iterations):
            intensity_transf_new=intensity_transf.copy()
            for i in range(il,intensity_transf.shape[0]-il):
                intensity_transf_new[i]=min(intensity_transf[i],0.5*(intensity_transf[i+il]+intensity_transf[i-il]))
            intensity_transf=intensity_transf_new
        #transform back
        intensity=np.power(np.exp(np.exp(intensity_transf)-1.)-1.,2.)-1.
        self.df[self.int_col]= intensity_old-intensity
        self.df[self.int_col+'_nobackgroundcorrection']=intensity_old
        self.df[self.int_col+'_background']=intensity
        
        
        if return_df:
            return self.df
        
        
    def _assign_peak_windows(self, prominence, rel_height, buffer, manual_peak_positions=None,baselinecorrection=False):
        """
        Breaks the provided chromatogram down to windows of likely peaks. 
        Parameters
        ----------
        prominence : float,  [0, 1]
            The promimence threshold for identifying peaks. Prominence is the 
            relative height of the normalized signal relative to the local
            background. Default is 1%.
        rel_height : float, [0, 1]
            The relative height of the peak where the baseline is determined. 
            Default is 95%.
        buffer : positive int
            The padding of peak windows in units of number of time steps. Default 
            is 100 points on each side of the identified peak window.
        manual_peak_positions : list of floats
            to provide peak position by hand instead of estimating peak positions via script
        Returns
        -------
        windows : pandas DataFrame
            A Pandas DataFrame with each measurement assigned to an identified 
            peak or overlapping peak set. This returns a copy of the chromatogram
            DataFrame with  a column  for the local baseline and one column for 
            the window IDs. Window ID of -1 corresponds to area not assigned to 
            any peaks
        """
        for param, param_name, param_type in zip([prominence, rel_height, buffer], 
                                     ['prominence', 'rel_height',  'buffer'],
                                     [float, float, int]):
            if type(param) is not param_type:
                raise TypeError(f'Parameter {param_name} must be of type `{param_type}`. Type `{type(param)}` was supplied.') 
        if (prominence < 0) | (prominence > 1):
            raise ValueError(f'Parameter `prominence` must be [0, 1].')
        if (rel_height < 0) | (rel_height > 1):  
            raise ValueError(f'Parameter `rel_height` must be [0, 1].')
        if (buffer < 0):
            raise ValueError('Parameter `buffer` cannot be less than 0.')

        # Correct for a negative baseline 
        df = self.df
        intensity = self.df[self.int_col].values
        norm_int = (intensity - intensity.min()) / (intensity.max() - intensity.min())



        # Identify the peaks and get the widths and baselines
        if manual_peak_positions==None:
            peaks, _ = scipy.signal.find_peaks(norm_int, prominence=prominence)
        else:
            timew = self.df[self.time_col].values
            deltatw= (timew[-1]-timew[0])/float(timew.shape[0])
            peaks = np.int_((np.array(manual_peak_positions)-timew[0])/deltatw) #peeak position in descrete step of time-window
            #print("peaks detected at times : "+str(self.df[self.time_col].values[peaks]))
            # 
            #print(timew[peaks]) #time of detected peaks
        if len(peaks)>0 and peaks[0]==0: #do not use first peak if it is exactly at left boundrary #maybe better solution
            peaks=peaks[1:]
        self.peaks_inds = peaks
        #need to fix, not needed for manual peak position
        out = scipy.signal.peak_widths(intensity, peaks, 
                                       rel_height=rel_height)
        _, heights, left, right = out
        widths, _, _, _ = scipy.signal.peak_widths(intensity, peaks, 
                                       rel_height=0.5)
        
        
        ###
        # Set up the ranges
        ranges = []
        for l, r in zip(left, right):
            if (l - buffer) < 0:
                l = 0
            else:
                l -= buffer
            if (r + buffer) > len(norm_int):
                r = len(norm_int)
            else:
                r += buffer
            ranges.append(np.arange(np.round(l), np.round(r), 1))

        # Identiy subset ranges and remove
        valid = [True] * len(ranges)
        for i, r1 in enumerate(ranges):
            for j, r2 in enumerate(ranges):
                if i != j:
                    if set(r2).issubset(r1):
                        valid[j] = False
        
        # Keep only valid ranges and baselines
        ranges = [r for i, r in enumerate(ranges) if valid[i] is True]
        baselines = [h for i, h in enumerate(heights) if valid[i] is True] 

        # Copy the dataframe and return the windows
        
        window_df = df.copy(deep=True)
        window_df.sort_values(by=self.time_col, inplace=True)
        window_df['time_idx'] = np.arange(len(window_df))
        for i, r in enumerate(ranges):
            window_df.loc[window_df['time_idx'].isin(r), 
                                    'window_idx'] = int(i + 1)
            window_df.loc[window_df['time_idx'].isin(r), 
                                    'baseline'] = baselines[i]
        window_df.dropna(inplace=True) 
        
        # Convert this to a dictionary for easy parsing
        window_dict = {}
        
        time_step = np.mean(np.diff(self.df[self.time_col].values))
        for g, d in window_df.groupby('window_idx'):
                _peaks = [p for p in peaks if (p in d['time_idx'].values) and (d[d['time_idx']==p][self.int_col].values[0] >0)] #ignores peaks where intensity is smaller zero
                peak_inds = [x for _p in _peaks for x in np.where(peaks == _p)[0]]
                if baselinecorrection:
                    _dict = {'time_range':d[self.time_col].values,
                         'intensity': d[self.int_col] - baselines[i], #? is this the good correction to make
                         'intensity_nobaselinecorrection': d[self.int_col], #added
                         'num_peaks': len(_peaks),
                         'amplitude': [d[d['time_idx']==p][self.int_col].values[0] - baselines[i] for p in _peaks],
                         'amplitude_nobaselinecorrection': [d[d['time_idx']==p][self.int_col].values[0] for p in _peaks],
                         'location' : [d[d['time_idx']==p][self.time_col].values[0] for p in _peaks],
                         'width' :    [widths[ind] * time_step for ind in peak_inds]
                         }
                else:
                    _dict = {'time_range':d[self.time_col].values,
                         'intensity': d[self.int_col], #added
                         'num_peaks': len(_peaks),
                         'amplitude': [d[d['time_idx']==p][self.int_col].values[0] for p in _peaks],
                         'location' : [d[d['time_idx']==p][self.time_col].values[0] for p in _peaks],
                         'width' :    [widths[ind] * time_step for ind in peak_inds]
                         }
                window_dict[g] = _dict
        self.window_props = window_dict
        return window_df  

    def _compute_skewnorm(self, x, *params):
        R"""
        Computes the lineshape of a skew-normal distribution given the shape,
        location, and scale parameters
        Parameters
        ----------
        x : float or numpy array
            The time dimension of the skewnorm 
        params : list, [amplitude, loc, scale, alpha]
            Parameters for the shape and scale parameters of the skewnorm 
            distribution.
                amplitude : float; > 0
                    Height of the peak.
                loc : float; > 0
                    The location parameter of the distribution.
                scale : float; > 0
                    The scale parameter of the distribution.
                alpha : float; > 
                    THe skew shape parater of the distribution.
        Returns
        -------
        scaled_pdf : float or numpy array, same shape as `x`
            The PDF of the skew-normal distribution scaled with the supplied 
            amplitude.
        Notes
        -----
        This function infers the parameters defining skew-norma distributions 
        for each peak in the chromatogram. The fitted distribution has the form 
            
        .. math:: 
            I = 2I_\text{max} \left(\frac{1}{\sqrt{2\pi\sigma^2}}\right)e^{-\frac{(t - r_t)^2}{2\sigma^2}}\left[1 + \text{erf}\frac{\alpha(t - r_t)}{\sqrt{2\sigma^2}}\right]
        where :math:`I_\text{max}` is the maximum intensity of the peak, 
        :math:`t` is the time, :math:`r_t` is the retention time, :math:`\sigma`
        is the scale parameter, and :math:`\alpha` is the skew parameter.
        """
        amp, loc, scale, alpha = params
        _x = alpha * (x - loc) / scale
        normfactor=1#np.sqrt(2 * np.pi * scale**2)**-1 * 
        norm = normfactor*np.exp(-(x - loc)**2 / (2 * scale**2))
        cdf = (1 + scipy.special.erf(_x / np.sqrt(2))) 
        return amp * norm * cdf

    def _fit_skewnorms(self, x, *params):
        R"""
        Estimates the parameters of the distributions which consititute the 
        peaks in the chromatogram. 
        Parameters
        ----------
        x : float
            The time dimension of the skewnorm 
        params : list of length 4 x number of peaks, [amplitude, loc, scale, alpha]
            Parameters for the shape and scale parameters of the skewnorm 
            distribution. Must be provided in following order, repeating
            for each distribution.
                amplitude : float; > 0
                    Height of the peak.
                loc : float; > 0
                    The location parameter of the distribution.
                scale : float; > 0
                    The scale parameter of the distribution.
                alpha : float; > 
                    THe skew shape parater of the distribution.
        Returns
        -------
        out : float
            The evaluated distribution at the given time x. This is the summed
            value for all distributions modeled to construct the peak in the 
            chromatogram.
        """
        # Get the number of peaks and reshape for easy indexing
        n_peaks = int(len(params) / 4)
        params = np.reshape(params, (n_peaks, 4))
        out = 0
        
        # Evaluate each distribution
        for i in range(n_peaks):
            out += self._compute_skewnorm(x, *params[i])
        return out
        
    def _estimate_peak_params(self, boundpars=None, verbose=True,baselinecorretion=False):
        R"""
        For each peak window, estimate the parameters of skew-normal distributions 
        which makeup the peak(s) in the window.  
        Parameters
        ----------
        boundpars : list
            Used to provide boundaries for peak fittng function
            8 boundaries provided (default value in parenthesis)
                1. lower boundary amplitude (0) 
                2. lower boundary time window; difference to peak position
                3. lower boundary peak width (0.0)
                4. lower boundary skew parameter (-np.inf)
                5. upper boundary amplitude (np.inf) 
                6. upper boundary time window; difference to peak position
                7. upper boundary peak width (np.inf)
                8. upper boundary skew parameter (np.inf)
                
        
        verbose : bool
            If `True`, a progress bar will be printed during the inference.
        """ 
        if self.window_props is None:
            raise RuntimeError('Function `_assign_peak_windows` must be run first. Go do that.')
        if verbose:
            iterator = tqdm.tqdm(self.window_props.items(), desc='Fitting peak windows...')  
        else:
            iterator = self.window_props.items()
        peak_props = {}
        for k, v in iterator:
            window_dict = {}
            # Set up the initial guess
            p0 = [] 
            bounds = [[],  []] 
            for i in range(v['num_peaks']):
                
                #if v['amplitude'][i]<0: #ignore peaks with negative amplitude
                #    print("peak "+str(i)+"amplitude: "+str(v['amplitude'][i])+" at position: "+str(v['location'][i]))
                #    print("Warning: negative amplitude. Value before reset: "+str(v['amplitude'][i])+" at time "+str(v['location'][i]))
                #    v['amplitude'][i]=0
                #if v['amplitude_nobaselinecorrection'][i]<0: #ignore peaks with negative amplitude
                #    print("peak "+str(i)+"amplitude: "+str(v['amplitude'][i])+" at position: "+str(v['location'][i]))
                #    print("Warning: negative amplitude. Value before reset: "+str(v['amplitude_nobaselinecorrection'][i])+" at time "+str(v['location'][i]))
                #    v['amplitude_nobaselinecorrection'][i]=0
                p0.append(v['amplitude'][i]) #before w/ baseline correction
                
                    
                p0.append(v['location'][i]),
                p0.append(max(v['width'][i] / 4.,0.3)) # scale parameter
                p0.append(0) # Skew parameter, starts with assuming Gaussian

                #REMOVE OPTION NOBASELINE..
                # Set boundaries of fitting
                if boundpars==None: #standard values
                    #lower boundaries
                    bounds[0].append(v['amplitude'][i]*0.5) #min amplitude
                    bounds[0].append(v['time_range'].min()) #min peak position
                    bounds[0].append(0) #min width
                    bounds[0].append(-np.inf) #min skew parameter
                    #upper boundaries
                    bounds[1].append(v['amplitude'][i]*2+2) #max amplitude #?before with baselinecorrection
                    bounds[1].append(v['time_range'].max()) #max peak position
                    bounds[1].append(np.inf) #max width
                    bounds[1].append(np.inf) #skew parameter
                else:
                        #lower boundaries
                        bounds[0].append(v['amplitude'][i]*boundpars[0])
                        bounds[0].append(v['location'][i]-boundpars[1])
                        bounds[0].append(boundpars[2])
                        bounds[0].append(boundpars[3])
                        #upper boundaries
                        bounds[1].append(v['amplitude'][i]*boundpars[4]+2)
                        bounds[1].append(v['location'][i]+boundpars[5]) 
                        bounds[1].append(boundpars[6]) 
                        bounds[1].append(boundpars[7]) 
                        
                        
                
            # Perform the inference
            if len(p0)>0:
                try:
                    #print("p0")
                    #print(p0)
                    #print("bounds")
                    #print(bounds)
                    #print(np.array(p0)-np.array(bounds[0]))
                    #print(np.array(bounds[1])-np.array(p0))
                    popt, _ = scipy.optimize.curve_fit(self._fit_skewnorms, v['time_range'],
                                                   v['intensity'], p0=p0, bounds=bounds,
                                                   maxfev=int(1E6))
                    # Assemble the dictionary of output 
                    if v['num_peaks'] > 1:
                            popt = np.reshape(popt, (v['num_peaks'], 4)) 
                    else:
                            popt = [popt]
                    for i, p in enumerate(popt):
                            window_dict[f'peak_{i + 1}'] = {
                                        'retention_time_firstguess' : v['location'][i],
                                        'amplitude': p[0],
                                        'retention_time': p[1],
                                        'std_dev': p[2],
                                        'alpha': p[3],
                                        'area':self._compute_skewnorm(v['time_range'], *p).sum()*float(v['time_range'][-1]-[v['time_range'][0]])/float(v['time_range'].shape[0]-1)}

                    peak_props[k] = window_dict
                except RuntimeError:
                    print('Warning: Parameters could not be inferred for one peak') #? or there is no peak in that window
                    print(p0)
                    print(v['time_range'].max())
                    print(v['time_range'].min())
                    print(v['intensity'])
            else:
                pass
                #print("Warning: Window without any peak to search for") 
        self.peak_props = peak_props
        return peak_props

    
    
    def quantify(self, time_window = None, prominence = 1E-3, rel_height = 1.0,
                 buffer = 100, manual_peak_positions=None, boundpars=None, peakpositionsonly=False, verbose = True):
        R"""
        Quantifies peaks present in the chromatogram
        Parameters
        ----------
        time_window: list [start, end], optional
            The retention time window of the chromatogram to consider for analysis.
            If None, the entire time range of the chromatogram will be considered.
        prominence : float,  [0, 1]
            The promimence threshold for identifying peaks. Prominence is the 
            relative height of the normalized signal relative to the local
            background. Default is 1%.
        rel_height : float, [0, 1]
            The relative height of the peak where the baseline is determined. 
            Default is 95%.
        buffer : positive int
            The padding of peak windows in units of number of time steps. Default 
            is 100 points on each side of the identified peak window. 
        manual_peak_positions : list
            Manual peak positions of chromatogram. When not provided, autodetection of peaks
        peakpositonsonly : bool
            If ture, only the peak positions will be provided and no full analys of peaks is included
        verbose : bool
            If True, a progress bar will be printed during the inference.
        Returns
        -------
        peak_df : pandas DataFrame
            A dataframe containing information for each detected peak. If peakpositononly=True, this is just a list of peak positions
        Notes
        -----
        This function infers the parameters defining skew-norma distributions 
        for each peak in the chromatogram. The fitted distribution has the form 
            
        .. math:: 
            I = 2I_\text{max} \left(\frac{1}{\sqrt{2\pi\sigma^2}}\right)e^{-\frac{(t - r_t)^2}{2\sigma^2}}\left[1 + \text{erf}\frac{\alpha(t - r_t)}{\sqrt{2\sigma^2}}\right]
        where :math:`I_\text{max}` is the maximum intensity of the peak, 
        :math:`t` is the time, :math:`r_t` is the retention time, :math:`\sigma`
        is the scale parameter, and :math:`\alpha` is the skew parameter.
        """
        
        if time_window is not None:
            dataframe = self.df
            self.df = dataframe[(dataframe[self.time_col] >= time_window[0]) & 
                              (dataframe[self.time_col] <= time_window[1])].copy(deep=True) 
        
        # Assign the window bounds (contains peak autodetection)
        _ = self._assign_peak_windows(prominence, rel_height, buffer, manual_peak_positions=manual_peak_positions)

        #stop script here if only peak positions are wanted
        if peakpositionsonly:
            peakpositions=[]
            #print(self.window_props)
            iterator = self.window_props.items()
            peak_props = {}
            for k, v in iterator: #go through every window
                window_dict = {}
                for i in range(v['num_peaks']):
                    peakpositions.append(v['location'][i])
            return peakpositions
            
        
        # Infer the distributions for the peaks
        peak_props = self._estimate_peak_params(boundpars=boundpars,verbose=verbose)
        
        # Set up a dataframe of the peak properties
        peak_df = pd.DataFrame([])
        iter = 0 
        for _, peaks in peak_props.items():
            for _, params in peaks.items():
                _dict = {'retention_time': params['retention_time'],
                         'retention_time_firstguess': params['retention_time_firstguess'],
                         'scale': params['std_dev'],
                         'skew': params['alpha'],
                         'amplitude': params['amplitude'],
                         'area': params['area'],
                         'peak_idx': iter + 1}     
                iter += 1
                peak_df = pd.concat([peak_df,pd.Series(_dict).to_frame().T],ignore_index=True)
                
                peak_df['peak_idx'] = peak_df['peak_idx'].astype(int)
        self.peak_df = peak_df

        # Compute the mixture
        time = self.df[self.time_col].values
        out = np.zeros((len(time), len(peak_df)))
        iter = 0
        for _k , _v in self.peak_props.items():
            for k, v in _v.items():
                params = [v['amplitude'], v['retention_time'], 
                          v['std_dev'], v['alpha']]
                #print(params)
                #print(time)
                #print(self._compute_skewnorm(time, *params))
                
                out[:, iter] = self._compute_skewnorm(time, *params)
                iter += 1
        self.mix_array = out
        return peak_df

    def show(self):
        """
        Displays the chromatogram with mapped peaks if available.
        """
        sns.set()

        # Set up the figure    
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        ax.set_xlabel(self.time_col)
        ax.set_ylabel(self.int_col)

        # Plot the raw chromatogram
        ax.plot(self.df[self.time_col], self.df[self.int_col], 'k-', lw=2,
                label='raw chromatogram') 

        # Compute the skewnorm mix 
        if self.peak_df is not None:
            time = self.df[self.time_col].values
            # Plot the mix
            convolved = np.sum(self.mix_array, axis=1)
            ax.plot(time, convolved, 'r--', label='inferred mixture') 
            for i in range(len(self.peak_df)):
                ax.fill_between(time, self.mix_array[:, i], label=f'peak {i+1}', 
                                alpha=0.5)
        ax.legend(bbox_to_anchor=(1,1))
        fig.patch.set_facecolor((0, 0, 0, 0))
        return [fig, ax]

def batch_plot(file_paths, time_window=None, cols={'time':'time_min', 'intensity':'intensity_mV'},plot_upper_limit=None,plot_lower_limit=None, plot_time_window=None,detection_peakpositions=None, plot_verticallines=None,backgroundsubstraction=False, backgroundsubstraction_iterations=80,plot_characteristictimes=False, plot_output=False):
    """
    T plot chromatograms (without any analysis). Use to check chromatograms and identify peak poisitions by hand. Data must first
    be converted to a tidy long-form CSV file by using `cremerlab.hplc.convert`

    Parameters
    ----------
    file_paths : list of str
        A list of the file paths.
    time_window : list, optional
        The upper and lower and upper time bounds to consider for the analysis.
        Default is `None`, meaning the whole chromatogram is considered.
    cols: dict, keys of 'time', and 'intensity', optional
        A dictionary of the retention time and intensity measurements of the
        chromatogram. Default is `{'time':'time_min', 'intensity':'intensity_mV'}`.
    plot_upper_limit: upper limit of y-axis, auto-scaling when not provided
    plot_lower_limit: lower limit of y-axis, auto-scaling when not provided. upper_limit must be defined when used.
    plot_time_window: time range 
    detection_peakpositions: Set to 'peakautodetection' to detect peaks
    plot_verticallines: List of retentiontimes to plot as verticallines
    backgroundsubstraction : bool
        Substract background signal following Morhac & Matousek (2008). 
        The original and the corrected chromatogram are shown. If peak detection is done it is based on the 
        corrected chromatogram
    backgroundsubstraction_iterations : int
        Number of neighboring time-points to use for background substraction
    plot_characteristictimes : [list_name_charact, list_time_range_characteristics]
    Returns
    --------
    chrom_df : pandas DataFrame
        A pandas DataFrame  of all of the chromatograms, indexed by file name
    fig : matplotlib.figure.Figure
        Matplotlib figure object for the chromatograms
    ax :  matplotlib AxesSubplot
        The axes of the matplotlib figure
    """

    # Instantiate storage lists
    chrom_dfs, mixes, peakpositions_list= [], [], []

    # Perform the processing for each file
    for i, f in enumerate(tqdm.tqdm(file_paths, desc='Processing files...')):
        # Generate the sample id
        if '/' in f:
            file_name = f.split('/')[-1]
        else:
            file_name = f

        # Check for common file name extension
        for pat in ['.csv', '.txt']:
            if pat in file_name:
                file_name = file_name.split(pat)[0]
                continue

        # Parse teh chromatogram and quantify the peaks
        chrom = Chromatogram(f, cols=cols, time_window=time_window)
        
        if backgroundsubstraction:
            chrom.backgroundsubstraction(num_iterations=backgroundsubstraction_iterations)
        
        if detection_peakpositions == 'peakautodetection':
            peakpositionsc = chrom.quantify(verbose=False, peakpositionsonly=True)
            peakpositions_list.append(peakpositionsc)
        
        # Set up the dataframes for chromatograms and peaks
        _df = chrom.df
        _df['sample'] = file_name
        
        chrom_dfs.append(chrom.df)
        #mixes.append(chrom.mix_array)

    # Concateante the dataframe
    chrom_df = pd.concat(chrom_dfs, sort=False)
    
    #peak_df = pd.concat(peak_dfs, sort=False)

    # Determine the size of the figure
    num_traces = len(chrom_df['sample'].unique())
    num_cols = int(3)
    num_rows = int(np.ceil(num_traces / num_cols))
    unused_axes = (num_cols * num_rows) - num_traces

    # Instantiate the figure
    fig, ax = plt.subplots(num_rows, num_cols, figsize=(6 * num_cols, 4 * num_rows))
    ax = ax.ravel()
    

    for a in ax:
        a.xaxis.set_tick_params(labelsize=6)
        a.yaxis.set_tick_params(labelsize=6)
        a.set_ylabel(cols['intensity'], fontsize=6)
        a.set_xlabel(cols['time'], fontsize=6)
    for i in range(unused_axes):
        ax[-(i + 1)].axis('off')

    # Assign samples to axes.
    mapper = {g: i for i, g in enumerate(chrom_df['sample'].unique())}

    # Plot the chromatogram
    for g, d in chrom_df.groupby(['sample']):
        ax[mapper[g]].set_title(' '.join(g.split('_')), fontsize=12)
        if backgroundsubstraction:
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']], 'b-', lw=1.5,
                           label='after BG correct.')
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']+"_nobackgroundcorrection"], 'y--', lw=1.5,
                           label='original')
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']+"_background"], color='m',ls=':', lw=1.5,
                           label='background')
            ax[mapper[g]].axvline(d[cols['time']].iloc[backgroundsubstraction_iterations],color='b',ls=':',lw=1.5)
        else:
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']], 'b-', lw=1.5,
                           label='original')
        
        cintensity=d[cols['intensity']].to_numpy()
        curtime=d[cols['time']].to_numpy()
        #find maximum
        if plot_upper_limit==None:
            if plot_time_window==None:
                plot_upper_limit=1.1*np.nanmax(d[cols['intensity']])
            else:
                indexfc= np.where((curtime >= plot_time_window[0]) & (curtime <= plot_time_window[1]))
                plot_upper_limit=1.1*np.nanmax(cintensity[indexfc])
                
        #find minimum
        if plot_lower_limit==None:
            if plot_time_window==None:
                plot_lower_limit= np.nanmin(cintensity)
            else:
                indexfc= np.where((curtime >= plot_time_window[0]) & (curtime <= plot_time_window[1]))
                plot_lower_limit=1.1*np.nanmin(cintensity[indexfc])
                
            if plot_lower_limit>0:
                plot_lower_limit=0.9*plot_lower_limit
            else:
                plot_lower_limit=1.1*plot_lower_limit
        if plot_time_window != None:
            ax[mapper[g]].set_xlim(plot_time_window)
            
            
        
        
        
        if detection_peakpositions != None:
            plcvc=-1
            for plcv in peakpositions_list[mapper[g]]:
                plcvc=plcvc+1
                colorlist=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']
                colorc=colorlist[plcvc%10]
                ax[mapper[g]].axvline(plcv,ls='--',c=colorc)
            #display("peaks "+str(g)+" "+str(peakpositions_list[mapper[g]]))
               
                
        if plot_verticallines != None and len(plot_verticallines)>0:
            if type(plot_verticallines[0])==list: #list of lists
                curval=plot_verticallines[mapper[g]]
            else:
                curval=plot_verticallines
            plcvc=-1
            for plcv in curval:
                plcvc=plcvc+1
                colorlist=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']
                colorc=colorlist[plcvc%10]
                ax[mapper[g]].axvline(plcv,ls=':',c=colorc,lw=2)
                
    #mark characteristic times in the plots
    if plot_characteristictimes != False:
        charl=plot_characteristictimes[0]
        timelc=plot_characteristictimes[1]
        timelc2=plot_characteristictimes[2]
        for ill in range(0,len(charl)):
            for a in ax:
                a.text((timelc[ill]+timelc2[ill])/2.,.99,charl[ill],transform=matplotlib.transforms.blended_transform_factory(a.transData, a.transAxes),color='red')
              
    plt.tight_layout()
    fig.patch.set_facecolor((0, 0, 0, 0))
    ax[0].legend(fontsize=6)
    if plot_output != False:
        fig.savefig(plot_output)
    return [chrom_df, [fig, ax], peakpositions_list]

def batch_process(file_paths, time_window=None,  show_viz=False,
                  cols={'time':'time_min', 'intensity':'intensity_mV'},
                 plot_lower_limit=None, plot_upper_limit=None, plot_time_window=None, manual_peak_positions=None,backgroundsubstraction=True,backgroundsubstraction_iterations=80,plot_characteristictimes=False,plot_output=False,
                  **kwargs):
    """
    Performs complete quantification of a set of HPLC data. Data must first 
    be converted to a tidy long-form CSV file by using `cremerlab.hplc.convert`
    Parameters
    ----------
    file_paths : list of str
        A list of the file paths.
    time_window : list, optional
        The upper and lower and upper time bounds to consider for the analysis.
        Default is `None`, meaning the whole chromatogram is considered.
    show_viz : bool 
        If `True`, the plot is printed to the screen. If `False`, a plot is 
        generated and returned, just not shown.
    cols: dict, keys of 'time', and 'intensity', optional
        A dictionary of the retention time and intensity measurements of the 
        chromatogram. Default is `{'time':'time_min', 'intensity':'intensity_mV'}`.
    manual_peak_positions : list
        Manual position of peaks (when not provided, autodetection of peak positions)
        Alternative a list of a list of peaks can be provided - used for the different samples processed during the batch run
    plot_upper_limit: float
        for plotting; upper limit of y-axis, auto-scaling when not provided
    plot_lower_limit: float
        for plotting; lower limit of y-axis, auto-scaling when not provided. upper_limit must be defined when used.
    plot_time_window: float
        for plotting; lower limit of y-axis, auto-scaling when not provided. upper_limit must be defined when used.
    backgroundsubstraction : bool
        Substract background signal following Morhac & Matousek (2008). To show comparison of original and corrected chromatogram use
        batch_plot function. 
    backgroundsubstraction_iterations : int
        Number of neighboring time-points to use for background substraction
    plot_characteristictimes : [list_name_charact, list_time_range_characteristics]
    kwargs: dict, **kwargs
        **kwargs for the peak quantification function `cremerlab.hplc.Chromatogram.quantify`
    Returns 
    --------
    chrom_df : pandas DataFrame
        A pandas DataFrame  of all of the chromatograms, indexed by file name
    peak_df : pandas DataFrame
        A pandas DataFrame of all identified peaks, indexed by file name 
    fig : matplotlib.figure.Figure
        Matplotlib figure object for the chromatograms
    ax :  matplotlib AxesSubplot
        The axes of the matplotlib figure
    """

    # Instantiate storage lists
    chrom_dfs, peak_dfs, mixes = [], [], []

    # Perform the processing for each file
    for i, f in enumerate(tqdm.tqdm(file_paths, desc='Processing files...')):
        # Generate the sample id
        if '/' in f:
            file_name = f.split('/')[-1]
        else:
            file_name = f

        # Check for common file name extension
        for pat in ['.csv', '.txt']:
            if pat in file_name:
                file_name = file_name.split(pat)[0]
                continue

        # Parse teh chromatogram and quantify the peaks
        chrom = Chromatogram(f, cols=cols, time_window=time_window)
        
        if backgroundsubstraction:
            chrom.backgroundsubstraction(num_iterations=backgroundsubstraction_iterations)
        
        if manual_peak_positions==None or type(manual_peak_positions[0]) != list:
            peaks = chrom.quantify(verbose=False, manual_peak_positions=manual_peak_positions, **kwargs)
        else:
            peaks = chrom.quantify(verbose=False, manual_peak_positions=manual_peak_positions[i], **kwargs)
        # Set up the dataframes for chromatograms and peaks
        _df = chrom.df
        _df['sample'] = file_name
        peaks['sample'] = file_name
        peak_dfs.append(peaks)
        
        #display(peaks)
        chrom_dfs.append(chrom.df)
        mixes.append(chrom.mix_array)

    # Concateante the dataframe
    chrom_df = pd.concat(chrom_dfs, sort=False)
    peak_df = pd.concat(peak_dfs, sort=False) 

    # Determine the size of the figure
    num_traces = len(chrom_df['sample'].unique())
    num_cols = int(3)
    num_rows = int(np.ceil(num_traces / num_cols))
    unused_axes = (num_cols * num_rows) - num_traces

    # Instantiate the figure
    fig, ax = plt.subplots(num_rows, num_cols, figsize=(6 * num_cols, 4 * num_rows))
    
    ax = ax.ravel()
    for a in ax:
        a.xaxis.set_tick_params(labelsize=6)
        a.yaxis.set_tick_params(labelsize=6)
        a.set_ylabel(cols['intensity'], fontsize=6)
        a.set_xlabel(cols['time'], fontsize=6)
    for i in range(unused_axes):
        ax[-(i + 1)].axis('off')
        
    # Assign samples to axes.
    mapper = {g: i for i, g in enumerate(chrom_df['sample'].unique())} 

    # Plot the chromatogram
    for g, d in chrom_df.groupby(['sample']): 
        if backgroundsubstraction:
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']], 'b-', lw=1.5,
                           label='after BG correct.')
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']+'_nobackgroundcorrection'], 'y--', lw=1.5,
                           label='original')
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']+"_background"], color='m',ls=':', lw=1.5,
                           label='background')
        else:
            ax[mapper[g]].plot(d[cols['time']], d[cols['intensity']], 'b-', lw=1.5,
                           label='original')
        ax[mapper[g]].set_title(' '.join(g.split('_')), fontsize=12)
        
    # Plot the mapped peaks
    curgroupby=peak_df.groupby(['sample'])
    for g, d in peak_df.groupby(['sample']):
        mix = mixes[mapper[g]] 
        #display(g)
        #display(d)
        #display(mix)
        
        convolved = np.sum(mix, axis=1)
        for i in range(len(d)): 
            _m = np.array(mix[:, i])
            colorlist=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']
            colorc=colorlist[i%10]
            time = np.linspace(time_window[0], time_window[1], len(_m))
            ax[mapper[g]].fill_between(time, 0, _m, alpha=0.5, label=f'peak {i + 1}',color=colorc)
            try:
                ax[mapper[g]].axvline(d.at[i,'retention_time'],ls='--',alpha=1,color=colorc) #print retention time of fit
                #ax[mapper[g]].axhline(d.at[i,'amplitude'],ls=':',alpha=0.5,color=colorc) #print retention time of first peak est.
                
            except:
                pass
        time = np.linspace(time_window[0], time_window[1], len(convolved))
        
        #find maximum
        if plot_upper_limit==None:
            if plot_time_window==None:
                plot_upper_limit=1.1*np.nanmax(convolved)
            else:
                indexfc= np.where((time >= plot_time_window[0]) & (time <= plot_time_window[1]))
                plot_upper_limit=1.1*np.nanmax(convolved[indexfc])
                
        #find minimum
        if plot_lower_limit==None:
            if plot_time_window==None:
                plot_lower_limit= np.nanmin(convolved)
            else:
                indexfc= np.where((time >= plot_time_window[0]) & (time <= plot_time_window[1]))
                plot_lower_limit=1.1*np.nanmin(convolved[indexfc])
                
            if plot_lower_limit>0:
                plot_lower_limit=0.9*plot_lower_limit
            else:
                plot_lower_limit=1.1*plot_lower_limit
        if plot_time_window != None:
            ax[mapper[g]].set_xlim(plot_time_window)
        
        ax[mapper[g]].plot(time, convolved, '--', color='red', lw=2, 
                            label=f'inferred mixture')
        
        ax[mapper[g]].set_ylim(plot_lower_limit,plot_upper_limit)


    #label peaks in the plots
    if plot_characteristictimes != False:
        
        
        
        charl=plot_characteristictimes[0]
        timelc=plot_characteristictimes[1]
        timelc2=plot_characteristictimes[2]
        
        if plot_time_window: #remove peak labels out of retention time range to plot
            indc=[i for i in range(len(timelc)) if ( timelc[i] > plot_time_window[0] and timelc2[i] < plot_time_window[1])]
            print(indc)
            
            charl=[charl[i] for i in indc]
            timelc=[timelc[i] for i in indc]
            timelc2=[timelc2[i] for i in indc]
            
            
        for ill in range(0,len(charl)):
            for a in ax:
                a.text((timelc[ill]+timelc2[ill])/2.,.99,charl[ill],transform=matplotlib.transforms.blended_transform_factory(a.transData, a.transAxes),color='red')
                
                
    plt.tight_layout()
    fig.patch.set_facecolor((0, 0, 0, 0))
    ax[0].legend(fontsize=6)
    
    if plot_output != False:
        fig.savefig(plot_output)
    
    if show_viz == False:
        plt.close()
    return [chrom_df, peak_df, [fig, ax]]


#######
#Additional functions to analyze growth curves, plot OD versues concentration etc - not to become part of the core package
#######


def read_hplcsettings(samplenames,filename_settings,calibration_settings=None):
    #read in hplc settings from csv file and returns settings as dictonary
    #samplenames: name of single hplc run, or list of hplc runs (if several runs, setting of first one is provided and warning is given if settings for other settings vary
    #calibration_settings="glucose" searches for concentration information on glucsose
    
    #ignore sample names which are nans.
    samplenames = [x for x in samplenames if str(x) != 'nan']
    samplenames = [str(x).strip() for x in samplenames]
    if type(samplenames) != list:
        samplenames=[samplenames]
    analysis_settings=pd.read_csv(filename_settings,skiprows=1,index_col=0)
    try:
        analysis_settings=analysis_settings.loc[samplenames]
    except:
        raise TypeError(f'Samples {str(samplenames)} not found in {filename_settings}.')
        
    
    #check if settings are unique, provide warning if not.
    for uniquereq in ['foldername','peak_assignment','hplc_machine','peaks_manual','peaks_autofitting','background_substraction','setting_background_substraction','fitting_boundaries']:
        if analysis_settings[uniquereq].unique().shape[0]>1:
            print("Warning: no unique values of key settings. Assume values for first sample.")
        #display(analysis_settings[uniquereq].unique())
    #prepare and return dict with settings
    settings={}  
    try:
        rowfirst=analysis_settings.iloc[0]
    except:
        raise TypeError(f'Problem reading {str(samplenames)} from {filename_settings}. Make sure setting spreadsheet contains information for these samples.')
        
    for key in ['foldername','peak_assignment','hplc_machine','setting_background_substraction']:
        settings[key]=rowfirst[key]
    for key in ['peaks_autofitting','background_substraction']:
        settings[key]=bool(rowfirst[key])
    for key in ['fitting_boundaries','time_window']:
        curlist=rowfirst[key][1:-1].split(",")
        curlist=list(np.float_(curlist))
        settings[key]=curlist
    for key in ['peaks_manual']:
        if rowfirst[key] in [None,"None","[]"]:
            settings[key]=[]
        else:
            curlist=rowfirst[key][1:-1].split(",")
            curlist=list(np.float_(curlist))
            settings[key]=curlist
    if calibration_settings==None:
        return settings
    else:
        unitc=analysis_settings.filter(regex='unit').columns
        for u in unitc:
            settings[u]=rowfirst[unitc]
        try: 
            conc=analysis_settings.filter(regex='concentration')
        except:
            conc=analysis_settings["concentration"]
        return settings, conc


def run_calibration(foldername,substrate,retentiontime_range, filename_settings="hplc_settings_calibration.csv", analysis_window=None, peakpositions=None, list_samplenames=None, name_prefix=None, date=None, output_json='calibration.json',save_full_calibrationcurve=None, AUC_treshhold=None,plot_time_window=None, plot_upper_limit=None,plot_lower_limit=None,display_fullpeaklist=False):
    """
    Generate calibration curves from collection of samples with different concentrations of a substrate or several different substrates. 
    Parameters
    ----------
    foldername : str
        Name of folder in which hplc row chromatograms are stored. E.g. data_hplcrawdata/calibrationrun
    substrate : str or list
        Name of substrate or substrates provided in samples and for which calibration curves should be generated. E.g. substrate='glucose' or substrate=['glucose','acetate']
    retentiontime_range : list or list of list
        Min and max of retention time range for which peak should be associated with substrate(s). E.g.  retentiontime_range=[15.5,17] or retentiontime_range=[[15.5,17],[20,22]]
    filename_settings : str
        Path to file with settings for each hplc setting processed in calibration run. This file also provide the concentrations of substrates which have been added to each sample.
    list_samplenames : list or None
        List of samples to be taken. If none, all hplc samples in folder will be analyzed.
    name_prefix : str or None
        Prefix to name calibration curve. Final name of calibration curve is assembled as {name_prefix}_{substrate}_{date}.
    date : str or int
        Date to name calibration curve. Final name of calibration curve is assembled as {name_prefix}_{substrate}_{date}.
    analysis_window : list or None
        Range of retention times included in analysis. If None, values will be taken from settings table (colum time window)
    peakpositions : list or None
        List of peak positions which should be included. If None, settings will be taken form settings table (column peaks_manual)
    output_json : str
        Path to json file to store parameters of calibration curve (slope and intercept). The json file can handle several calibration curves.  
    save_full_calibrationcurve : str or None
        Folder name to save calibration curves and plots of calibration curves. If None, only parameters of calibration files will be saved in json file but not the curves themselves. 
    AUC_treshhold : Float or None
        Minimum area under the curve required to consider fit. If None, all peaks will be used.
    plot_time_windows : None or list
        Range of chromatograph to be plotted. If None, entire analyzed range will be plotted.
    plot_upper_limit : None or list. 
        Upper limit of chromatograph intesity. If None, automatic adjustment of range to show all peaks. 
    plot_lower_limit : None or list. 
        Upper limit of chromatograph intesity. If None, automatic adjustment. 
    Output to harddrive
    -------------------
    json file. Saved to harddrive
    

    """
    
    #if list of samples is not provided, take all samples in folder. 
    if list_samplenames==None:
        #get list of samples in the folder
        raw_calib_files = glob.glob(os.path.join(foldername,'*.txt')) #get all files
        if len(raw_calib_files)==0:
            raise TypeError(f'No hplc files found in {foldername}. Make sure folder name is correct and hplc raw reads are stored ending with txt.')
        
        list_samplenames = [os.path.basename(x) for x in raw_calib_files]
        list_samplenames = [x[:-4] for x in list_samplenames]
        
    #read in stored hplc settings
    hplc_setting, concentrations_calibrations=read_hplcsettings(list_samplenames,filename_settings,calibration_settings=True)
    
    #make folder for converted files if folder does not exist
    pathconverted=os.path.join(foldername,"converted")
    if os.path.exists(pathconverted)==False:
        os.makedirs(pathconverted)
    
    files=[]
    #go through each sample and look for converted file. It it does not exist, run conversion. 
    for samplename in list_samplenames:
        samplepath_converted=os.path.join(pathconverted,samplename+"_chromatogram.csv")
        if os.path.exists(samplepath_converted)==False:
            print("conversion with")
            print(glob.glob(os.path.join(foldername,'*.txt')))
            convert(glob.glob(os.path.join(foldername,'*.txt')),hplc_machine=hplc_setting["hplc_machine"])
            if os.path.exists(samplepath_converted)==False:
                print(samplepath_converted)
                conversion_needed
        else:
            pass
        files.append(samplepath_converted)
        
    
    time_window_analysis=hplc_setting["time_window"]
    if analysis_window !=None:
        time_window_analysis=analysis_window
    peaks_manualc=hplc_setting["peaks_manual"]
    if peakpositions !=None:
        peaks_manualc=peakpositions
 
    
    #plot peaks
    cal_chroms, cal_plot, cal_rettimes = batch_plot(files,time_window=time_window_analysis,plot_verticallines=peaks_manualc,backgroundsubstraction=hplc_setting["background_substraction"],backgroundsubstraction_iterations=hplc_setting["setting_background_substraction"],plot_time_window=plot_time_window,plot_upper_limit=plot_upper_limit,plot_lower_limit=plot_lower_limit,detection_peakpositions='peakautodetection')
    
    #combine detected with manual peak
    if type(peaks_manualc)==list and len(peaks_manualc)>0:
        peakpositions=add_missing_peaks(cal_rettimes,peaks_manualc)
    else:
        peakpositions=cal_rettimes 
    
    if type(peakpositions) != list or len(peakpositions)==0:
        print(peaks_manualc)
        raise TypeError(f'No peak positions define. {str(peakpositions)}. Make sure to set peak positions by hand or use peak autofitting.')
        
        
    
    #print(backgroundsubstraction_iterations+2)
    #print(chromogram_characteristics)
    #print(plot_output_chromatogram)
    #plot peak positions, including added peaks
    
    
    #prepare for saving plots
    if save_full_calibrationcurve !=None:
        name_output_simple=name_prefix+"_"+str(date) #unique name of calibration curve
        plot_outputprocess=os.path.join(save_full_calibrationcurve,"calibration"+name_output_simple+"_chromatograms.pdf")
    else:
        plot_outputprocess=False
        
    
    
    #fit peaks
    cal_chroms, cal_peaks, cal_plot = batch_process(files,time_window=time_window_analysis,backgroundsubstraction=hplc_setting["background_substraction"],backgroundsubstraction_iterations=hplc_setting["setting_background_substraction"],plot_output=plot_outputprocess,manual_peak_positions=peakpositions,plot_time_window=plot_time_window,plot_upper_limit=plot_upper_limit,plot_lower_limit=plot_lower_limit,show_viz=True)
    
    if display_fullpeaklist:
        display(cal_peaks)
    
    #generate calibration curve
    if type(substrate) != list:
        substrate=[substrate]
        retentiontime_range=[retentiontime_range]
    for iS in range(0,len(substrate)):
        sub=substrate[iS]
        retentiontimec=retentiontime_range[iS] 
        print("******Calibration curve*****")
        print("Calibration for substrate "+sub)
        print(retentiontimec)
        
        name_output=name_prefix+"_"+sub+"_"+str(date) #unique name of calibration curve
        

        #select peaks according to retention time
        calibration_peak = cal_peaks.loc[(cal_peaks['retention_time'] >= retentiontimec[0]) & (cal_peaks['retention_time'] <= retentiontimec[1])]
        print("Peaks associated with substrate")
        display(calibration_peak)
        
        concs = []
        for index, row in calibration_peak.iterrows():
            samplec=row["sample"].replace("_chromatogram","")
            try:
                concadd=concentrations_calibrations.at[samplec,"concentration_"+sub]
                if str(concadd) in ["","nan"]:
                    concadd=concentrations_calibrations.at[samplec,"concentration"]
            except:
                concadd=concentrations_calibrations.at[samplec,"concentration"]
            concs.append(concadd) 
            
        
        if len(concs)>len(list_samplenames):
            print("Warning: More than 1 peak was assigned to substrate "+sub+". Consider changing peak assignment and providing narrower time windows.")
        calibration_peak.loc[:,'concentration'] = list(concs)
        
        if AUC_treshhold != None:
            calibration_peak=calibration_peak.loc[calibration_peak['area']>AUC_treshhold]
        
        
        # Perform the linear regression and get slope and intercept
        output = scipy.stats.linregress(calibration_peak['area'], calibration_peak['concentration'])
        slope = output[0]
        intercept = output[1]
        print("Slope: "+str(slope))
        print("Intercept: "+str(intercept))

        #Save fit coordinates to file (providing unique name for calibration)
        if os.path.exists(output_json):
            with open(output_json, "r") as jsonFile:
                data = json.load(jsonFile)
            data[name_output+'_slope'] = slope
            data[name_output+'_intercept']=intercept
        else:
            data={} 
        
        try:
            unitout=hplc_setting["unit_"+sub].values[0]
            
            if str(unitout) in ["","nan"]:
                unitout=hplc_setting["unit"].values[0]
        except:
            unitout=hplc_setting["unit"].values[0]
        data[name_output+'_unit']=unitout
        with open(output_json, "w") as jsonFile:
            json.dump(data, jsonFile)
        print("Calibration curve saved as: "+output_json)

        #plot the calibration curve 
        peak_range = np.linspace(0, 1.1*calibration_peak['area'].max(),300) # Set up a range of concentrations to plot 
        fit = intercept + slope * peak_range # Compute the calibration curve for plotting
        
        fig, ax = plt.subplots(1,1,figsize=(5,3)) # its for the plots: 1 row and 3 panels in 1 row (figure size in cm)
        ax.set_title(sub)
        ax.plot(peak_range, fit, 'k-', label='fit')
        ax.plot(calibration_peak['area'],calibration_peak['concentration'],'o')
        ax.set_ylabel('concentration ['+unitout+']')
        ax.set_xlabel('AUC of fitted peak [mV min]')
        fig.tight_layout() #optimizing arangement of axes in figures  or the plots to look better in terms of layout
        plt.show()
        plotfilename="calibration"+name_output
        if save_full_calibrationcurve != None:
            fig.savefig(os.path.join(save_full_calibrationcurve,plotfilename+".pdf"))
            calibration_peak[['area','concentration']].to_csv(os.path.join(save_full_calibrationcurve,plotfilename+".csv"))

def add_missing_peaks(peakpositions,manuallist,prioritize_manual=False,mindistance=0.5): #standard 0.3
            #go through every sample and add values if missing
            
            if type(manuallist)==list and len(manuallist)>0:
                outputlist = []
                
                
                for si in range(len(peakpositions)):
                    if prioritize_manual:
                        outputlist.append(manuallist.copy())
                        majorlist=np.abs(np.array(manuallist))
                        for pk2 in peakpositions[si]:


                            curl=pk2*np.ones(majorlist.shape)
                            idx = (np.abs(majorlist-curl)).argmin()
                            diff=np.abs(curl-majorlist)[idx]
                            if diff>mindistance:
                                outputlist[si].append(pk2)  
                        outputlist[si].sort()
                    else:
                        outputlist.append(peakpositions[si].copy())
                        majorlist=np.abs(np.array(peakpositions[si]))
                        for pk2 in manuallist:
                            curl=pk2*np.ones(majorlist.shape)
                            idx = (np.abs(majorlist-curl)).argmin()
                            diff=np.abs(curl-majorlist)[idx]
                            if diff>mindistance:
                                outputlist[si].append(pk2)  
                        outputlist[si].sort()
                return outputlist      
            else:
                    return peakpositions
            
def analyze_chromatograms(list_samplenames,foldername="data_hplcrawdata/test",filename_settings="hplc_setting.csv",analysis_window=None, peakpositions=None, plot_output_chromatogram='chromatogram.pdf',plot_output_chromatogramanalysis='analysis.pdf',calibrationfile="calibration.json",chromogram_characteristics=[],peakassignment_file="peak_assignment_standard.csv",folder_peakassignment="data_peakassignment",plot_time_window=None, plot_upper_limit=None,plot_lower_limit=None,display_fullpeaklist=False):
    
    """
    Analyses series of chromatograms. 
    Parameters
    ----------
    list_samplenames : str or list
    foldername : str
        Name of folder in which hplc row chromatograms are store3d
    filename_settings : str
        path of file in which HPLC settings are stored. 
    plot_output_chromatogram : 
    
    plot_output_chromatogramanalysis :
    
    analysis_window : list or None
        Range of retention times included in analysis. If None, values will be taken from settings table (colum time window)
    peakpositions : list or None
        List of peak positions which should be included. If None, settings will be taken form settings table (column peaks_manual)
    
    
    calibrationfile:
    
    chromogram_characteristics:
    
    peakassignment_file
    
    folder_peakassignment
    
    Returns
    -------
    hplc number : json file. Saved to harddrive

    """
    
    
    

    #get settings from settings file
    hplc_settings=read_hplcsettings(list_samplenames,filename_settings)
    #get peak assignment information
    num_substrates,peaks, peaks_short, calibration, peak_peak_retentiontimes_min, peak_peak_retentiontimes_max = read_peakassignment(peakassignment_file,folder_peakassignment)
    
    #peak position
    time_window_analysis=hplc_settings["time_window"]
    if analysis_window !=None:
        time_window_analysis=analysis_window
    peaks_manualc=hplc_settings["peaks_manual"]
    if peakpositions !=None:
        peaks_manualc=peakpositions
        
        
        
        
    
    #filter out nans from sample
    list_samplenames = [x for x in list_samplenames if str(x) != 'nan']
    #check if hplc raw read files are there and converted
    files=[]
    
    #make folder for converted files if folder does not exist
    pathconverted=os.path.join(foldername,"converted")
    if os.path.exists(pathconverted)==False:
        os.makedirs(pathconverted)
    #go through each sample and look for converted file. It it does not exist, run conversion. 
    
    for samplename in list_samplenames:
        samplepath_converted=os.path.join(pathconverted,samplename+"_chromatogram.csv")
        samplepath=os.path.join(foldername,samplename+".txt")
        if os.path.exists(samplepath_converted)==False:
            convert(samplepath,hplc_machine=hplc_settings["hplc_machine"])
            
            if os.path.exists(samplepath_converted)==False:
                print(hplc_settings["hplc_machine"])
                print(samplepath)
                print(samplepath_converted)
                conversion_needed
        else:
            pass
        files.append(samplepath_converted)

    backgroundsubstraction=hplc_settings["background_substraction"]
    backgroundsubstraction_iterations=hplc_settings["setting_background_substraction"]
    manuallist=peaks_manualc
    manuallist.sort()
    
    #merge manual peaks and peaks from autofitting
    if hplc_settings["peaks_autofitting"]:
        chroms, plot, peakpositions = batch_plot(files,time_window=time_window_analysis, plot_time_window=plot_time_window, plot_upper_limit=plot_lower_limit,plot_lower_limit=plot_lower_limit, detection_peakpositions='peakautodetection',backgroundsubstraction=backgroundsubstraction,backgroundsubstraction_iterations=backgroundsubstraction_iterations)
        peakpositions=add_missing_peaks(peakpositions,manuallist)
    else:
        peakpositions=manuallist
    
    #print(backgroundsubstraction_iterations+2)
    #print(chromogram_characteristics)
    #print(plot_output_chromatogram)
    #plot peak positions, including added peaks
    
    chroms, plot, _ = batch_plot(files,time_window=time_window_analysis, plot_time_window=plot_time_window, plot_upper_limit=plot_upper_limit,plot_lower_limit=plot_lower_limit, plot_verticallines=peakpositions,backgroundsubstraction=backgroundsubstraction,backgroundsubstraction_iterations=backgroundsubstraction_iterations,plot_characteristictimes=chromogram_characteristics,plot_output=plot_output_chromatogram)
    boundpars=hplc_settings["fitting_boundaries"]
    #1. min. amplitude of peak height (multiplication of amplitude of detected peak position)
    #2. min. position of peak (offset of tected peak position to the left)
    #3. min. width of fitted peak
    #4. min. value of skew parameter
    #5. max. amplitude of peak height (multiplication of amplitude of detected peak position)
    #6. max. position of peak (offset of detected peak position to the right)
    #7. max. width of fitted peak
    #8. max. value of skew parameter

    #set options background correction
    data_chroms, data_peaks, data_plot = batch_process(files, time_window=time_window_analysis, show_viz=True,plot_time_window=plot_time_window, plot_upper_limit=plot_upper_limit,plot_lower_limit=plot_lower_limit,manual_peak_positions=peakpositions,boundpars=boundpars,buffer=100,backgroundsubstraction=backgroundsubstraction,backgroundsubstraction_iterations=backgroundsubstraction_iterations,plot_characteristictimes=chromogram_characteristics,plot_output=plot_output_chromatogramanalysis)
    
    #open calibration file
    with open(calibrationfile, "r") as jsonFile:
            cal_data = json.load(jsonFile)

    #### fix here....script needs information on peak anmes, positions, etc.
            
    #go through each sample.....
    samplelist=data_peaks['sample'].unique()
    for sample in samplelist:
        dictout={}    
        data_selection=data_peaks.loc[data_peaks['sample']==sample]
        
        if display_fullpeaklist:
            display(data_selection)

        #check for  peaks we want to look at in that particular sample (e.g. peaks=["acetate",glucose])
        #statements below are relating the calibration data with the HPLC analysis
        for il in range(0,len(peaks)):
            peak=peaks[il]
            cal_name=calibration[il]
            slope=cal_data[cal_name+"_slope"]
            intercept=cal_data[cal_name+"_intercept"]
            peak_retentiontime_min=peak_peak_retentiontimes_min[il]
            peak_retentiontime_max=peak_peak_retentiontimes_max[il]

            data_selection_currentpeak=data_selection[(data_selection['retention_time'] >= peak_retentiontime_min) & (data_selection['retention_time'] <= peak_retentiontime_max)]
            #display(data_selection_currentpeak)
            if data_selection_currentpeak.shape[0]==0:
                print("Peak "+peak+" not found in sample "+sample)
            elif data_selection_currentpeak.shape[0]>0:
                dictout[peak]=data_selection_currentpeak['area'].values[0]*slope+intercept
                dictout[peak+"_cal"]=cal_name
                dictout['sample']=sample

            elif data_selection_currentpeak.shape[0]>1:
                morethanonepeak #check why you get more than one peak within the retention time range
        #display(dictout)
        with open('data_chromatogramanalysis/'+sample+'.json', 'w') as fp:
            json.dump(dictout, fp)
        # Isolate the peak used for calibration
    #print(data_chroms)
    return(data_chroms)


def read_peakassignment(filename_peakassignment,folder_peakassignment="data_peakassignment"):
    if len(filename_peakassignment)<5 or filename_peakassignment[-4:] != ".csv":
        peakinfo=pd.read_csv(os.path.join(folder_peakassignment,filename_peakassignment)+".csv")
    else:
        peakinfo=pd.read_csv(os.path.join(folder_peakassignment,filename_peakassignment))
    num_substrates=peakinfo["substrate"].shape[0]
    return num_substrates, peakinfo["substrate"].tolist(), peakinfo["substrate_label"].tolist(), peakinfo["calibration"].tolist(), peakinfo["retention_time_min"].tolist(), peakinfo["retention_time_max"].tolist()
    
def process_growth_curve(file_growthdata,folder_hplcdata=None,folder_growthdata="data_growthcurves",culture_list="All",filename_settings="hplc_settings.csv",folder_peakassignment="data_peakassignment", analysis_window=None, peakpositions=None, folder_output="output_growthanalysis",plot_time_window=None, plot_upper_limit=None,plot_lower_limit=None,calibrationfile="calibration.json",display_fullpeaklist=False,peak_assignment=None):
    
    """
    Process growth curve and associated hplc samples. Plot substrat vs OD figures. Store obtained trends for further processing. 
    
    Parameters
    ----------
    list_samplenames : str or list
    foldername : str
        Name of folder in which hplc row chromatograms are stored
    filename_settings : str
        path of file in which HPLC settings are stored. 
    culture_list : str or list
        
    plot_output_chromatogram : 
        If "All" all growth curves in csv file will be analyzed. Otherwise, provide list with specific names of hplc runs.
    folder_hplcdata : str or None
        Folder in which hplc runs are store. If None, it is assumed data are in a subfolder of "data_hplcrawdata" named file_growthdata
    calibrationfile: str
        Path of json file with calibration information.
    folder_peakassignment: str
        Path of folder in which peak assignments are stored.z
    plot_time_window=None, plot_upper_limit=None,plot_lower_limit=None
    
    analysis_window : list or None
        Range of retention times included in analysis. If None, values will be taken from settings table (colum time window)
    peakpositions : list or None
        List of peak positions which should be included. If None, settings will be taken form settings table (column peaks_manual)
    peak_assignment : str or None
        Name of peak assignment file to use. If None, values will be taken from settings table (colum peak_assignment)
        
    
    Returns
    -------
    hplc number : json file. Saved to harddrive

    """
    
    #generate folder if needed
    if not os.path.exists(folder_output):
        os.makedirs(folder_output)
    
    grdata=pd.read_csv(os.path.join(folder_growthdata,file_growthdata+".csv"))
    
    grdata.dropna(axis=1, how='all',inplace=True)
    growthrates=pd.DataFrame() #prepare dataframe for results

    #read in all samples in growth file
    namesamplelistc=[]
    for il in range(0,int(grdata.shape[1]/5)):
        namesamplelistc.append(grdata.columns[5*il+2])
    print("All cultures in growth rate file "+os.path.join(folder_growthdata,file_growthdata+".csv")+":")
    print(namesamplelistc)
    
    strainnames=[]
    for il in range(0,int(grdata.shape[1]/5)):
        namesamplec=grdata.columns[5*il+2]
        if culture_list=="All" or namesamplec in culture_list:
            
            
            #read in growth data
            time=grdata.iloc[:,5*il+1].to_numpy().flatten()
            od=grdata.iloc[:,5*il+2].to_numpy().flatten()
            od2=grdata.iloc[:,5*il+3].to_numpy().flatten()
            hplcname=grdata.iloc[:,5*il+4].to_numpy().flatten() #list of hplcsamples as named in excel column
            hplcdata=hplcname.copy()
          
            #read in hplc settings
            hplc_settings=read_hplcsettings(hplcname,filename_settings)
            
            #get from these settings which peak assignment file to use and read in peak assignment information
            if peak_assignment==None:
                peakassignmentc=hplc_settings['peak_assignment']
            else:
                peakassignmentc=peak_assignment
            num_substrates,peaks, peaks_short, calibration, peak_peak_retentiontimes_min, peak_peak_retentiontimes_max = read_peakassignment(peakassignmentc,folder_peakassignment=folder_peakassignment)
            
            
           
            chromogram_characteristics=[peaks_short,peak_peak_retentiontimes_min,peak_peak_retentiontimes_max]
            
            #prepare plots
            fig, ax = plt.subplots(1,num_substrates+1,figsize=((num_substrates+1)*4,3)) # its for the plots: 1 row and 3 panels in 1 row (figure size in cm)
            substrate_ods = [[] for x in range(num_substrates)]
            substrate_values = [[] for x in range(num_substrates)]
            
            #if folder name is not provided use same name as filename with growth rate data
            if folder_hplcdata == None:
                folder_hplcdata=os.path.join("data_hplcrawdata",file_growthdata)
            #### run analysis of chromatograms
            
            curchrom=analyze_chromatograms(hplcname,folder_hplcdata,analysis_window=analysis_window, peakpositions=peakpositions, plot_output_chromatogram=os.path.join(folder_output,namesamplec+"_"+file_growthdata+'_chromatogram.pdf'),plot_output_chromatogramanalysis=os.path.join(folder_output,'chromatogram_analysis_'+namesamplec+"_"+file_growthdata+'.pdf'),filename_settings=filename_settings,chromogram_characteristics=chromogram_characteristics,peakassignment_file=peakassignmentc,calibrationfile=calibrationfile,plot_time_window=plot_time_window, plot_upper_limit=plot_upper_limit,plot_lower_limit=plot_lower_limit,display_fullpeaklist=display_fullpeaklist)
            
            
            
            #######################################
            #### plot growth curve and analyze trends in growth
            #######################################
            for il3 in range(hplcname.shape[0]):
                curname=hplcname[il3] # name from the excel file we give the name and od (below), as we are interested in these two data for hplc analysis and even if we dont have anything in the hplc column, the script runs as it is for growth curves without showing error
                curod=od[il3]
                curtime=time[il3]
                #check if hplc file is defined and can be openend
                if str(curname) != "nan":
                    with open(os.path.join('data_chromatogramanalysis',curname+'_chromatogram.json'), "r") as jsonFile:
                        curconcentrations = json.load(jsonFile)
                        sC=-1
                        for substr in peaks:
                            sC=sC+1
                            try: 
                                substrate_values[sC].append(curconcentrations[substr])
                                substrate_ods[sC].append(curod)
                                #substrate_times[sC].append(curtime)
                            except:
                                print("Warning: substrate information for "+substr+" not found")
            #fit growth rate
            #use only data without nan values to do fit
            indd=np.argwhere(np.isfinite(od2))[:,0] #list of all rows where od column has a finite value
            timefit=time[indd]
            odfit=od2[indd] 
            fitc=np.polyfit(timefit,np.log(odfit),1) #what is np.polyfit doing. What is np.log doing? polyfit fits the data to a polynomial curve. And np.log ensures that we are transforming the y-axis data to support a linear fit on the final semilog plot. In general, these commands are

            #trying to fit the observed exponential growth curve to an equation
            namec=grdata.columns[5*il+2]
            strainnames.append(namec.split("_")[0])
            timemax=np.nanmax(time)
            odmax=np.nanmax(od)
            odmin=np.nanmin(od)

            #plot growth curve
            ax[0].plot(time,od,color='b',marker='o',ls='',alpha=0.3) #uncorrected OD values
            ax[0].plot(time,od2,color='b',marker='o',ls='') #OD data used
            #plot exponential fit to growth curve
            trange=np.linspace(0,4,100)
            ax[0].plot(trange,np.exp(fitc[1])*np.exp(fitc[0]*trange),ls='-',color='b',alpha=0.5)
            ax[0].set_title(strainnames[-1]+"; "+str(round(fitc[0],2))+" 1/h")
            #set time range
            ax[0].set_xlim([0,1.1*timemax])
            ax[0].set_ylim([0.9*odmin,1.1*odmax])
            #set log scale and labels
            ax[0].set_yscale("log")
            ax[0].set_xlabel("time (h)")
            ax[0].set_ylabel("density (OD)")
            ax[0].set_yticks([0.04,0.08,0.16,0.32,0.64])
            ax[0].set_yticklabels([0.04,0.08,0.16,0.32,0.64])

            fit_substrate=[]
            fit_substrate_error=[]
            #generate substrate concentration vs OD plots and get slopes
            for sC in range(0,num_substrates): #go through all substrates
                #plot obtained values for each substrate
                ax[sC+1].plot(substrate_ods[sC],substrate_values[sC],marker='o',ls='',color='r')
                ax[sC+1].set_xlabel("density (OD)")
                ax[sC+1].set_ylabel(peaks[sC]+" (mM)")
                try:
                    curfit,curfiterror=np.polyfit(substrate_ods[sC],substrate_values[sC],1, cov=True)
                    fit_substrate.append(curfit)
                    fit_substrate_error.append(curfiterror)
                    odrange=np.linspace(0.1,0.5,200)
                    ax[sC+1].plot(odrange,odrange*fit_substrate[-1][0]+fit_substrate[-1][1],ls='--',alpha=0.5,color='r')
                    ax[sC+1].set_title(str(round(fit_substrate[-1][0],2))+" mM/OD")
                except:
                    fit_substrate.append([np.nan,np.nan])
                    fit_substrate_error.append([[np.nan]])

            dictc={'experiment': file_growthdata, 'culture':namec, 'growth-rate':fitc[0],'y0':np.exp(fitc[1])}
            for sC in range(0,num_substrates):
                          dictc[peaks[sC]+"-turnover"]=fit_substrate[sC][0]
                          dictc[peaks[sC]+"_error"]=fit_substrate_error[sC][0][0]
                          dictc[peaks[sC]+"_val"]=str([substrate_ods[sC],substrate_values[sC]])     
                            
            growthrates = growthrates.append(dictc, ignore_index=True)
            
            
            
            

            fig.tight_layout() #optimizing arangement of axes in figures  or the plots to look better in terms of layout
            plt.show()
            fig.savefig(os.path.join(folder_output,strainnames[-1]+"_"+file_growthdata+"_growthplots.pdf"))

    substrates=peaks
    listorder=["culture","growth-rate","y0"]

    for sub in substrates:
        listorder.append(sub+"-turnover")
        listorder.append(sub+"_error")
        listorder.append(sub+"_val")

    growthrates = growthrates[listorder]
   
    filename_output=os.path.join(folder_output,'analysis_'+file_growthdata+'.csv')
    growthrates.to_csv(filename_output)
    return(growthrates)


def generate_setting_template(foldername,output="hplc_settings_template.csv",hplc_machine="Shimadzu_ProminenceLC2030C"):
    #to start adding settings for samples to list of hplc analysis settings, you can use the following:
    raw_calib_files = glob.glob(os.path.join(foldername,'*.txt')) #this gives list of the files
    samplename = []
    for i in raw_calib_files:
        try:
            samplename.append(os.path.basename(i)[:-4])
        except:
            pass
    dfout=pd.DataFrame(samplename,columns=["sample"])
    dfout["foldername"]=os.path.basename(foldername)
    dfout["peak_assignment"]=""
    dfout["hplc_machine"]=hplc_machine
    dfout["peaks_manual"]="[]"
    dfout["peaks_autofitting"]="True"
    dfout["background_substraction"]="True"
    dfout["setting_background_substraction"]="80"
    dfout["fitting_boundaries"]='[0.1,0.3,0.0,-0.1,10,0.3,0.8,0.1]'
    dfout["time_window"]='[12,40]'
    dfout.sort_values(by=['sample'],inplace=True)
    display(dfout)
    dfout.to_csv(output,index=False)
    display("Adjust "+output+" and copy rows to main settings table, e.g. hplc_settings.csv")
