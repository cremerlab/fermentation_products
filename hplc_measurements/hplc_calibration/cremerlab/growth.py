import numpy as np 
import pandas as pd
from .bayes import *
import tqdm
import statsmodels.tools.numdiff as smnd
import matplotlib.pyplot as plt 
import seaborn as sns
import scipy.optimize
import scipy.stats
sns.set()

def infer_growth_rate(data, 
                      od_bounds=None, 
                      convert_time=True, 
                      groupby=None,
                      viz=False,
                      cols={'time':'clock_time', 'od':'od_600nm'}, 
                      return_opts=False,
                      print_params=True,
                      **kwargs):                  
    """
    Infers the maximal a posteriori (MAP) parameter set for the steady state growth 
    rate given measurements of optical density. This is performed via optimization 
    by minimization.

    Parameters
    ----------
    data : pandas DataFrame
        A tidy long-form pandas DataFrame with columns corresponding to the 
        measurement time and optical density.
    od_bounds : list of floats
        The lower and upper bounds of the optical density range to be considered.
        The default bounds assumed are [0.04, 0.41] inclusive.
    convert_time : bool
        If `True`, the provided time column needs to be converted to elapsed 
        time. In this case, the provided time column is assumed to be 
        the clock time of the measurement and is converted to minutes elapsed.
    groupby : list of str, optional
        The column names for the groupby operation to operate upon. For example,
        if there are multiple strains measured in the data set, a groupby of 
        `['strain']` will yield a growth rate estimate for each strain in 
        the data. A groupby of `['strain', 'replicate']` will return a growth 
        rate estimate for each strain and biological replicate.
    viz : bool
        If `True`, a vizualization of the inference will be produced and a 
        figure and axis object will also be returned.
    cols : dict, keys 'time', and 'od'
        The column names of the time and optical density measurements in the 
        DataFrame. 
    return_opts : bool
        If `True`, the approximated covariance matrix, optimal parameters, and 
        approximate hessian matrix for each grouping is returned as a dictionary.
    print_params : bool 
        If `True`, the estimated parameters will be printed to the screen when 
        the estimation is finished. 
    
    Returns
    -------
    data_df : pandas DataFrame
        A pandas DataFrame with the converted time measurements cropped to the 
        provided optical density bounds.
    param_df : pandas DataFrame
        A pandas DataFrame containing the parameters, values, and their 95% credible 
        intervals for each obejct in the provided groupby.
    opts : dict 
        If `return_opts = True`, the estimated covariance matrix, optimal parameters, 
        and approximate Hessian matrix is returned. 
    fig : matplotlib.figure.Figure
        Matplotlib figure object for the inference. Returned only if `viz = True`.
    ax :  matplotlib AxesSubplot
        The axes of the matplotlib figure. Returned only if `viz = True`.

    Notes
    -----
    This function infers the "maximal a posteriori" parameter set using a 
    Bayesian definition of probability. This function  calls the posterior 
    defined by `cremerlab.bayes.steady_state_log_posterior` which contains 
    more information about the construction of the statistical model.
    """

    # TODO: Include type checks

    if (groupby is not None) & (type(groupby) is not list):
        groupby = [groupby]

    # Unpack the time col
    time_col = cols['time'] 
    od_col = cols['od']

    # Determine the OD bounds to consider
    if od_bounds is not None:
        data = data[(data[od_col] >= od_bounds[0]) & (data[od_col] <= od_bounds[1])]

    faux_groupby = False
    if groupby is None:            
        faux_groupby = True
        data['__group_idx'] = 1   
        groupby=['__group_idx']
        iterator = data.groupby(groupby)    
    else:
        iterator = tqdm.tqdm(data.groupby(groupby), desc='Estimating parameters...')

    # Iterate through each grouping
    data_dfs = []
    param_dfs = []
    opts = {'groupby':groupby}
    iter = 0  # Iterator for opts
    output = """\n
============================================================
Parameter Estimate Summary
============================================================
\n
"""

    for g, d in iterator:
        
        # Convert time if necessary
        if convert_time:
            d[time_col] = pd.to_datetime(d[time_col])
            d.sort_values(by=time_col, inplace=True)
            d['elapsed_time_hr'] = d[time_col].values - d[time_col].values[0]
            d['elapsed_time_hr'] = (d['elapsed_time_hr'].astype('timedelta64[m]')
                                    )/60
            _time = d['elapsed_time_hr'].values
            _od = d[od_col].values
        else:
            _time  = d[time_col].values
            _od = d[od_col].values

        # Define the arguments and initial guesses of the parameters
        # lam_guess = np.mean(np.diff(np.log(_od)) / np.diff(_time))
        params = [1, _od.min(), 0.1]
        args = (_time, _od)
    
        # Compute the MAP
        res = scipy.optimize.minimize(steady_state_growth_rate_log_posterior,
                                      params, args=args, method="powell")
        # Get the optimal parameters
        popt = res.x

        # Compute the Hessian and covariance matrix
        hes = smnd.approx_hess(popt, steady_state_growth_rate_log_posterior, args=args)
        cov = np.linalg.inv(hes)

        # Extract the MAP parameters and CI
        lam_MAP, od_init_MAP, sigma_MAP = popt
        lam_CI =  1.96 * np.sqrt(cov[0, 0]) 
        od_init_CI =  1.96 * np.sqrt(cov[1, 1]) 
        sigma_CI =  1.96 * np.sqrt(cov[2, 2]) 

        if print_params:
            if faux_groupby == False: 
                header =  f"""Parameter Estimates for grouping {groupby}: {g}
------------------------------------------------------------
""" 
            else:
                header =  """Parameter Estimates
------------------------------------------------------------
"""
            output +=  header + f"""growth rate,  λ = {lam_MAP:0.2f} ± {lam_CI:0.3f} [per unit time]
initial OD, OD_0 = {od_init_MAP:0.2f} ± {lam_CI:0.3f} [a.u.]
homoscedastic error, σ = {sigma_MAP:0.2f} ± {sigma_CI:0.3f} [a.u.]
\n
"""
        # Assemble the data dataframe
        _data_df = pd.DataFrame([])
        if convert_time:
            _data_df['elapsed_time_hr'] = _time
        else:
            _data_df[time_col] = _time
        _data_df[od_col] = _od

        # Include other columns that were not carried through
        colnames = [k for k in d.keys() if k not in [time_col, od_col]]
        for c in colnames:
            _data_df[c] = d[c].values
        if '__group_idx' in _data_df.keys():
            _data_df.drop(columns=['__group_idx'], inplace=True)
        _data_df.rename(columns={'od':od_col})
    
        # Assemble the parameter dataframe
        _param_df = pd.DataFrame([])
        for title, MAP, CI in zip(['growth_rate', 'od_init', 'sigma'],
                                  [lam_MAP, od_init_MAP, sigma_MAP],
                                  [lam_CI, od_init_CI, sigma_CI]): 
            _param_df = _param_df.append({'parameter':title,
                                          'map_val': MAP,
                                          'cred_int': CI},
                                          ignore_index=True)

        # Add grouping identifier if provided
        if groupby is not None:
            if type(g) is not list:
                _g = [g]
            for title, value in zip(groupby, _g):
                _data_df[title] = value 
                _param_df[title] = value
        
        #  Append the dataframes to the storage lists
        param_dfs.append(_param_df)
        data_dfs.append(_data_df)
        opts[iter] = {'groupby': g, 'cov':cov, 'popt':popt, 'hessian':hes}  
        iter += 1

    # Concatenate the dataframes for return
    if len(data_dfs) == 1:
        data_df = data_dfs[0]
        param_df = param_dfs[0]
    else:
        data_df = pd.concat(data_dfs, sort=False)
        param_df = pd.concat(param_dfs, sort=False)
    
    if print_params:
        print(output)

    if return_opts:
        return_obj = [data_df, param_df, opts]
    else:
        return_obj = [data_df, param_df]

    if viz:
        print('Generating plot...')
        fig, ax = viz_growth_rate_inference(data_df, opts, **kwargs)
        return_obj.append([fig, ax])

    return return_obj


def viz_growth_rate_inference(data, 
                             opts, 
                             levels=None, 
                             cols={'time':'elapsed_time_hr', 'od':'od_600nm'},
                             time_units='hr',
                             show_viz=True):
    """
    Generates a vizualization of the best-fit for the growth rate as well as 
    the contours of the approximate marginal posterior for the growth rate λ and 
    initial od.

    Parameters
    ----------
    data : pandas DataFrame
        The data from which the parameters are inferred. This is generally the 
        output from `cremerlab.growth.infer_growth_rate_MAP`.
    opts : dict 
        A dictionary of the output from the inference with the covariance matrix,
        optimal parameters, and approximate Hessian for each inferred dataset. 
    levels : list, optional
        The levels to show in the contour plot. If `None`, default are the limits 
        of the 1-σ, 2-σ, 3-σ, 4-σ levels of a 2D gaussian distribution.
    cols : dict, optional
        The column names for the time (`time`) and optical density (`od`) dimensions 
        of the data. 
    time_unts : str, optional
        The units of the time dimension to be shown on the resulting plots
    show_viz : bool
        If `True`, the plot will be displayed to the screen. Default is True.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Matplotlib figure object for the chromatograms
    ax :  matplotlib AxesSubplot
        The axes of the matplotlib figure

    Notes
    -----
    The countour plots generated by this function show only the Normal approximation 
    of the posterior around the maximal a posterior parameter set. This is 
    typically OK, but may not accurately reflect the whole posterior. A more 
    formal analysis of the growth rate should use Markov Chain Monte Carlo to 
    sample the posterior distribution, but this is outside the scope of the 
    function.
    """
    sns.set()

    # Define the columns for easy referencing
    time_col = cols['time']
    od_col = cols['od']

    #TODO: type checks
    def _compute_POST_norm(LAM, OD0, lam_range, od_init_range, sig_range, popt, cov):
        LAM, OD0, SIG = np.meshgrid(lam_range, od_init_range, sig_range)

        # Marginalize over sigma
        LAM, OD0 = LAM[:, :, 0], OD0[:, :, 0]
        _popt, _cov = popt[:-1], cov[:-1, :-1]

        # Compute the posterior
        POST = scipy.stats.multivariate_normal.pdf(
                    np.dstack((LAM, OD0)), _popt, _cov)
        return POST / POST.max()

    # Get the number of unique groups.
    if opts['groupby'][0] == '__group_idx':
        n_groups = 1
    else:
        n_groups = data.groupby(opts['groupby']).ngroups
    
    # Set up the figure canvas.
    fig, ax = plt.subplots(n_groups, 2, figsize=(7, n_groups * 3))

    for i in range(n_groups):
        if n_groups == 1:
            _ax0 = ax[0]
            _ax1 = ax[1]
        else:
            _ax0 = ax[i, 0]
            _ax1 = ax[i, 1]
        for a in [_ax0, _ax1]:
            a.xaxis.set_tick_params(labelsize=8)
            a.yaxis.set_tick_params(labelsize=8)
            
        # Set tick labels
        _ax0.set_xlabel(f'elapsed time [{time_units}]', fontsize=8)
        _ax0.set_ylabel(f'optical density [a.u.]', fontsize=8)
        _ax0.set_yscale('log')  
        _ax1.set_ylabel(f'initial OD ($OD_0$)', fontsize=8)
        _ax1.set_xlabel(f'growth rate (λ) [inv {time_units}', fontsize=8)

    # Set up the iteration to account for different groupings
    if n_groups == 1:
        data['__group_idx'] = 1 
        groupby = ['__group_idx']
    else:
        groupby = opts['groupby']
    iterator = data.groupby(groupby)


    # Iterate and generate the plots
    iter = 0
    for g, d in iterator:
        _od = d[od_col].values
        _time = d[time_col].values

        # Get the parameters
        popt = opts[iter]['popt'] 
        cov = opts[iter]['cov']
        lam_MAP, od_init_MAP = popt[:2]
        lam_CI =  1.96 * np.sqrt(cov[0, 0]) 
        od_init_CI =  1.96 * np.sqrt(cov[1, 1]) 

        # Define the fit
        t_range = np.linspace(0, _time.max() + 0.2, 5)
        fit = od_init_MAP * np.exp(lam_MAP * t_range)

        # Plot the best fit
        if n_groups == 1:
            _ax0 = ax[0]
            _ax1 = ax[1]
        else:
            _ax0 = ax[iter, 0]
            _ax1 = ax[iter, 1]
        _ax0.plot(t_range, fit, '-', lw=2, color='dodgerblue', alpha=0.6,
                        label='best fit')

        # Plot the data
        _ax0.plot(d[time_col], _od, 'o', color='dodgerblue', ms=5, 
                        label='data', markeredgewidth=0.5, markeredgecolor='w')

        # Growth rate MAP neighborhood
        lam_min = lam_MAP - 1.25 * lam_CI
        lam_max = lam_MAP + 1.25 * lam_CI
        lam_range = np.linspace(lam_min, lam_max, 300)

        # Initial OD MAP neighborhood
        od_init_min = od_init_MAP - 1.25 * od_init_CI 
        od_init_max = od_init_MAP + 1.25 * od_init_CI
        od_init_range = np.linspace(od_init_min, od_init_max, 300)
        sig_range = np.linspace(0, 1, 300)

        # Compute the marginalized posterior
        LAM, OD0, _ = np.meshgrid(lam_range, od_init_range, sig_range)
        POST_norm = _compute_POST_norm(LAM, OD0, lam_range, od_init_range, 
                                       sig_range, popt, cov)

        # plot the contours
        if levels is None:
            levels = [0]
            _levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)
            for l in _levels:
                levels.append(l)

        _ax1.contourf(lam_range, od_init_range, POST_norm, levels=levels,
                            cmap='Blues') 
        _ax1.contour(lam_range, od_init_range, POST_norm, levels=levels,
                            colors='k', linewidths=0.75, linestyles='-')

        # Plot the MAP position
        _ax1.plot(lam_MAP, od_init_MAP, 'o', color='crimson', ms=5, 
                        markeredgewidth=0.5, markeredgecolor='w', zorder=1000,
                        label='MAP value')

        # Apply special tick formatting
        lower = np.round(_od.min(), decimals=2) - 0.01
        upper = np.round(_od.max(), decimals=2) + 0.01
        if lower <= 0:
            lower = _od.min()
        y_ticks = np.arange(lower, upper, 0.06)

        # Reset the ticks and labels 
        _ax0.set_yticks(y_ticks)
        _ax0.set_yticklabels([str(np.round(y, decimals=2)) for y in y_ticks])

        # Add titles as necessary
        _ax0.set_title(
           f"growth curve\n grouping {opts['groupby']} = {opts[iter]['groupby']}\n λ = {lam_MAP:0.3f} ± {lam_CI:0.3f} [inv {time_units}]",
                    fontsize=8)
        _ax1.set_title(
           f"posterior distributions\ngrouping {opts['groupby']}  = {opts[iter]['groupby']}",
                    fontsize=8)
        iter += 1

    # Add legends, tidy, and return
    if n_groups == 1:
        _ax0 = ax[0]
        _ax1 = ax[1]
    else:
        _ax0 = ax[0, 0]
        _ax1 = ax[0, 1]
    _ax0.legend(fontsize=8) 
    _ax1.legend(fontsize=8) 
    plt.tight_layout()

    if show_viz is not True:
        plt.close()

    return [fig, ax]

    
    





