%% Calibration Toolkit %%
For broadband scientific echosounders using standard targets
(single beam and 4-quadrand split beam)

Written by Rachel Kahn

First, make sure you have the EK80 processing package in your path. For some
reason this code has trouble reading the filter info in pulse compressed
data files without it. 

This code reads in pulse compressed .mat files containing EK80 data.
If you have not unpacked your data, please do so before proceeding.

This package contains four 'main' files that run the calibration:
1) 'main_cal' - for 4-quadrant split-beam transducer and 1 single data file
    with all the pulse compressed calibration data.
2) 'main_multi_cal' - for 4-quadrant split-beam transducer and more than one
    data file to read in.
3) 'main_cal_singlebeam' - for single-beam transducer and 1 single data file
    to read in.
4) 'main_multi_singlebeam' - for single-beam transducer and more than one
    data file to read in.

If using more than one transducer for a single experiment, I recommend creating
a copy of the correct 'main' file for each individual transducer so that you 
can quickly rerun the calibration without having to change all the settings.

Be sure to set all the correct channel info, nominal center frequency, 
environmental parameters, standard target radius and depth, etc.
These are near the beginning of the main script.
The 'offset' is arbitrary and just for plotting purposes. 
You may tweak the inputs in the 'calcTS' function if you wish.

The code creates a 'cal' data structure with all the relevant info and data
for the calibration, and a 'gain' script with the cal curve.
cal.G is the raw gain curve and might be noisy, cal.Gsmooth is the smoothed
gain curve that you can apply to data.
