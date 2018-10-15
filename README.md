# Model

The file GenM_L_101018.m is the matlab file for the current vascular model. The file name stands for: Generational Main Long. This is the one I use for plots, hence why the function output is long as I need several variables for plotting. 

There are variations of this function for fewer outputs that feed into the compartmental model, and another one that feeds into the autoregulation function. All have the same lines of code but different outputs. The one feeding the autoregulation function has the following parameters suppressed: ica_fac, va_ratio, N, gamma, lambda. That is since the autoregulation function needs to find the ideal combination of parameters.
