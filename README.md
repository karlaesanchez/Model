# Model

The file GenM_L_101018.m is the code for the current vascular model. The file name stands for: Generational Main Long. This is the one I use for plots, hence why the function output is long as I need several variables for plotting. 

There are variations of this function for fewer outputs that feed into the compartmental model, and another one that feeds into the autoregulation function. All have the same lines of code but different outputs. The one feeding the autoregulation function has the following parameters suppressed: ica_fac, va_ratio, N, gamma, lambda. That is since the autoregulation function needs to find the ideal combination of parameters.


The file ARGLcombM_101018.m is the code for finding the appropriate set of combinations for ica_fac, va_ratio, N, gamma, lambda. The file name stands for: Autoregulation Gamma Lamba Main. As we have set ica_fac, va_ratio, N, lambda, in previous analysis, these are all suppressed and only the changing gamma is not. This file keeps changing specifically on the acceptance criteria (if pA and pV, etc.), but structurally is the same as here.

The file RegionsM_0918.m is the code for finding the acceptable region of gamma and lambda. This is separate from the combination code above to reduce code time as I can store the resulted matrix of acceptable gamma and lambda and then run against this matrix on the combination code.
