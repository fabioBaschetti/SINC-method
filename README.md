# SINC-method [FT and FFT] Option pricing with the SINC approach: experiments under the rough Heston model

The repository contains everything you need for pricing European options and their digital components with the SINC approach. 
Matlab script 'main.m' performs the pricing exercise with both the FT and the FFT versions of the method in the rough Heston model 
(forward variance specification). It therefore uses functions
- SINC_discFT to price one single strike with full precision (this is actually vectorized and allows dealing with multiple strikes at 
  the same time)
- SINC_fastFT to price one entire smile by exploiting the computational power of the FFT algorithm.

On the other hand, 'cmp.m' compares two different strategies for FT pricing that correspond to functions
- SINC_discFT    : uses N evaluations of the characteristic function, always
- SINC_discFT_sbs: uses as many evaluations of the CF as they are needed to reach satisfactory accuracy (and up to N).

This impacts on CPU time. You can change N in the script and study their behavior.
Both the FT and the FFT versions of the SINC approach require the PDF of the asset log-return is truncated. Function 'truncMeasures.m' 
serves this purpose. The truncation rule depends on the cumulants and a multiplier L which we set at 100 to ensure maximum precision. 
You can also change it if you are happy with lower accuracy. Reducing the multiplier clearly boosts convergence (i.e. you will 
typically need much smaller N).
The experiments are under the rough Heston model. 'dh_pade33_coeff.m', 'dh_pade33.m' and 'phirHeston.m' are needed to compute its 
characteristic function.
You can obviously change the call to the rough Heston CF with any other model the characteristic function of which is formally known, 
and use SINC formulas with those models as well. Recall that they require a factor 2*\pi to be included in the definition of the CF 
for a martingale.
