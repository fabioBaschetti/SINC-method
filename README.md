# SINC method [source code from - The SINC way: A fast and accurate approach to Fourier pricing]

The repository contains everything you need for pricing European options and their digital components with the SINC approach and other standard Fuorier methods from the literature: Carr-Madan (1999), Lewis (2000) and Fang-Oosterlee (2008). FT\FFT\frFFT versions are implemented for all those methods - when they exist - and a number of tests carried out to assess the relative performance of the SINC as opposed to its standard competitors. GBM, Heston, CGMY and the rough Heston model (forward variance form) are employed for such purposes.

Matlab functions SINC_discFT\SINC_fastFT\SINC_fracFT are commented in full detail so that the user can easily jump from the code to the text, and viceversa. Similarly, the scripts refer to the figures and tables in the paper that they reproduce:
- main_convMODEL.m studies the convergence of the FT version of both SINC and COS to the 'true' option price under MODEL. When MODEL=rHeston this should be read with main_convrHeston_fig.m, where a pictorial representation for the convergence of the methods is given
- main_drawPDF.m   plots the densities in the text (for CGMY and rHeston) {other models can be easily implemented} (this uses function chf2pdf to reconstruct the pdf from the characteristic function via the FFT) 
- main_FFTdate.m   prices an entire volatility surface with the FFT\frFFT versions of SINC\Carr-Madan\Lewis

We also upload the surfaces we use and .txt files where every maturity is given the corresponding value of Xc for determining the truncation range.

The content of any other function\script should be clear from the title when following the logic described above.
