# GCC-large-FFT-size-by-part
Computation of GCC on a very large length buffer by using small frames.

<p>The code is avalaible in MATLAB and C++.</p>
<p>It is based on https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/66959/versions/2/previews/fft_algorithm.m/index.html.</p>
<p>The Cooley-Tukey FFT is slightly modified to work frame per frame (both side FFT and inverse FFT) in order to compute the Generalized Cross Correlation (GCC).
I have also implemented different weight for the GCC : PHAT, CSPm and rhoCSP.</p>
<p>These weights are described in this paper:
Marinescu, R. S., Buzo, A., Cucu, H., & Burileanu, C. (2013). Applying the accumulation of cross-power spectrum technique for traditional generalized cross-correlation time delay estimation. International Journal On Advances in Telecommunications–IARIA, submitted invited paper.</p>

<p>For C++, it requires Eigen library. The include directory must be updated in the project.</p>
