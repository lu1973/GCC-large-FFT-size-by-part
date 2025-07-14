## -*- texinfo -*-
## @deftypefn  {} {[@var{prm}, @var{states}] =} longGCC_init(@var{fs}, @var{bufferLength}, @var{bufferCoeff}, @var{smoothingFactor}, @var{fLow}, @var{fHigh})
## Initialiaze the parameters @var{prm} and the internal states @var{states} for the long GCC computation.
## The sampling frequency of the signal in Hz is given by @var{fs}. The signal can be masked in the frequency domain with @var{fLow} and @var{fHigh} (in Hz).
## The length of the innovation buffer is @var{bufferLength} and the large FFT will be compute following @var{bufferCoeff} accumulation of input buffers.
## Thus the FFT size will be @code{@var{bufferLength} * @var{bufferCoeff}}.
## The @var{smoothingFactor} will impact the reactivity of the system to detect changes in the delay (the higher value, the longer it will take to commute to the new value).
##
## Examples:
## @example
## @group
## [prm, states] = longGCC_init(16000, 300, 6000, 512, 42);
## @end group
## @end example
##
## @end deftypefn

function [prm, states] = longGCC_init(fs, bufferLength, bufferCoeff, smoothingFactor, fLow, fHigh)

    prm.buffer = bufferLength;          % buffer size in samples
    prm.fs = fs;                        % sampling frequency in Hz

    if nargin==5 || nargin<3
        error('longGCC_init: Invalid number of arguments')
    endif

    if nargin<=3
        prm.lambda = 0.8;       % smoothing factor
    else
        prm.lambda = smoothingFactor;       % smoothing factor
    endif

    if nargin<=4
        prm.f0 = 0;                     % lower bound of the frequency mask in Hz
        prm.f1 = fs/2;                  % higher bound of the frequency mask in Hz
    else
        prm.f0 = fLow;                  % lower bound of the frequency mask in Hz
        prm.f1 = fHigh;                 % higher bound of the frequency mask in Hz
    endif

    prm.N = bufferCoeff*prm.buffer;     % length of the large FFT (multiple of the buffer length)
    prm.N1 = prm.buffer;                % buffer length

    prm.maxlag = 1000;                  % maximum correlation lag for peak search (used as positive and as negative value)

    prm.rho = 0.7;                      % whitening parameter for rho-CSP GCC
    prm.m = 4;                          % root power for CSP-m GCC
    prm.gccWeight = 'rhoCsP';           % GCC method : CSPm | rhoCSP | PHAT

    %% Derived params

    prm.N2 = prm.N/prm.N1;
    prm.N2_pt_DFT_matrix = exp(-1i*2*pi/prm.N2).^([0:prm.N2-1].'*[0:prm.N2-1]);
    prm.twid = exp(-1i*2*pi/prm.N).^((0:prm.N1-1).'*[0:prm.N2-1]);
    prm.N1_pt_DFT_matrix = exp(-1i*2*pi/prm.N1).^([0:prm.N1-1].'*[0:prm.N1-1]);

    % Definition of the frequency mask based on Hann windows
    prm.fMask = zeros(1, prm.N);
    bMin = max(1, round(prm.f0*prm.N/prm.fs));
    bMax = round(prm.f1*prm.N/prm.fs);

    prm.bMax_N1 = ceil(bMax / prm.N2);

    windowMask = hann(bMax-bMin+1);
    prm.fMask(1:length(windowMask)) = windowMask;
    prm.fMask = circshift(prm.fMask, bMin);
    prm.fMask(prm.N/2+1:end) = fliplr(prm.fMask(1:prm.N/2));

    prm.bMin_N1 = max(floor(bMin / prm.N2), 1);
    prm.N1_pt_DFT_matrix_reduced_up = prm.N1_pt_DFT_matrix(prm.bMin_N1:prm.bMax_N1, :);
    prm.N1_pt_DFT_matrix_reduced_low = prm.N1_pt_DFT_matrix(prm.N1-prm.bMax_N1+1:prm.N1-prm.bMin_N1+1, :);

    % for IFFT only
    prm.bInv = ceil(prm.maxlag / prm.N2);
    prm.N1_pt_IDFT_matrix_reduced_up = prm.N1_pt_DFT_matrix(1:prm.bInv, :);
    prm.N1_pt_IDFT_matrix_reduced_low = prm.N1_pt_DFT_matrix(prm.N1-prm.bInv:end, :);

    %% Init States
    states.cur_frame_idx = 1;

    states.matrix_sig_out_reduced = zeros(prm.N1, prm.N2);
    states.matrix_sigref_out_reduced = zeros(prm.N1, prm.N2);
    states.matrix_ifft_out_reduced = zeros(prm.N1, prm.N2);

    states.gcc = zeros(1, prm.N); % GCC in frequency domain
    states.delay = 0; % current delay

endfunction

