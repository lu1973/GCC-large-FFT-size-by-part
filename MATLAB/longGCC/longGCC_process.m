## -*- texinfo -*-
## @deftypefn  {} {[@var{delay_out}, @var{states}] =} longGCC_process(@var{sig}, @var{refsig}, @var{prm}, @var{states})
## compute the GCC between @var{sig} and @var{refsig} following the parameters defined in @var{prm}.
## returns the delay estimated in @var{delay_out}
##
## based on https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/66959/versions/2/previews/fft_algorithm.m/index.html
##
## @end deftypefn

function [delay_out, states] = longGCC_process(sig, refsig, prm, states)
    delay_out = 0;
    sig = sig(1:prm.N1, 1);
    refsig = refsig(1:prm.N1, 1);

    % compute partial FFT for the signal and the reference
    states = longGCC_partial_sig_fft(sig, prm, states);
    states = longGCC_partial_sigref_fft(refsig, prm, states);

    % compute partial IFFT of the GCC
    states = longGCC_partial_ifft(states.gcc((states.cur_frame_idx-1)*prm.N1+1:states.cur_frame_idx*prm.N1), prm, states);

    % compute full FFT and GCC in the frequency domain
    if states.cur_frame_idx == prm.N2
        % convert matrix into vector with GCC
        ifft_flat = states.matrix_ifft_out_reduced(:, :).';

        gcc_time = transpose(ifft_flat(:)) / prm.N;
        gcc_time_flipped = [gcc_time(1), fliplr(gcc_time(2:end))];
        gcc_time_abs = abs(gcc_time_flipped);
        gcc_time_abs = [gcc_time_abs(end-prm.maxlag:end), gcc_time_abs(1:prm.maxlag)];

        % Find index of the max
        [gccMax,idxMax] = max(gcc_time_abs);
        if gccMax > 0
            % Lagrange polynomial peak extraction
            g1 = gcc_time_abs(idxMax-1);
            g2 = gccMax;
            g3 = gcc_time_abs(idxMax+1);
            f = (g1-g3)/(2*g1-4*g2+2*g3);
            % bias correction
            if abs(f) > eps
                f = (sqrt(32*f^2+1)-1)/(8*f);
            endif
            % compensate offset
            states.delay = idxMax+f - prm.maxlag;
        endif
        states.matrix_ifft_out_reduced = prm.lambda * states.matrix_ifft_out_reduced;

        mic_flat = states.matrix_sig_out_reduced.';
        xFFT = transpose(mic_flat(:));

        xFFT_std = xFFT .* prm.fMask;

        ref_flat = states.matrix_sigref_out_reduced(:, :).';
        rFFT = transpose(ref_flat(:));

        % introduce here freq masking
        rFFT = rFFT .* prm.fMask;
        est_corr = xFFT_std .* conj(rFFT);
        abs_corr_1 = max(abs(xFFT_std), eps);
        abs_corr_2 = max(abs(rFFT), eps);
        abs_corr = max(abs(est_corr), eps);
        if strcmp(prm.gccWeight, 'rhoCSP')
            states.gcc(1, :) = est_corr./(abs_corr.^prm.rho);
        elseif strcmp(prm.gccWeight, 'CSPm')
            states.gcc(1, :) = est_corr./(abs_corr_1.*abs_corr_2).^(1/prm.m);
        else
            states.gcc(1, :) = est_corr./abs_corr;
        endif

        delay_out = states.delay;
        states.cur_frame_idx = 0;

        states.matrix_sig_out_reduced = prm.lambda * states.matrix_sig_out_reduced;
        states.matrix_sigref_out_reduced = prm.lambda * states.matrix_sigref_out_reduced;
    endif
    states.cur_frame_idx = states.cur_frame_idx + 1;
endfunction


