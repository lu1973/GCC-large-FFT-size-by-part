% Compute partial FFT for ref signal
function states = longGCC_partial_sigref_fft(y, prm, states)
    % stepl: performing N2-point DFTs of each one the N1 rows
    % optimization remark: in Matlab, it is better to do this large
    % multiplication than using the fact the N2 matrix is
    % symmetrical and doing a memory operation (flip+conj)
    matrixl = y*(prm.N2_pt_DFT_matrix(:, states.cur_frame_idx).');
    % step2: multiplying elements of matrix with twiddle factor W_N^{n1k2}
    matrix2 = matrixl.*prm.twid;
    % step3: performing N1-point DFTs of each column of array
    prod_up = prm.N1_pt_DFT_matrix_reduced_up*matrix2;
    states.matrix_sigref_out_reduced(prm.bMin_N1:prm.bMax_N1, :) = states.matrix_sigref_out_reduced(prm.bMin_N1:prm.bMax_N1, :) + (1-prm.lambda) * prod_up;
    prod_down = prm.N1_pt_DFT_matrix_reduced_low*matrix2;
    states.matrix_sigref_out_reduced(prm.N1-prm.bMax_N1+1:prm.N1-prm.bMin_N1+1, :) = states.matrix_sigref_out_reduced(prm.N1-prm.bMax_N1+1:prm.N1-prm.bMin_N1+1, :) + (1-prm.lambda) * prod_down;
endfunction
