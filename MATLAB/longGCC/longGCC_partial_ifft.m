% Compute partial IFFT
function states = longGCC_partial_ifft(y, prm, states)
    % stepl: performing N2-point DFTs of each one the N1 rows
    matrixl = (y.') * (prm.N2_pt_DFT_matrix(:, states.cur_frame_idx).');
    % step2: multiplying elements of matrix with twiddle factor W_N^{n1k2}
    matrix2 = matrixl.*prm.twid;
    % step3: performing N1-point DFTs of each column of array
    prod_up = prm.N1_pt_IDFT_matrix_reduced_up*matrix2;
    states.matrix_ifft_out_reduced(1:prm.bInv, :) = states.matrix_ifft_out_reduced(1:prm.bInv, :) + (1-prm.lambda) * prod_up;
    prod_down = prm.N1_pt_IDFT_matrix_reduced_low*matrix2;
    states.matrix_ifft_out_reduced(prm.N1-prm.bInv:end, :) = states.matrix_ifft_out_reduced(prm.N1-prm.bInv:end, :) + (1-prm.lambda) * prod_down;
endfunction
