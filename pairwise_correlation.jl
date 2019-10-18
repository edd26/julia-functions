
"""
    get_pairwise_correlation_matrix(vectorized_video, tau_max=25)

Computes pairwise correlation of the input signals accordingly to the formula
presented in paper "Clique topology reveals intrinsic geometric structure
in neural correlations" by Chad Giusti et al.

The Computations are done only for upper half of the matrix, the lower half is
a copy of upper half. Computation-wise the difference is at level of 1e-16, but
this causes that inverse is not the same as non-inverse matrix.

"""
function get_pairwise_correlation_matrix(vectorized_video, tau_max=25)
    number_of_signals = size(vectorized_video,1)
    T = size(vectorized_video,2)

    C_ij = zeros(number_of_signals,number_of_signals);
    # log_C_ij = zeros(number_of_signals,number_of_signals);

     # this is given in frames
    lags = -tau_max:1:tau_max


    for row=1:number_of_signals
        for column=row:number_of_signals
            signal_ij = vectorized_video[row,:];
            signal_ji = vectorized_video[column,:];

            # cross_corelation
            ccg_ij = crosscov(signal_ij, signal_ji, lags);
            ccg_ij = ccg_ij ./ T;

            A = sum(ccg_ij[tau_max+1:end]);
            B = sum(ccg_ij[1:tau_max+1]);
            r_i_r_j = 1;
            C_ij[row, column] = max(A, B)/(tau_max*r_i_r_j);
            C_ij[column, row] = C_ij[row, column]
            # log_C_i_j[row, column] = log10(abs(C_ij[row, column]));
        end
    end

    return C_ij
end
