using LinearAlgebra

"""
    shift_to_non_negative(matrix)

Returns a matrix in which values are non-negative. This is done by finding the
minimal value in the input matrix and adding its absolute value to the matix
elements.
"""
function shift_to_non_negative(matrix)

    min_val = findmin(matrix)[1]
    if min_val < 0
        return matrix .-= min_val
    else
        return matrix
    end
end


"""
    normalize_to_01(matrix, norm_factor=256)

Returns a matrix which values are in range [0, 1]. If the values in the input
matrix are below 0, then they are shifted so that only positive numbers are in
the matrix. If the values in the matrix of shifted matrix excced value of the
@norm_factor parameter, then the matrix is normalized to the maximal value from
the matrix.
"""
function normalize_to_01(matrix; norm_factor=256)
    normalized_matrix = copy(matrix)
    min_val = findmin(normalized_matrix)[1]
    max_val = findmax(normalized_matrix)[1]

    if min_val < 0
        normalized_matrix .+= abs(min_val)
    end

    if max_val > norm_factor
        @warn "Values normalized to maximal value, not notmalization factor."
        normalized_matrix = normalized_matrix./max_val
    else
        normalized_matrix = normalized_matrix./norm_factor
    end
    return normalized_matrix
end


"""
    get_ordered_matrix(input_matrix)

Takes a `input_matrix` and returns ordered form of this matrix.
The values of returned matrix represent position in descending ordering of
of input matrix.

If `input_matrix` is symmetric, then the returned values are ordered above the
diagonal. Lower diagonal is symetrically copied from values above diagonal.

# Examples
```julia-repl
julia> a = [0 11 12; 11 0 13; 12 13 0];
julia> get_ordered_matrix(a)
3×3 Array{Int64,2}:
 0  3  2
 3  0  1
 2  1  0
```
"""
function get_ordered_matrix(input_matrix)
    data_copy = copy(input_matrix)
    mat_size = size(data_copy,1)
    ordered_matrix = zeros(Int, mat_size)

    if issymmetric(data_copy)
        symetry_order = true
    else
        symetry_order = false
        @warn "Doing non-symetric ordering"
    end

    # ====
    # Get all cartesian indices to be sorted
    matrix_indices = CartesianIndices((1:mat_size, 1:mat_size))
    if symetry_order
        matrix_indices = findall(x->x[1]>x[2], matrix_indices)
    end

    # Get all values which will be sorted
    sorting_values = data_copy[matrix_indices]

    # Sort indices by values (highest to lowest)
    ordered_indices = sort!([1:size(matrix_indices,1);],
                        by=i->(sorting_values[i],matrix_indices[i]), rev=true)

    # Put evrything together
    if symetry_order
        # how many elements are above diagonal
        repetitions = Int(ceil((size(data_copy)[1] * (size(data_copy)[1]-1))/2))
    else
        # how many elements are in whole matrix
        repetitions = Int(ceil((size(data_copy)[1] * size(data_copy)[1])))
    end

    for k=1:repetitions
        next_position = ordered_indices[k]
        matrix_index = matrix_indices[next_position]
        ordered_matrix[matrix_index] = k
        ordered_matrix[matrix_index[2], matrix_index[1]] = k
    end

    # ====
    max_orig = (findmax(input_matrix)[2])
    max_new = (findall(x->x==1,ordered_matrix)[1])
    @debug "Original maximal value was at position: " max_orig
    @debug "After ordering the first index value is at position: " max_new
    return ordered_matrix
end


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



"""
    get_local_correlations(video_array, centers, sub_img_size, shift)

Computes the correlation between the subimages and subimages shifted by values
from range -@shift:@shift and returns array with frames of size
length(@centers) x length(@centers) with the number of frames equal to the
number of rames in @video_array.

Each of the subimage is center around values stored in  @centers
"""
function get_local_correlations(video_array, centers, sub_img_size, shift)
    half_size = ceil(Int,(sub_img_size-1)/2)
    half_range = half_size + shift
    h, w, len = get_video_dimension(video_array)
    extracted_pixels = zeros(sub_img_size, sub_img_size, len)

    for frame = 1:len
        img = video_array[frame]
        for index_x = 1:size(centers,2)
            c_x = centers[2, index_x]
            for index_y = 1:size(centers,2)
                c_y = centers[1, index_y]
                subimage = img[(c_x-half_range):(c_x+half_range),
                                (c_y-half_range):(c_y+half_range)]
                center = img[(c_x-half_size):(c_x+half_size), (c_y-half_size):(c_y+half_size)]

                for left_boundary = 1:(2*shift+1)
                    for lower_boundary = 1:(2*shift+1)
                        corelation = center .* subimage[left_boundary:left_boundary+sub_img_size-1, lower_boundary:lower_boundary+sub_img_size-1]
                        corelation = sum(corelation)
                        extracted_pixels[index_x, index_y, frame] += corelation
                    end
                end
                extracted_pixels[index_x, index_y, frame] /= 256*(sub_img_size^2)*(shift*2)^2
            end
        end
    end
    return extracted_pixels
end



"""
    reduce_arrs_to_min_len(arrs)

Takes vector of vectors of different length and returns array of arrays which
are of the same length. Length in the output is the shortest vector length from
the input- values above this size are discarded.
"""
function reduce_arrs_to_min_len(arrs)
    new_arr = copy(arrs)

    simulation = size(new_arr,1)
    min_size = Inf
    for m=1:simulation
        @debug "Simulation number" m
        current_size = size(new_arr[m],1)
        @debug "Current size: " current_size
        if convert(Float64,current_size) < min_size
            min_size = current_size
            @debug "min size changed to: " min_size
        end
    end
    # min_size = Int.(min_size)
    @debug "Concatenating"
    for m=1:simulation
        new_arr[m] = new_arr[m][1:min_size,:]
    end
    min_size = Inf
    return new_arr
end
