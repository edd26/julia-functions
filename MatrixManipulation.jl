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
    ordered_matrix = zeros(Int, size(data_copy))

    if issymmetric(data_copy)
        symetry_order = true
    else
        symetry_order = false
        @warn "Doing non-symetric ordering"
    end

    # ====
    # Get all cartesian indices to be sorted
        if symetry_order
        matrix_indices = findall(x->x>0, UpperTriangular(data_copy))
    else
        matrix_indices = findall(x->x>0, data_copy)
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
