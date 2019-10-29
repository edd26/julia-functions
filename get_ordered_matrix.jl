using LinearAlgebra

"""
    get_ordered_matrix(input_matrix)

Takes a `input_matrix` and returns ordered form of this matrix.
The values of returned matrix represent position in descending ordering of
of input matrix.

If `input_matrix` is symmetric, then the returned values are ordered above the
diagonal. Lower diagonal is symetrically copied from values aboce diagonal.

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
    input_matrix = matrix_set[m][1:5,1:5]
    data_copy = copy(input_matrix)
    if issymmetric(data_copy)
        symetry_order = true
    else
        symetry_order = false
        @warn "Doing non-symetric ordering"
    end
    # data_copy .+= abs(findmin(data_copy)[1])
    # data_copy ./= abs(findmax(data_copy)[1])
    ordered_matrix = zeros(Int, size(data_copy))
    min_value = -Inf

    # ====
    # Sorting instead of subsequent findmax
    # Get all cartesian indices to be sorted
    ordered_matrix_new = zeros(Int, size(data_copy))
    upper_indices = findall(x->x>0, UpperTriangular(data_copy))

    # Get all values above diagonal
    upper_values = data_copy[upper_indices]

    # Sort indices by values
    order = sort!([1:size(upper_indices,1);], by=i->(upper_values[i],upper_indices[i]), rev=true)

    # Put evrything together
    repetitions = Int(ceil((size(data_copy)[1] * (size(data_copy)[1]-1))/2))
    for k=1:repetitions
        next_position = order[k]
        ordered_matrix_new[upper_indices[next_position]] = k
        ordered_matrix_new[upper_indices[next_position][2], upper_indices[next_position][1]] = k
    end
    # ====
    # ordered_matrix .*= (-1)
    # findmax(ordered_matrix)
    @info "Original maximal value was at position: " (findmax(input_matrix)[2])
    @info "After ordering the first index value is at position: " (findall(x->x==1,ordered_matrix)[1])
    @info "\n"

    return ordered_matrix
end
