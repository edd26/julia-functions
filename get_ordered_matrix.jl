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

    # repetitions is equal to number of elements above diagonal= (n*(n-2))/2
    repetitions = Int(ceil((size(data_copy)[1] * (size(data_copy)[1]-1))/2))
    for k=1:repetitions
        value, index = findmax(data_copy)
        data_copy[index[1], index[2]] = min_value
        data_copy[index[2], index[1]] = min_value
        ordered_matrix[index[1], index[2]] = k
        if symetry_order
            ordered_matrix[index[2], index[1]] = k
        end
    end
    # ordered_matrix .*= (-1)
    # findmax(ordered_matrix)
    @info "Original maximal value was at position: " (findmax(input_matrix)[2])
    @info "After ordering the first index value is at position: " (findall(x->x==1,ordered_matrix)[1])
    @info "\n"

    return ordered_matrix
end
