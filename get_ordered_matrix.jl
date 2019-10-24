using LinearAlgebra

""" Takes a symmetrc matrix and returns ordered form of this matrix.
The values of returned matrix represent position in descending ordering of
of input matrix."""
function get_ordered_matrix(symmetric_matrix)
    data_copy = copy(symmetric_matrix)
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
    repetitions = ceil((size(data_copy)[1] * (size(data_copy)[1]-1))/2)
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
    @info "Original maximal value was at position: " (findmax(symmetric_matrix)[2])
    @info "After ordering the first index value is at position: " (findall(x->x==1,ordered_matrix)[1])
    @info "\n"

    return ordered_matrix
end
