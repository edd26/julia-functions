"""
    normalize_to_01(matrix, norm_factor=256)

Returns a matrix which values are in range [0, 1]. If the values in the input
matrix are below 0, then they are shifted so that only positive numbers are in
the matrix. If the values in the matrix of shifted matrix excced value of the
@norm_factor parameter, then the matrix is normalized to the maximal value from
the matrix.
"""
function normalize_to_01(matrix, norm_factor=256)
    normalized_matrix = copy(matrix)
    min_val = findmin(normalized_matrix)[1]
    max_val = findmax(normalized_matrix)[1]

    if min_val < 0
        normalized_matrix .-= min_val
    end

    if max_val > norm_factor
        @warn "Values normalized to maximal value, not notmalization factor."
        normalized_matrix = normalized_matrix./max_val
    else
        normalized_matrix = normalized_matrix./norm_factor
    end
    return normalized_matrix
end
