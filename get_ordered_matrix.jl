"""Takes a symmetrc matrix and returns ordered form of this matrix with the
entries representing order of values within all matrix"""
function get_ordered_matrix(symmetric_matrix)
    data_copy = copy(symmetric_matrix)
    data_copy .+= findmin(data_copy)[1]
    # data_copy ./= findmax(data_copy)[1]
    ordered_matrix = zeros(Int, size(data_copy)) .+2

    repetitions = ceil((size(data_copy)[1] * (size(data_copy)[1]-1))/2)
    for k=1:repetitions
        value, index = findmax(data_copy)
        data_copy[index[1], index[2]] = 0
        data_copy[index[2], index[1]] = 0
        ordered_matrix[index[1], index[2]] = k
        ordered_matrix[index[2], index[1]] = k
    end
    # ordered_matrix .*= (-1)
    # findmax(ordered_matrix)
    @info "Original maximal value was at position: " (findmax(symmetric_matrix)[2])

    @info "After ordering the maximal value is at position: " (findmax(ordered_matrix)[2])

    return ordered_matrix
end
