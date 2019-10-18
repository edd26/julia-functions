

"""
    plotimg(matrix_to_plot)

Display an image as a plot. The values from the input matrix are adjusted to the
value range of [0, 1].

If @cut_off is true then the matrix values above 256 are set to 256 and then all
values are normalized to the value 256. If @cut_off is false, then values are
normalized to maximal value.
"""
function plotimg(matrix_to_plot, cut_off=false)
    matrix_type = typeof(matrix_to_plot)
    min_val = findmin(matrix_to_plot)[1]
    int_types_arr = [Matrix{UInt8}; Matrix{UInt16}; Matrix{UInt32};
                    Matrix{UInt64}; Matrix{UInt128}; Matrix{Int8};
                    Matrix{Int16}; Matrix{Int32}; Matrix{Int64};
                    Matrix{Int128}]
    float_types_arr = [Matrix{Float16} Matrix{Float32} Matrix{Float64}]

    if min_val<0
        matrix_to_plot = shift_to_non_negative(matrix_to_plot)
    end

    max_val = findmax(matrix_to_plot)[1]

    if max_val > 256 && cut_off
        matrix_to_plot[findall(x -> x>256, matrix_to_plot)] = 256
    end

    if in(matrix_type, int_types_arr)
        matrix_to_plot = normalize_to_01(matrix_to_plot)
    elseif in(matrix_type, float_types_arr)
        matrix_to_plot = normalize_to_01(matrix_to_plot, max_val)
    end

    return colorview(Gray, matrix_to_plot)
end
