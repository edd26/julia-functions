using Plots
"""
    plot_square_heatmap(matrix, tick_step, tick_end;
                                plt_title, img_size=(900, 800), img_dpi=300)

Takes matrix and plots it as a heatmap. Funtion returns the handler to the
heatmap.
"""
function plot_square_heatmap(matrix, tick_step, tick_end;
                                plt_title="", yflip_matrix=true,
                                 plot_params= (dpi=300,
                                				size=(900,800),
                                				lw=1,
                                				thickness_scaling=1,
                                				top_margin= 0,
                                				left_margin=[0 0],
                                				bottom_margin= 0
                                				),
                                color_palete=:lightrainbow,
                                add_labels=true)
    heat_map = heatmap(matrix,  color=color_palete,
                    title=plt_title,
                    size=plot_params.size, dpi=plot_params.dpi,
                    ticks=0:tick_step:tick_end);
    yflip_matrix && plot!( yflip = true,);

    if add_labels
        xlabel!("Matrix index")
        ylabel!("Matrix index")
    end

    return heat_map
end



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



"""
   plot_image_analysis(plots_set; description::NamedTuple, original_img, kwargs...)

Takes set of plots and puts them in 2 coulm layout. If 'description' is given,
adds entry with the data processing description. If 'original_img' is given, it
is also displayed next to the descrtions field.

'kwargs' are plot properties.
"""
function plot_image_analysis(plots_set; description::NamedTuple, original_img, kwargs...)
    kwarg_keys = kwargs.keys()

    (!isempty(original)) ? (orig_img_flag = true) : (orig_img_flag = false)
    (!isempty(description)) ? (desc_flag = true) : (desc_flag = false)

    l = @layout [a{0.2w} [grid(3,3) b{0.2h}]]

    total_plot_sets = 7
    total_cols = 2
    total_rows = ceil(Int,total_plot_sets/total_cols)
    if orig_img_flag || desc_flag
        total_rows +=1
    end

    height_unit = 1/total_rows

    matrix =   [1 2 3;
                4 5 6;
                7 8 9]

    l = @layout [a{0.4w,} b{0.4w,};
                    # grid(1,4);
                    # grid(1,4);
                    # grid(1,4);
                    # grid(1,4);
                    ]
                    # [grid(2,2)  grid(2,2)]]
                    # [a [grid(4,2) b]]]
    data = [rand(10, 4), rand(11, 4)]

    l = @layout [a{0.4w} b
                c d e f
                c d e f
                c d e f
                c d e f
                c d e f]
     ref = plot(grid=false,
            axis=false,
            layout = l,
            legend = false,
            # seriestype = [:scatter :path],
            dpi=300,
            size=(900,1200),
            )

    ref.series_list

    p2 = plot!(ref.subplots[18],rand(10, 1),seriestype = :scatter,axis=true,grid=true, title="")
    p2 = plot!(ref.subplots[21],rand(10, 10),seriestype = :heatmap, legend=true, xlabel="index", ylabel="index")
    annotate!(ref.subplots[0], 0, 0, "my text", :red)
    p1.subplots
    # color scheme

end
