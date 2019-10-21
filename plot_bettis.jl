using Plots
using Eirene
include("save_figures.jl")
using Measures



# Source: https://github.com/JuliaPlots/Plots.jl/issues/897
function set_default_plotting_params(;upscale=2)
    #8x upscaling in resolution
    fntsm = Plots.font("sans-serif", pointsize=round(12.0*upscale))
    fntlg = Plots.font("sans-serif", pointsize=round(18.0*upscale))
    default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
    default(size=(800*upscale,600*upscale)) #Plot canvas size
    default(dpi=500) #Only for PyPlot - presently broken
end


"""
Uses betticurve function to generate range of Betti curves.
"""
function get_bettis(results_eirene, max_dim)
    bettis  = Matrix{Float64}[]
    for d =1:(max_dim+1)
        result = betticurve(results_eirene, dim=d-1)
        push!(bettis, result)
    end
    return bettis
end

"""
Normalise the horizontal values of the Eirene betti numbers.
"""
function normalise_bettis(bettis)
    norm_bettis = copy(bettis)
    @debug "norm_bettis size :" size(norm_bettis)[1][1]

    max_dim = size(norm_bettis)[1]
    @debug "typeof(max_dim) :" typeof(max_dim[1])

    for d =1:(max_dim)
        if !isempty(norm_bettis[d])
            norm_bettis[d][:,1] /= findmax(norm_bettis[d][:,1])[1]
        end
    end
    return norm_bettis
end

"""
Creates a plot for set of betti numbers.
"""
function plot_bettis(bettis, plot_title)#; plot_size = (width=1200, height=800),
                                        #                        base_dpi = 500)
    # set_default_plotting_params()
    cur_colors = get_color_palette(:auto, plot_color(:white), 17)
    colors_set =  [cur_colors[7], cur_colors[5], [:red], cur_colors[1], cur_colors]

    final_title = "Eirene betti curves, "*plot_title*" data, size "

    plot_ref = plot(title=final_title);
    max_dim = size(bettis)[1]
    for p = 1:(max_dim)
        plot!(bettis[p][:,1], bettis[p][:,2], label="\\beta_"*string(p-1),
                                            # size=plot_size, dpi = base_dpi,
                                                margin=30mm,   lc=colors_set[p]);
    end
    ylabel!("Number of cycles")
    return plot_ref
end


"""
Plot Betti curves from 0 up to max_dim using results from Eirene library and
returns handler for figure. Optionally, saves the figure or normalise the
    horizontal axis to maximal value
"""
function plot_and_save_bettis(eirene_results, plot_title::String,
                                results_path::String;data_size::String="",
                               do_save=true, do_normalise=true, max_dim=3)
    bettis = get_bettis(eirene_results, max_dim);
    norm_bettis = normalise_bettis(bettis);
    plot_ref = plot_bettis(bettis, plot_title);

    if do_save
        savefig(plot_ref, "betti_curves_"*plot_title*data_size*".png")
    end
    return plot_ref
end


"""
Plots Set of betti numbers at the same graph as well as each of them separatelay.
"""
function plot_decomposed_bettis(results_eirene, dataset_name)
    max_dim = 3;
    bettis = get_bettis(results_eirene, max_dim)

    # p0 = plot(betti_0[:,1], betti_0[:,2], label="\\beta_0", linecolor=:blue); #, ylims = (0,maxy)
    p1 = plot(bettis[1][:,1], bettis[1][:,2], label="\\beta_"*string(1),
                                        linecolor=:orange, legend=:topleft);
    p2 = plot(bettis[2][:,1], bettis[2][:,2], label="\\beta_"*string(1),
                                        linecolor=:red, legend=:topleft);
    p3 = plot(bettis[3][:,1], bettis[3][:,2], label="\\beta_"*string(1),
                                        linecolor=:steelblue, legend=:topleft);

    p0 = plot(bettis[1][:,1], bettis[1][:,2], label="\\beta_1", linecolor=:orange,
        legend=:topleft, title="All Bettis for "*dataset_name*" data, Eirene");
    plot!(bettis[2][:,1], bettis[2][:,2], label="\\beta_2", linecolor=:red);
    plot!(bettis[3][:,1], bettis[3][:,2], label="\\beta_3", linecolor=:steelblue);
    return plot_ref = plot(p0, p1, p2, p3, layout=4)
end
