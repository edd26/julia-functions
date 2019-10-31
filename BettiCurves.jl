using Plots
using Eirene
using Measures



# Source: https://github.com/JuliaPlots/Plots.jl/issues/897
function setdefaultplottingparams(;upscale=2)
    #8x upscaling in resolution
    fntsm = Plots.font("sans-serif", pointsize=round(12.0*upscale))
    fntlg = Plots.font("sans-serif", pointsize=round(18.0*upscale))
    default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
    default(size=(800*upscale,600*upscale)) #Plot canvas size
    default(dpi=500) #Only for PyPlot - presently broken
end


"""
	get_bettis(results_eirene, max_dim)

Uses betticurve function to generate Betti curves up to `max_dim` diemsion from
the `results_eirene` dictionary.
"""
function get_bettis(results_eirene<:Dict, max_dim<:Int)
    bettis  = Matrix{Float64}[]
    for d =1:(max_dim+1)
        result = betticurve(results_eirene, dim=d-1)
        push!(bettis, result)
    end
    return bettis
end


"""
	normalise_bettis(bettis)

Normalise the number of steps for every Eirene betti number in 'bettis' variable.
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
	plot_bettis(bettis, plot_title; legend_on=true)

Creates a plot for set of betti numbers stored in `bettis`.
`plot_title` is used for the title of the plot.
"""
function plot_bettis(bettis, plot_title; legend_on=true)#; plot_size = (width=1200, height=800),
                                        #                        base_dpi = 500)
    # set_default_plotting_params()
    cur_colors = get_color_palette(:auto, plot_color(:white), 17)
    colors_set =  [cur_colors[7], cur_colors[5], [:red], cur_colors[1], cur_colors]

    final_title = "Eirene betti curves, "*plot_title

    plot_ref = plot(title=final_title);
    max_dim = size(bettis)[1]
    for p = 1:(max_dim)
        plot!(bettis[p][:,1], bettis[p][:,2], label="\\beta_"*string(p-1),
                                                    lc=colors_set[p]);
        if legend_on
            plot!(legend=true)
        else
            plot!(legend=false)
        end

    end
    ylabel!("Number of cycles")
    return plot_ref
end


"""
	plot_and_save_bettis(eirene_results, plot_title::String,
								results_path::String; extension = ".png",
								data_size::String="", do_save=true,
								extend_title=true, do_normalise=true, max_dim=3,
								legend_on=true)

Plot Betti curves from 0 up to `max_dim` using `eirene_results` from Eirene library and
returns handler for figure. Optionally, if `do_save` is set, saves the figure
or if `do_normalise` is set, sets the steps range to be normalised to the
horizontal axis maximal value.
"""
function plot_and_save_bettis(eirene_results, plot_title::String,
								results_path::String; extension = ".png",
								data_size::String="", do_save=true,
								extend_title=true, do_normalise=true, max_dim=3,
								legend_on=true)
    bettis = get_bettis(eirene_results, max_dim);
    norm_bettis = normalise_bettis(bettis);
    plot_ref = plot_bettis(bettis, plot_title, legend_on=legend_on);

    if do_save
		if extend_title
			file_name = "betti_c_"*plot_title*data_size*extension;
		else
			file name = plot_title*extension
		end

        savefig(plot_ref, file_name)
        @info "Saved file as " file_name

    end
    return plot_ref
end


"""
	plot_decomposed_bettis(results_eirene, dataset_name)

Plots betti curves 0 up to 3 at the same plot. Beside, betti curves 1 to 3 are
plotted separately, in 3 additional plots at the same graph.
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


"""
	plot_ebetti_curves_heatmap(C, C_ij)

Plot betti curves obtained from Eirene and the correlation matrix.

Note: Exported from TestingPairwiseCorrelationmatrix.jl
"""
function plot_ebetti_curves_heatmap(C, C_ij)
    betti_0 = betticurve(C, dim=0)
    betti_1 = betticurve(C, dim=1)
    betti_2 = betticurve(C, dim=2)
    betti_3 = betticurve(C, dim=3)

    title = "Betti curves for pairwise corr. matrix"

    p1 = plot(betti_0[:,1], betti_0[:,1], label="beta_0", title=title)
    #, ylims = (0,maxy)
    plot!(betti_1[:,1], betti_1[:,2], label="beta_1")
    plot!(betti_2[:,1], betti_2[:,2], label="beta_2")
    plot!(betti_3[:,1], betti_3[:,2], label="beta_3")

    heat_map1 = heatmap(C_ij,  color=:lightrainbow,
					title="Cij, $(choice), number of points: $points_per_dim");
end
