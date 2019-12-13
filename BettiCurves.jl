using Plots
using Eirene
using Measures
include("MatrixProcessing.jl")


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
function get_bettis(results_eirene::Dict, max_dim::Integer)
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

Creates a plot for set of betti numbers stored in `bettis` and return the
handler to the plot.
`plot_title` is used for the title of the plot.
"""
function plot_bettis(bettis, plot_title; legend_on=true, min_dim=0)#; plot_size = (width=1200, height=800),
                                        #                        base_dpi = 500)
    # set_default_plotting_params()
    cur_colors = get_color_palette(:auto, plot_color(:white), 17)
    colors_set =  [cur_colors[7], cur_colors[5], [:red], cur_colors[1], cur_colors]

    # final_title = "Eirene betti curves, "*plot_title

    plot_ref = plot(title=plot_title);
    max_dim = size(bettis)[1]
    for p = (1+min_dim):(max_dim)
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
z
Plot Betti curves from 0 up to `max_dim` using `eirene_results` from Eirene library and
returns handler for figure. Optionally, if `do_save` is set, saves the figure
or if `do_normalise` is set, sets the steps range to be normalised to the
horizontal axis maximal value.
"""
function plot_and_save_bettis(eirene_results, plot_title::String,
								results_path::String; file_name="",
								extension = ".png", data_size::String="",
								do_save=true,
								extend_title=true, do_normalise=true, min_dim=0,
								max_dim=3, legend_on=true)
    bettis = get_bettis(eirene_results, max_dim);
    norm_bettis = normalise_bettis(bettis);
    plot_ref = plot_bettis(bettis, plot_title, legend_on=legend_on, min_dim=min_dim);

    if do_save
		if extend_title && isempty(file_name)
			file_name = "betti_c_"*plot_title*data_size*extension;
		elseif isempty(file_name)
			file_name = plot_title*extension
		elseif isempty(findall(x->x==extension[2:end], split(file_name, ".")))
			#check for the extension in file name
			file_name *= extension
		end

        savefig(plot_ref, results_path*file_name)
        @info "Saved file as " results_path*file_name

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


# TODO merge functions for getting betti curves
# Original function returns 2 different types of betti curves. If no default
# value parameters is given, it returns vector of matrices. If num of steps is
# given, then it return matrix maxdim x numsteps.
"""
bettis_eirene(matr, maxdim; mintime=-Inf, maxtime=Inf, numofsteps=Inf)

Takes the `matr` and computes Betti curves up to `maxdim`. Return matrix only
with betti curve values


Function taken from: https://github.com/alexyarosh/hyperbolic
"""
function bettis_eirene(matr, maxdim; mintime=-Inf, maxtime=Inf, numofsteps=Inf)
    c = eirene(matr, minrad = mintime, maxrad= maxtime, numrad= numofsteps, maxdim=maxdim)

    int_length = maxtime-mintime
    step_length= int_length/numofsteps

    if (mintime == -Inf) || (maxtime == Inf) || (numofsteps == Inf)
        # return [betticurve(c, dim=maxdim) for d=1:maxdim]
        return hcat([betticurve(c, dim=d)[:,2] for d=1:maxdim]...)
    end

    betts = zeros(numofsteps, maxdim)
    # For every dimension compute betti curve
    for dim=1:maxdim
        bet = betticurve(c, dim=dim)

        #for every element in betti curve return betti value if index is positive
        for i=1:size(bet,1)
            b = bet[i,:]
            ind = Int(ceil((b[1]-mintime)/step_length))
            if ind > 0
                betts[ind,dim]=b[2]
            else
                betts[1,dim]=b[2]
            end
        end
    end
    return betts
end


"""
	get_avg_bettis_from_JLD(data_sets; range=-1, maxsim=-1, steps=-1,
													subset_size=-1, maxdim=3)

Takes the 'data_sets' (vector of dictionaries from loading data with JLD) and
computes average betti curves with their std's.

"""
function get_avg_bettis_from_JLD(data_sets; range=-1,
                                maxsim=-1, steps=-1, subset_size=-1, maxdim=3)
# TODO change name- it does not use JLD file
    avg_betti = Array[]
    std_betti = Array[]

    if maxsim == -1
        maxsim=size(data_sets[1]["dist"], 1)
    end

    if range == -1 || range > size(data_sets,1)
        range = 1:size(data_sets,1)
    else
        range = 1:range
    end

    for k = range
        @debug "== Next k" k
        matrix_set = data_sets[k]["dist"]
        bettis_set = Array[]
        # steps = 2600;

        new_bettis= []
        for m=1:maxsim

            if subset_size == -1
                subset_range=1:size(matrix_set[1],1)
            else
                subset_range=1:subset_size
            end

            @debug "=== Next m" m
            ordered_matrix = get_ordered_matrix(matrix_set[m][subset_range,subset_range])
            new_bettis = bettis_eirene(ordered_matrix, maxdim)
            push!(bettis_set, new_bettis)
        end

        # Set size to be the same for every betti curve
        new_bettis_set = reduce_arrs_to_min_len(bettis_set)

        # Compute average betti curve
        push!(avg_betti, average_bettis(new_bettis_set))
        push!(std_betti, std_bettis(new_bettis_set))
    end

    return avg_betti, std_betti
end
