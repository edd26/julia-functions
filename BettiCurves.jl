using Plots
using Eirene
using Measures
using Plots.PlotMeasures
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
function plot_bettis(bettis, plot_title; legend_on=true, min_dim=0, plot_params= (dpi=300,
				size=(900,800),
				lw=1,
				thickness_scaling=1,
				top_margin= `:match`,
				left_margin= `:match`,
				bottom_margin= `:match`
				))
                                        #                        base_dpi = 500)
    # set_default_plotting_params()
    cur_colors = get_color_palette(:auto, plot_color(:white), 17)
    colors_set =  [cur_colors[7], cur_colors[5], [:red], cur_colors[1], cur_colors]

    # final_title = "Eirene betti curves, "*plot_title
	plot_ref = plot(title=plot_title);
    max_dim = size(bettis)[1]
    for p = (1+min_dim):(max_dim)
        plot!(bettis[p][:,1], bettis[p][:,2], label="\\beta_"*string(p-1),
                                                    lc=colors_set[p],
													dpi=plot_params.dpi,
													size=plot_params.size,
													lw=plot_params.lw,
													thickness_scaling=plot_params.thickness_scaling,
													top_margin=plot_params.top_margin,
								                    left_margin=plot_params.left_margin,
								                    bottom_margin=plot_params.bottom_margin
													);
        if legend_on
            plot!(legend=true)
        else
            plot!(legend=false)
        end

    end
    ylabel!("Number of cycles")
	xlabel!("Edge density")
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
								results_path::String; file_name="",
								extension = ".png", data_size::String="",
								do_save=true,
								extend_title=true, do_normalise=true, min_dim=0,
								max_dim=3, legend_on=true, plot_params)

    bettis = get_bettis(eirene_results, max_dim);
    norm_bettis = normalise_bettis(bettis);
    plot_ref = plot_bettis(bettis, plot_title, legend_on=legend_on, min_dim=min_dim, plot_params=plot_params);

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
bettis_eirene(matr, maxdim; mintime=-Inf, maxtime=Inf, numofsteps=Inf, mindim=1)

Takes the `matr` and computes Betti curves up to `maxdim`. Return matrix only
with betti curve values


Function taken from: https://github.com/alexyarosh/hyperbolic
"""
function bettis_eirene(matr, maxdim;
							mintime=-Inf, maxtime=Inf, numofsteps=Inf, mindim=1)
    c = eirene(matr, minrad = mintime, maxrad= maxtime, numrad= numofsteps, maxdim=maxdim)

    int_length = maxtime-mintime
    step_length= int_length/numofsteps

    if (mintime == -Inf) || (maxtime == Inf) || (numofsteps == Inf)
		@debug "Inf mintime, maxtime or number of steps."
        # return [betticurve(c, dim=maxdim) for d=1:maxdim]
		result = vectorize_bettis(c, maxdim, mindim)
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
            ordered_matrix = get_ordered_matrix(matrix_set[m][subset_range,subset_range];
													assing_same_values=false)
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


"""
	function get_curves_from_matrix(img_name)

Computes Betti curves for the image file indicated by @img_name. If the image is
	not symmetric, then it is the elements below diagonal are copied over the
	elmenents above the diagonal.
"""
function get_curves_from_matrix(img_name, plot_params; file_path="",
									plot_heatmaps = true, save_heatmaps=false,
								plot_betti_figrues = true)
  file_n = split(img_name, ".")[1]
  img1_gray = Gray.(load(file_path*img_name))
  img_size = size(img1_gray)

  C_ij = Float64.(img1_gray)

  if !issymmetric(C_ij)
    img1_gray = symmetrize_image(img1_gray)
    C_ij = Float64.(img1_gray)
  end
  img_size = size(C_ij,1)
  # C_ij =-C_ij
  # C_ij .+= 1


  # ==============================================================================
  # =============================== Ordered matrix ===============================
  if size(C_ij,1) > 80
    @warn "Running Eirene for big matrix: " img_size
    @warn "Eirene may have trobules with big matrices/images."
  end

  ordered_matrix = get_ordered_matrix(C_ij; assing_same_values=false)


  # ==============================================================================
  # ============================ Persistance homology ============================
  C = eirene(ordered_matrix,maxdim=3,model="vr")


  # ==============================================================================
  # ================================ Plot results ================================

  if plot_heatmaps

    full_ordered_matrix= get_ordered_matrix(C_ij; assing_same_values=false)
    heat_map2 = plot_square_heatmap(full_ordered_matrix, 10, img_size;
            plt_title = "Order matrix of $(file_n)", plot_params=plot_params)

    if save_heatmaps
        heatm_details = "_heatmap_$(file_n)"
        savefig(heat_map2, heatmaps_path*"ordering"*heatm_details)
    end
  end

  if plot_betti_figrues
    plot_title = "Betti curves of $(file_n), size=$(img_size) "
    figure_name = "betti_$(file_n)_n$(img_size)"
    ref = plot_and_save_bettis(C, plot_title, figure_path, ; 		file_name=figure_name, plot_params=plot_params,
                                    do_save=false, extend_title=false,
    								do_normalise=false, max_dim=3,legend_on=true,
                                    min_dim=1)
  end
  display(img1_gray)
  display(heat_map2)
  display(ref)
end

"""
multiscale_matrix_testing(sample_space_dims = 3,
                                    maxsim=5,
                                    min_B_dim = 1,
                                    max_B_dim = 3,
                                    size_start = 10,
                                    size_step = 5,
                                    size_stop = 80; do_random=true)

Function for testing the average number of cycles from geometric and random
    matrices.

It is possible to save intermidiate results- for that, @control_saving must be
set true.

Performance of computation of Betti curves can be monitored, if the
@perform_eavl is set too true. Bydefault, it is set to false.
"""
function multiscale_matrix_testing(sample_space_dims = 3,
                                    maxsim=5,
                                    min_B_dim = 1,
                                    max_B_dim = 3,
                                    size_start = 10,
                                    size_step = 5,
                                    size_stop = 50;
										do_random=true, control_saving=false,
										perform_eavl=false)
    num_of_bettis = length(collect(min_B_dim:max_B_dim))

    if length(sample_space_dims) > 1
        @warn "Can not do random processing for multiple dimensions"
        do_random =false
    end

    geom_mat_results = Any[]
    if do_random
        rand_mat_results = Any[]
        result_list = [geom_mat_results, rand_mat_results]
    else
        result_list = [geom_mat_results]
    end

    for sample_space_dim in sample_space_dims
		if !do_random; @info "Sampling space size: " sample_space_dim; end

        repetitions = size_start:size_step:size_stop
        for space_samples in repetitions
            @info "Generating data for: " space_samples
            # ==========================================
            # ============= Generate data ==============
            # ===
            # Generate random matrix
            if do_random
                symm_mat_rand = [generate_random_matrix(space_samples) for
		 															i=1:maxsim]
                ordered_mat_rand = [get_ordered_matrix(symm_mat_rand[i];
									assing_same_values=false) for i=1:maxsim]
            end

            # ===
            # Generate geometric matrix
            pts_rand = [generate_random_point_cloud(sample_space_dim,space_samples) for i=1:maxsim]
            symm_mat_geom = [generate_geometric_matrix(pts_rand[i]') for i=1:maxsim]
            ordered_mat_geom = [get_ordered_matrix(symm_mat_geom[i];
									assing_same_values=false) for i=1:maxsim]

            # ======================================================================
            # ========================= Do the Betti analysis ======================
            if do_random
                set = [ordered_mat_geom, ordered_mat_rand]
            else
                set = [ordered_mat_geom]
            end
            for matrix_set in set
                @debug("Betti analysis!")
                # ===
                # Generate bettis
				many_bettis = Array[]
				if perform_eavl
					many_timings = Float64[]
					many_bytes = Float64[]
					many_gctime = Float64[]
					many_memallocs = Base.GC_Diff[]
				end

                for i=1:maxsim
					if i%10 == 0
                    	@info "Computing Bettis for: " i
					end

					if perform_eavl
						results, timing, bytes, gctime, memallocs =
					 			@timed bettis_eirene(matrix_set[i],
												max_B_dim, mindim=min_B_dim)
	                    push!(many_bettis,results)
						push!(many_timings,timing)
						push!(many_bytes,bytes)
						push!(many_gctime,gctime)
						push!(many_memallocs,memallocs)
					else
						push!(many_bettis, bettis_eirene(matrix_set[i],
												max_B_dim, mindim=min_B_dim))
					end
                end

                # ===
                # Get maximal number of cycles from each Betti from simulations
                max_cycles = zeros(maxsim, max_B_dim)
                for i=1:maxsim,  betti_dim = 1:max_B_dim
                    @debug("\tFindmax in bettis")
                    max_cycles[i, betti_dim] = findmax(many_bettis[i][:, betti_dim])[1]
                end

                # ===
                # Get the statistics
                avg_cycles = zeros(1, length(min_B_dim:max_B_dim))
                std_cycles = zeros(1, length(min_B_dim:max_B_dim))
                k=1
                for betti_dim=min_B_dim:max_B_dim
                    avg_cycles[k] = mean(max_cycles[:, betti_dim])
                    std_cycles[k] = std(max_cycles[:, betti_dim])
                    k+=1
                end

                # ===
                # Put results into dictionary
                betti_statistics = Dict()
                if matrix_set == ordered_mat_geom
                    @debug("Saving ordered")
                    betti_statistics["matrix_type"] = "ordered"
                    betti_statistics["space_dim"] = sample_space_dim
                    result_list = geom_mat_results
                else
                    @debug("Saving radom")
                    betti_statistics["matrix_type"] = "random"
                    result_list = rand_mat_results
                end
                betti_statistics["space_samples"] = space_samples
                betti_statistics["simualtions"] = maxsim
                betti_statistics["min_betti_dim"] = min_B_dim
                betti_statistics["max_betti_dim"] = max_B_dim
                betti_statistics["avg_cycles"] = avg_cycles
                betti_statistics["std_cycles"] = std_cycles

				if perform_eavl
					betti_statistics["many_timings"] = many_timings
					betti_statistics["many_bytes"] = many_bytes
					betti_statistics["many_gctime"] = many_gctime
					betti_statistics["many_memallocs"] = many_memallocs
				end
                push!(result_list, betti_statistics)
            end # matrix type loop
            @debug("===============")
			if control_saving
				if do_random
					save("multiscale_matrix_testing_$(space_samples)_$(sample_space_dim).jld",
				 							"rand_mat_results", rand_mat_results,
	                                        "geom_mat_results", geom_mat_results)
	            else
					save("multiscale_matrix_testing_dimension_$(space_samples)_$(sample_space_dim).jld",
				 						  "geom_mat_results", geom_mat_results)
	            end
			end
        end # matrix_size_loop
    end # sampled space dimension

    if do_random
        return geom_mat_results, rand_mat_results
    else
        return geom_mat_results
    end
end


"""
	vectorize_bettis(eirene_results, maxdim, mindim)

Returns the betti curves in the form of matrix, which rows are Betti values, and
columns are Betti dimensions starting with @mindim up to @maxdim.
"""
function vectorize_bettis(eirene_results, maxdim, mindim)
	number_of_steps = length(betticurve(eirene_results, dim=0)[:,1])
	number_of_bettis = maxdim-mindim+1

	result = zeros(number_of_steps, number_of_bettis)
	try
		iter=1
		for d=mindim:maxdim
			bett_res = betticurve(eirene_results, dim=d)[:,2]
			if size(bett_res,1) == 0
				@warn "Computed betti curve had 0 elements, creating vector with zeros"
				bett_res = zeros(size(result,1))
			end
			result[:,iter] = bett_res
			iter+=1
		end
		@debug size(result)
		return result
	catch err
		if isa(err, DimensionMismatch)
			@error "Dimension mismatch error"
			for d=mindim:maxdim
				@error size(betticurve(c, dim=d))
			end
			@error hcat([betticurve(c, dim=d)[:,2] for d=mindim:maxdim]...)
			throw(err)
	   else
			@error "Unknown error occurred"
			throw(err)
	   end
	end
end


 # ===============================================
function get_curves_from_matrix2(img_name; file_path="",
									plot_heatmaps = true, save_heatmaps=false,
								plot_betti_figrues = true)
  file_n = split(img_name, ".")[1]
  img1_gray = Gray.(load(file_path*img_name))
  img_size = size(img1_gray)

  C_ij = Float64.(img1_gray)

  if !issymmetric(C_ij)
    img1_gray = symmetrize_image(img1_gray)
    C_ij = Float64.(img1_gray)
  end
  img_size = size(C_ij,1)
  # C_ij =-C_ij
  # C_ij .+= 1


  # ==============================================================================
  # =============================== Ordered matrix ===============================
  if size(C_ij,1) > 80
    @warn "Running Eirene for big matrix: " img_size
    @warn "Eirene may have trobules with big matrices/images."
  end

  ordered_matrix = get_ordered_matrix(C_ij; assing_same_values=false)


  # ==============================================================================
  # ============================ Persistance homology ============================
  C = eirene(ordered_matrix,maxdim=3,model="vr")


  # ==============================================================================
  # ================================ Plot results ================================

  if plot_heatmaps

    full_ordered_matrix= get_ordered_matrix(C_ij; assing_same_values=false)
    heat_map2 = plot_square_heatmap(full_ordered_matrix, 10, img_size;
            plt_title = "Order matrix of $(file_n)")

    if save_heatmaps
        heatm_details = "_heatmap_$(file_n)"
        savefig(heat_map2, heatmaps_path*"ordering"*heatm_details)
    end
  end

  if plot_betti_figrues
    plot_title = "Betti curves of $(file_n), size=$(img_size) "
    figure_name = "betti_$(file_n)_n$(img_size)"
    ref = plot_and_save_bettis2(C, plot_title, figure_path, ; 		file_name=figure_name, plot_params=plot_params,
                                    do_save=false, extend_title=false,
    								do_normalise=false, max_dim=3,legend_on=true,
                                    min_dim=1)
  end
  display(img1_gray)
  display(heat_map2)
  display(ref)
end

function plot_and_save_bettis2(eirene_results, plot_title::String,
								results_path::String; file_name="",
								extension = ".png", data_size::String="",
								do_save=true,
								extend_title=true, do_normalise=true, min_dim=0,
								max_dim=3, legend_on=true)

    bettis = get_bettis(eirene_results, max_dim);
    norm_bettis = normalise_bettis(bettis);
    plot_ref = plot_bettis2(bettis, plot_title, legend_on=legend_on, min_dim=min_dim);

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


function plot_bettis2(bettis, plot_title; legend_on=true, min_dim=0)#; plot_size = (width=1200, height=800),
                                        #                        base_dpi = 500)
    # set_default_plotting_params()
    cur_colors = get_color_palette(:auto, plot_color(:white), 17)
    colors_set =  [cur_colors[7], cur_colors[5], [:red], cur_colors[1], cur_colors]

    # final_title = "Eirene betti curves, "*plot_title
	plot_ref = plot(title=plot_title);
    max_dim = size(bettis)[1]
    for p = (1+min_dim):(max_dim)
        plot!(bettis[p][:,1], bettis[p][:,2], label="\\beta_"*string(p-1)
													);
        if legend_on
            plot!(legend=true)
        else
            plot!(legend=false)
        end

    end
    ylabel!("Number of cycles")
	xlabel!("Edge density")
    return plot_ref
end
