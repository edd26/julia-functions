# new_component
# using Eirene
using DelimitedFiles
 	using Plots
    using LinearAlgebra
    using Images
    using Distances
    using Images
    using JLD


# cd("../matrix-tests/")
    julia_func_path = "./julia-functions/"
    # include(julia_func_path*"MatrixToolbox.jl")
    # include(julia_func_path*"MatrixProcessing.jl")
    # include(julia_func_path*"BettiCurves.jl")
    # include(julia_func_path*"ImageProcessing.jl")
    include(julia_func_path*"PlottingWrappers.jl")

	include("PointsSubstitution.jl")

"""
	function expand_matrix(input_matrix, expansion_size, last_components;do_plot=false)

Takes 'input_matrix' (an ordering matrix used for creating cliques) and and adds
2×'expansion_size' number of rows. 'last_components' are the values in original
matrix that are added last to the clique.

Results may be be plotted by setting 'do_plot=true'.
"""
function expand_matrix(input_matrix, expansion_size, last_components;do_plot=false)
	new_comp = last_components
	matrix_size = size(input_matrix,1)
	for mat_sizes = matrix_size:2:(matrix_size+2expansion_size)
		input_matrix, new_comp = add_step_to_matrix(input_matrix, new_comp)
	end
	if do_plot
		expand_plt_ref = plot_square_heatmap(input_matrix, 1,size(input_matrix,1);
											plt_title = "Original, size:$(matrix_size)",
											color_palete=:lightrainbow)
		display(expand_plt_ref)
	end
	return input_matrix, expand_plt_ref
end

# Shuffle matrix entries
"""
	function shuffle_matrix(input_matrix, shuffles; do_plot=false)

Takes symmetric 'input_matrix' and randomly swaps rows 'shuffles' many times.

Results may be plotted by setting 'do_plot=true'.
"""
function shuffle_matrix(input_matrix, shuffles; do_plot=false)
	matrix_size = size(input_matrix,1)
	rows = randcycle(matrix_size)
	shuffled_ord_mat = copy(input_matrix)

	for k = 1:shuffles
		# global shuffled_ord_mat, rows
		srcs, trgts = rand(rows,2)
		swap_rows!(shuffled_ord_mat, srcs, trgts)
	end

	if do_plot
		shuff_plt_ref = plot_square_heatmap(shuffled_ord_mat, 1,size(shuffled_ord_mat,1);
											plt_title = "Shuffled, size:$(matrix_size)",
											color_palete=:lightrainbow)
		display(shuff_plt_ref)
	end
	return input_matrix, shuff_plt_ref
end



"""
	function organize_shuff_matrix(input_matrix; do_plots=false)

Reorganizes 'input_matrix' so that values highest values in a row are positioned
next to the diagonal.

Results may be plotted by setting 'do_plot=true'.
"""
function organize_shuff_matrix(input_matrix; do_plots=false)
	unscrambled_matrix = copy(input_matrix)
	matrix_size = size(input_matrix,1)
	for k = matrix_size:-2:2
		max_row_val = findmax(unscrambled_matrix[k,:])[2]
		# put to the previous last position
		swap_rows!(unscrambled_matrix, max_row_val, k-1)
		# skip 1 row and work on next one
	end
	if do_plots
		reorganized_plt_ref = plot_square_heatmap(unscrambled_matrix, 1,size(unscrambled_matrix,1);
									plt_title = "unscrambled_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		display(reorganized_plt_ref)
	end
	return unscrambled_matrix, reorganized_plt_ref
end


"""
		function order_max_vals_near_diagonal(input_matrix; do_plots=false, direction=:descending)

Orders values in 'input_matrix' so that values next to diagonal are descending
(by default).

TODO- not working-  Optionally, ascending order can be used by setting 'direction' to
':ascending'.

Results may be plotted by setting 'do_plot=true'.
"""
function order_max_vals_near_diagonal(input_matrix; do_plots=false, direction=:descending)
	# Find max values next to the diagonal
	matrix_size = size(input_matrix,1)

	if direction == :descending
		# ordering_function = findmax
		new_ord_value = -1
		iteration_values = matrix_size:-2:2
		# iteration_values = 2:2:matrix_size

	elseif direction == :ascending
		# ordering_function = findmin
		new_ord_value = findmax(input_matrix)[1]*2
		iteration_values = 2:2:matrix_size
	else
		# @error "Unknow ordering was given"
		throw("Unknow ordering was given")
	end

	reordered_matrix = copy(input_matrix)
	row_indices = 1:2:matrix_size
	col_indices = 2:2:matrix_size
	coord_set = [CartesianIndex(row_indices[k], col_indices[k]) for k=1:matrix_size÷2]
	diag_max_values = reordered_matrix[coord_set]

	for k = iteration_values
		max_val, max_ind = findmax(diag_max_values)
		(direction == :descending) ? (position = floor(k÷2)) : (position = floor(k÷2))
		diag_max_values[max_ind] = diag_max_values[position]
		diag_max_values[position] = new_ord_value
		max_ind *= 2

		swap_rows!(reordered_matrix, k, max_ind)
		swap_rows!(reordered_matrix, k-1, max_ind-1)
	end


	if do_plots
		reorganized_plt_ref = plot_square_heatmap(reordered_matrix, 1,size(reordered_matrix,1);
									plt_title = "reordered_matrix, size:$(matrix_size)",
									color_palete=:lightrainbow)
		display(reorganized_plt_ref)
	end
	return reordered_matrix, reorganized_plt_ref
end


"""
	function fine_tune_matrix(input_matrix; do_plots=false)

Check if velues next to the maximal values are organized in descending order.
"""
function fine_tune_matrix(input_matrix; do_plots=false)#, direction=:descending)
	# Find max values next to the diagonal
	matrix_size = size(input_matrix,1)
	fine_tune_matrix = copy(input_matrix)

	# if direction == :descending
	# 	# ordering_function = findmax
	# 	new_ord_value = -1
	# 	iteration_values = matrix_size:-2:2
	# 	# iteration_values = 2:2:matrix_size
	#
	# elseif direction == :ascending
	# 	# ordering_function = findmin
	# 	new_ord_value = findmax(input_matrix)[1]*2
	# 	iteration_values = 2:2:matrix_size
	# else
	# 	# @error "Unknow ordering was given"
	# 	throw("Unknow ordering was given")
	# end


	for k = 2:2:matrix_size-1
		if fine_tune_matrix[k-1,k+1] > fine_tune_matrix[k,k+1]
			swap_rows!(fine_tune_matrix, k, k-1)
		end
	end


	if do_plots
		fine_tuned_plt_ref = plot_square_heatmap(fine_tune_matrix, 1,size(reordered_matrix,1);
										plt_title = "fine_tuned, size:$(matrix_size)",
										color_palete=:lightrainbow)
		display(fine_tuned_plt_ref)
	end
	return fine_tune_matrix, fine_tuned_plt_ref
end
