"""
Save figures with given parameters.

Note: Exported from TestingPairwiseCorrelationmatrix.jl
"""
function save_figures_old(ref, path, choice, points_per_dim, tau, size_limiter)
    name = split(choice, ".")[1]
    name =  path * name *
            "_size$(size_limiter)_points$(points_per_dim)_tau$(tau).png"
    savefig(ref, name)

    @info "File saved: " name
end


"""
Save figures.
"""
function save_figure(plot_ref, results_path, plot_title; extention=".png" )
    full_path = results_path*plot_title*extention
    savefig(plot_ref, full_path)
    @info "File saved under: " full_path
end


"""
Save betti curves.
"""
function save_betti(plot_ref, results_path, plot_title)
    full_title =  "betti_curves_"*plot_title;
    save_figure(plot_ref, results_path, full_title)
end
