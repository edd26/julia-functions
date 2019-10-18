"""
Save figures with given parameters.

Note: Exported from TestingPairwiseCorrelationmatrix.jl
"""
function save_figures(ref, path, choice, points_per_dim, tau, size_limiter)
    name = split(choice, ".")[1]
    name =  path * name *
            "_size$(size_limiter)_points$(points_per_dim)_tau$(tau).png"
    savefig(ref, name)

    @info "File saved: " name
end
