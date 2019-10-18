"""
Plot betti curves obtained from Eirene and the correlation matrix.

Note: Exported from TestingPairwiseCorrelationmatrix.jl
"""
function plot_eirene_betti_curves(C, C_ij)
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

    heat_map1 = heatmap(C_ij,  color=:lightrainbow, title="Cij, $(choice), number of points: $points_per_dim");
end
