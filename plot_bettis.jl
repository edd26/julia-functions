using Plots

function plot_bettis(topology_res, dataset_name)
    betti_0 = betticurve(topology_res, dim=0)
    betti_1 = betticurve(topology_res, dim=1)
    betti_2 = betticurve(topology_res, dim=2)
    betti_3 = betticurve(topology_res, dim=3)

    # p0 = plot(betti_0[:,1], betti_0[:,2], label="\\beta_0", linecolor=:blue); #, ylims = (0,maxy)
    p1 = plot(betti_1[:,1], betti_1[:,2], label="\\beta_1",
                                        linecolor=:orange, legend=:topleft);
    p2 = plot(betti_2[:,1], betti_2[:,2], label="\\beta_2",
                                        linecolor=:red, legend=:topleft);
    p3 = plot(betti_3[:,1], betti_3[:,2], label="\\beta_3",
                                        linecolor=:steelblue, legend=:topleft);

    p0 = plot(betti_1[:,1], betti_1[:,2], label="\\beta_1", linecolor=:orange,
        legend=:topleft, title="All Bettis for "*dataset_name*" data, Eirene");
    plot!(betti_2[:,1], betti_2[:,2], label="\\beta_2", linecolor=:red);
    plot!(betti_3[:,1], betti_3[:,2], label="\\beta_3", linecolor=:steelblue);
    return plot_ref = plot(p0, p1, p2, p3, layout=4)
end

"""
Plot Betti curves from 0 up to max_dim using results from Eirene library and
returns handler for figure. Optionally, saves the figure or normalise the
    horizontal axis to maximal value
"""
function plot_and_save_bettis(eirene_results, plot_title,
                               results_path;
                               do_save=true, do_normalise=true, max_dim=3)

    bettis  = Matrix{Float64}[]
    for d =1:(max_dim+1)
        result = betticurve(eirene_results, dim=d-1)
        push!(bettis, result)
        if do_normalise && !isempty(bettis[d])
            bettis[d][:,1] /= findmax(bettis[d][:,1])[1]
        end
    end

    cur_colors = get_color_palette(:auto, plot_color(:white), 17)
    colors_set =  [cur_colors[7], cur_colors[5], [:red], cur_colors[1], cur_colors]

    final_title = "Eirene betti curves, "*plot_title*" data, size "

   plot_ref = plot(title=final_title);
   for p = 1:(max_dim+1)
       plot!(bettis[p][:,1], bettis[p][:,2], label="\\beta_"*string(p-1), lc=colors_set[p]);
   end
    ylabel!("Number of cycles")

    if do_save
        current_path = pwd()
        cd(results_path)
        savefig(plot_ref, "betti_curves_"*plot_title*data_size*".png")
        cd(current_path)
    end
    return plot_ref
end
