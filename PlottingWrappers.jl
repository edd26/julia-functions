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
                                				),)
    heat_map = heatmap(matrix,  color=:lightrainbow,
                    title=plt_title,
                    size=plot_params.size, dpi=plot_params.dpi,
                    ticks=0:tick_step:tick_end);
    yflip_matrix && plot!( yflip = true,);

    return heat_map
end
