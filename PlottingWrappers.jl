"""
    plot_square_heatmap(matrix, tick_step, tick_end;
                                plt_title, img_size=(900, 800), img_dpi=300)

Takes matrix and plots it as a heatmap. Funtion returns the handler to the
heatmap.
"""
function plot_square_heatmap(matrix, tick_step, tick_end;
                                plt_title, plot_params= (dpi=300,
                                				size=(900,800),
                                				lw=1,
                                				thickness_scaling=1,
                                				top_margin= 0px,
                                				left_margin=[0px 0px],
                                				bottom_margin= 0px
                                				),)
    heat_map = heatmap(matrix,  color=:lightrainbow,
                    title=plt_title,
                    size=plot_params.size, dpi=plot_params.dpi,
                    ticks=0:tick_step:tick_end)

    return heat_map
end
