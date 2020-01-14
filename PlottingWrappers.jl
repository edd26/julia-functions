"""
    plot_square_heatmap(matrix, tick_step, tick_end;
                                plt_title, img_size=(900, 800), img_dpi=300)

Takes matrix and plots it as a heatmap. Funtion returns the handler to the
heatmap.
"""
function plot_square_heatmap(matrix, tick_step, tick_end;
                                plt_title, img_size=(900, 800), img_dpi=300)
    heat_map = heatmap(matrix,  color=:lightrainbow,
                    title=plt_title,
                    size=img_size,
                    dpi=img_dpi,
                    # lw=2,
                    thickness_scaling=3,
                    top_margin= -33px,
                    left_margin=[-50px 0px],
                    right_margin=[15px 0px],
                    bottom_margin= -35px,
                    ticks=0:tick_step:tick_end);

    return heat_map
end
