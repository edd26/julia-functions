"""
    plot_square_heatmap(matrix, tick_step, tick_end;
                                title, img_size=(900, 800), img_dpi=300)

Takes matrix and plots it as a heatmap. Funtion returns the handler to the
heatmap.
"""
function plot_square_heatmap(matrix, tick_step, tick_end;
                                title, img_size=(900, 800), img_dpi=300)
    heat_map = heatmap(matrix,  color=:lightrainbow,
                    title="Order matrix of $(file_n), tiles size:$(sub_img_size)",
                    size=img_size, dpi=img_dpi,
                    ticks=0:tick_step:tick_end)

    return heat_map
end
