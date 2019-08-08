function make_plotly_figure(path, fig_file)
fig = strcat(path, fig_file);
new_fig = plotlyfig(gcf);
fig2plotly(fig)
end
