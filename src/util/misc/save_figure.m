function save_figure(FName)
set(gcf, 'Color', 'w');
set(gcf, 'PaperPosition', [-0.5 -0.25 6 5.5]); %Position the plot further to the left and down. Extend the plot to fill entire paper.
set(gcf, 'PaperSize', [5 5]); %Keep the same paper size
hold off;
saveas(gcf,FName,'fig')
%export_fig test -eps -pdf -jpg
%export_fig(FName,'-painters','-eps','-pdf','-jpg');
export_fig(FName,'-painters','-jpg');
end
