function my_savefig(hin,save_name)
% my_savefig(hin,save_name)
% enlarge the figure with handle "hin" to fullscreen and sves it under the name "save_name" in format png,eps and pdf




set(gcf, 'Units','Inches');
pos = get(hin,'Position');
set(hin,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hin,[save_name,'.pdf'], '-dpdf','-opengl')
print(hin,[save_name,'.png'], '-dpng')
print(hin,[save_name,'.svg'], '-dsvg')
end
