function PlotContrast(cw_t,cw_exp,gate_t,gate_exp)
if isempty(gate_t)
    gate_t = zeros(size(gate_exp,1),8,8);
end
clim(1) = min(min(cw_t(:),cw_exp(:)));
clim(2) = max(max(cw_t(:),cw_exp(:)));
h = figure;h.NumberTitle = 'off';h.Name = ['Contrast teo cw'];
imagesc(squeeze(cw_t),clim), axis image; colormap pink;
% imagesc(squeeze(cw_t)), axis image; colormap pink;
colorbar('location','southoutside');
h = figure;h.NumberTitle = 'off';h.Name = ['Contrast exp cw'];
imagesc(squeeze(cw_exp),clim), axis image; colormap pink;
% imagesc(squeeze(cw_exp)), axis image; colormap pink;
colorbar('location','southoutside');


ngate = size(gate_t,1);
P=numSubplots(ngate);
h = figure;h.NumberTitle = 'off';h.Name = ['Contrast teo Gated'];
clim(1) = min(min(gate_t(:),gate_exp(:)))-eps;
clim(2) = max(max(gate_t(:),gate_exp(:)))+eps;
for in = 1:ngate
    subplot(P(1),P(2),in)
    imagesc(squeeze(gate_t(in,:,:)),clim), axis image; colormap pink;
%     imagesc(squeeze(gate_t(in,:,:))), axis image; colormap pink;
    colorbar('location','southoutside');
    title(['Gate ' num2str(in)])
end
h = figure;h.NumberTitle = 'off';h.Name = ['Contrast exp Gated'];
for in = 1:ngate
    subplot(P(1),P(2),in)
    imagesc(squeeze(gate_exp(in,:,:)),clim), axis image; colormap pink;
%     imagesc(squeeze(gate_exp(in,:,:))), axis image; colormap pink;
    colorbar('location','southoutside');
    title(['Gate ' num2str(in)])
end

end