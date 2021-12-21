function [h] = SubPlotMap(X,titlefig,nfigure,Nr,Nc,isub,res_fact)
%SUBPLOTMAP Plot a montage of the volume data X
if nargin<7
    res_fact = [1,1,1];
end
    
%h = figure(nfigure);
subplot(Nr,Nc,isub),imagesc(GenerateMontageMatrix(...
X,res_fact)),title(titlefig),...
colorbar, axis equal, axis tight, axis off;
end

