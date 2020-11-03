function [] = SubPlotMap(X,titlefig,nfigure,Nr,Nc,isub,res_fact)
%SUBPLOTMAP Plot a montage of the volume data X

figure(nfigure),subplot(Nr,Nc,isub),imagesc(GenerateMontageMatrix(...
X,res_fact)),title(titlefig),...
colorbar, axis equal, axis tight, axis off;
end

