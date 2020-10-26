% load('C:\Users\monia\Desktop\PROVE\15-2m\0.01spectral_tk1-ROItuttalacurva_fattorediconversione\Test_Standard_REC')
function [bmua, bmusp, bConc] = remove_voxels(bmua, bmusp, bConc, dim, mua0, musp0, conc0)
for j = 1 : (dim(1)*dim(2))
    bConc(j,:) = conc0(:);
end
for j = 1 : (dim(1)*dim(2))
    bmua(j,:) = mua0(:);
    bmusp(j,:) = musp0(:);
end
end

% openfig('fig.fig')
% roi = drawcircle('Color','k','FaceAlpha',0.4);
% mask = createMask(roi);
% imshow(mask)

