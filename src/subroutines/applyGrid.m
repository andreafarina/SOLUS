function [Mua,Musp] = applyGrid(grid,muaB,muspB)
nL = numel(muaB);
if iscolumn(muaB), muaB = muaB'; end
if iscolumn(muspB), muspB = muspB'; end
Mua = ones(prod(grid.dim),nL).* muaB;
Musp = ones(prod(grid.dim),nL).* muspB;
Mua = squeeze(reshape(Mua,[grid.dim nL]));
Musp = squeeze(reshape(Musp,[grid.dim nL]));
end