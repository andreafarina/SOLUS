function [Mua,Musp,Conc,A,B] = applyGrid(grid,muaB,muspB,solv_type,concB,aB,bB)
nL = numel(muaB);
if iscolumn(muaB), muaB = muaB'; end
if iscolumn(muspB), muspB = muspB'; end
Mua = ones(prod(grid.dim),nL).* muaB;
Musp = ones(prod(grid.dim),nL).* muspB;
Mua = squeeze(reshape(Mua,[grid.dim nL]));
Musp = squeeze(reshape(Musp,[grid.dim nL]));

if contains(solv_type,'spectral')
    nC = numel(concB);
    if iscolumn(concB), concB = concB'; end
    Conc = ones(prod(grid.dim),nC).* concB;
    A = ones(prod(grid.dim),1).* aB;
    B = ones(prod(grid.dim),1).* bB;
    Conc = squeeze(reshape(Conc,[grid.dim nC]));
    A = squeeze(reshape(A,[grid.dim 1]));
    B = squeeze(reshape(B,[grid.dim 1]));
else
   Conc = [];A =[];B = []; 
end
end
