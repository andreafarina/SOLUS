function idxQM = findMeasIndex(dmask)
% return index of measurements given by dmask as vector of cells.
% Each cell contains the index of measurements in the second dimension of
% the measurements matrix by lambda.

    
nQM_lam = sum(sum(dmask,1),2);
nQM_lam = cat(1,0, cumsum(nQM_lam(:)));

for inl=1:size(dmask,3)
    idxQM{inl} = (nQM_lam(inl)+1):nQM_lam(inl+1);
end

return
end