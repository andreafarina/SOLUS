function dm = findIndexfromSQ(dmask,d,s)
% dm = findIndexfromSQ(dmask,d,s)
% returns dm, the index of a measure in an array of fixed wavelength, TxSQ
% dmask: mask of detector-source measurements for a given lambda
% d: detector number
% s: source number
dmask = double(dmask);
dmask(d,s) = 2;
dummy = dmask(dmask >= 1);
dm = find(dummy == 2);
return

end
