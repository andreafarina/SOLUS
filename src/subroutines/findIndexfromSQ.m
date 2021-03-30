function dm = findIndexfromSQ(dmask,s,d)

dmask = double(dmask);
dmask(s,d) = 2;
dummy = dmask(dmask >= 1);
dm = find(dummy == 2);
return

end
