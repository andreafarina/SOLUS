function lind = linIndex(d, ind)
  %Compute linearized index
    dprod = d(1);
    lind = ind(:,1);
    for jd = 2:length(d)
      lind = lind + (ind(:,jd)-1)*dprod;
      dprod = dprod*d(jd);
    end
  end

