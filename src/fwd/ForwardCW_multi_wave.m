function phi=ForwardCW_multi_wave(grid,SourcePos,DetectorPos, dmask,...
    muaB,muspB,Mua,Musp,A, geom, TYPE_FWD, radiometry)
phi = zeros(sum(dmask(:)),radiometry.nL);

for inl = 1:radiometry.nL
    phi(:,inl) = ForwardCW(grid,SourcePos, DetectorPos, dmask, ...
        muaB(inl), muspB(inl), squeeze(Mua(:,:,:,inl)), ...
        squeeze(Musp(:,:,:,inl)), A, geom, TYPE_FWD);
end
end