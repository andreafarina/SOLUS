
% CW FLUENCE INSIDE AN INFINITE MEDIUM 
function phi = SemiInfinite_CW(rs,rd,mua,mus,A)
%%
% rs source position
% rd detector position

if rs(3) > 0,
    z0=rs(3);     
 elseif rs(3)==0,
    z0=1/mus;
    rs(3)=z0;
end
D = 1/(3*mus);
ze=2*A*D;
rs_min=rs; rs_min(3)=-2*ze-z0;

phi=(Infinite_CW(rs,rd,mua,mus)-Infinite_CW(rs_min,rd,mua,mus))./(2*A);
phi(isnan(phi))=0;
return