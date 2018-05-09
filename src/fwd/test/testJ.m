mua = 0;
mus = 1;
n = 1.4;
nr = 1;
A = A_factor(nr);
cs = 0.2999/n;
vi = 1;
t = (1:3500)';
ri = [-20 0 0];
rj = [20 0 0];
XX = 0;
YY = 0;
ZZ = 20;
x = -30:2:30;
y = -30:2:30;
z = 0:2:50;
[XX,YY,ZZ] = ndgrid(x,y,z);
 tic
 J= J_Semi_Inf_PCBC_TR_3D_Edo(t,mua,mus,cs,A,ri,rj,XX,YY,ZZ,vi,'ad');
 toc
 %plot(J')
 
