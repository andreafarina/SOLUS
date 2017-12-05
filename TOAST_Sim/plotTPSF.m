%% plot TPSF for homegeonus case using diffusion equation
mua = 0.01;  % mm-1
musP = 1.0; 
n_in = 1.4;
n_ext = 1.0; % air
c0 = 0.299;  % mm/ps
Tmax = 2000; % ps
pos_Q = [0,0,0];
pos_M = [0,30,0];
%% 
t = 1:2000;
A = A_factor(n_in./n_ext);
y = SemiInfinite_TR(t,pos_Q,pos_M,mua,musP,c0./n_in,A);
figure,semilogy(t,y),ylim([max(y)./1e4 max(y)*1.1]),
xlabel('time (ps)'),ylabel('mm^{-2}ps^{-1}')

