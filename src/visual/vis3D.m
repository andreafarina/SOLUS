% load(['/scratch0/NOT_BACKED_UP/gdisciac/',...
%     'SOLUS_spectralfitUCL_bis/SOLUS/'...
%     'example/example_200202_series/US_prior_exp_t0.001/ex_bulk10_0.2_incl10_0.4/Test_Standard_REC.mat'])

% load('/cs/research/medim/gdisciac/SOLUS/example/example_202103_series/firstgate_ex_bulk10_0.1_incl10_0.2/Test_Standard_typeUSprior_tau2e+00_selfnorm_REC.mat')
% mask3D = REC.solver.prior.refimage;
% 
% m = mask3D;
function h = vis3D(m,REC, offset)
m = padarray(m,[1,1,1]);
m = (m - min(m(:)))/( max(m(:)) - min(m(:)) );

xlim_sq = find(squeeze( sum(sum(m,3),2))~=0);
xlim_sq(1) = xlim_sq(1)-1;
xlim_sq(end) = xlim_sq(end)+1;
ylim_sq = find(squeeze( sum(sum(m,3),1))~=0);
ylim_sq(1) = ylim_sq(1)-1;
ylim_sq(end) = ylim_sq(end)+1;
zlim_sq = find(squeeze( sum(sum(m,2),1))~=0);
zlim_sq(1) = zlim_sq(1)-1;
zlim_sq(end) = zlim_sq(end)+1;

m = m(xlim_sq(1):xlim_sq(end),ylim_sq(1):ylim_sq(end),zlim_sq(1):zlim_sq(end));

a = -1 * REC.grid.dx+ REC.grid.dx*(0:(xlim_sq(end)- xlim_sq(1)));
b = -1 * REC.grid.dy + REC.grid.dy*(0:(ylim_sq(end)- ylim_sq(1)));
c = -1 * REC.grid.dz+ REC.grid.dz*(0:(zlim_sq(end)- zlim_sq(1)));
h = figure;
if exist('offset','var')
    a = a +0.5*(offset(1)-max(a(:)));
    b = b +0.5*(offset(2)-max(b(:)));
    c = c +0.5*(offset(3)-max(c(:)));
end
isosurface(a,b,c,flip(permute(m, [2,1,3]),3), 0.5,'noshare'), axis image
xlabel('x(mm)'),ylabel('y(mm)'), zlabel('z(mm)')
colormap([0.25,0.25,0.25])
camlight left
axis image

%print('Prior_Large.pdf','-dpng')
end
