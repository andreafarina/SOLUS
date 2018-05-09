clear all
close all

t = linspace(0,3500,3500)';
mua = 0;
mus = 1;
XX = [20 20 20 20]; YY = [0 0 0 0]; ZZ = [5 10 20 30]; ampli = [1.2 15 300 3000];
ri = [0 0 0];
rj = [40 0 0];
nB = 1.4;
v = 0.2999/nB;
nE = 1.4;
A = A_factor(nB/nE);
start = tic;

for ii =1:1
    figure
    tic
    J = J_Semi_Inf_PCBC_TR_3D_muaD(t,mua,mus,v,A,ri,rj,XX,YY,ZZ,1,'A');
    toc
    %plot(J.*100), grid, ylim([-2e-11 2e-11]);
    [phi,~]=ForwardTD([],ri,rj,1,mua,mus,nB,[],[],A,1,3500,0,'semi-inf','homo');
    %plot(phi)
    plot(J./repmat(phi,[1 4]),'LineWidth',2)
    
    ylim([-4e-2 1e-2]);
    xlabel Time(ps)
    ylabel('Relative Absorption Perturbation')
    set(gca,'FontSize',22)
    legend(strcat(repmat({'Z2 = '},4,1), num2str([5 10 20 30]')),'Location','best')
    ah = gca;
    ah.XAxis.MinorTick = 'on';
    ah.XAxis.MinorTickValues = 250:250:3250;
    ah.YAxis.MinorTick = 'on';
    ah.YAxis.MinorTickValues = -0.035:0.01:0.01;
    grid
    ah.YMinorGrid = 'on'; ah.XMinorGrid = 'on';
    ah.GridAlpha = 0.4; ah.MinorGridAlpha = 0.7;
    
    figure
    tic
    J = J_Semi_Inf_PCBC_TR_3D_muaD(t,mua,mus,v,A,ri,rj,XX,YY,ZZ,1,'D');
    toc
    [phi,area]=ForwardTD([],ri,rj,1,mua,mus,nB,[],[],A,1,3500,0,'semi-inf','homo');
    plot(J.*ampli./repmat(phi,[1 4]),'LineWidth',2)
    ylim([-4e-2 6e-2]);
    xlabel Time(ps)
    ylabel('Relative Scattering Perturbation')
    set(gca,'FontSize',22)
    legend(strcat(repmat({'Z2 = '},4,1), num2str([5 10 20 30]')),'Location','best')
    ah = gca;
    ah.XAxis.MinorTick = 'on';
    ah.XAxis.MinorTickValues = 250:250:3250;
    ah.YAxis.MinorTick = 'on';
    ah.YAxis.MinorTickValues = -0.03:0.02:0.06;
    grid
    ah.YMinorGrid = 'on'; ah.XMinorGrid = 'on';
    ah.GridAlpha = 0.4; ah.MinorGridAlpha = 0.7;
end
toc(start);
