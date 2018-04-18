function PlotHeteQM(DOT,data,bkg)
%% Plots the configuration of optodes and represents the heterogeneities


%% plot sources
plot3(DOT.Source.Pos(:,1),DOT.Source.Pos(:,2),DOT.Source.Pos(:,3),'r*','MarkerSize',20),grid,
xlabel('x'),ylabel('y'),zlabel('z'),hold on
%% plot detectors
plot3(DOT.Detector.Pos(:,1),DOT.Detector.Pos(:,2),DOT.Detector.Pos(:,3),'bx','MarkerSize',20),
xlabel('x(mm)'),ylabel('y(mm)'),zlabel('z(mm)')

%% plot heterogeneities
RES_FACT = 1;
data = data - bkg;
data = imresizen(data,RES_FACT);
xx = imresizen(DOT.grid.x',RES_FACT);
yy = imresizen(DOT.grid.y',RES_FACT);
zz = imresizen(DOT.grid.z',RES_FACT);
[x,y,z] = ind2sub(size(data),find(data>0.5*max(data(:))));
plot3(xx(x), yy(y), zz(z), '.k','MarkerSize',15);
set(gca,'zdir','reverse'),axis equal,
xlim([DOT.grid.x1 DOT.grid.x2]),...
    ylim([DOT.grid.y1 DOT.grid.y2]),...
    zlim([DOT.grid.z1 DOT.grid.z2])
%legend('source','detector','pert');




% for i = 1:num_hete
%     h_str = ['hete',num2str(i)];
% plot3(DOT.opt.(h_str).c(1),DOT.opt.(h_str).c(2),DOT.opt.(h_str).c(3),'bo'),
% [a,b,c]= sphere(100);
% surf(a*DOT.opt.(h_str).sigma + DOT.opt.(h_str).c(1), ...
%     b*DOT.opt.(h_str).sigma + DOT.opt.(h_str).c(2),...
%     c*DOT.opt.(h_str).sigma + DOT.opt.(h_str).c(3))
% end
% set(gca,'zdir','reverse'),axis equal,
% xlim([DOT.grid.x1 DOT.grid.x2]),...
%     ylim([DOT.grid.y1 DOT.grid.y2]),...
%     zlim([DOT.grid.z1 DOT.grid.z2])



%figure(2),
%plot3(DOT.Source.Pos(:,1),DOT.Source.Pos(:,2),DOT.Source.Pos(:,3),'r*'),grid,
%xlabel('x'),ylabel('y'),zlabel('z'),hold on
% cylinder 1
%[b,c,a]= cylinder(ones(1000,1),1000);
%surf(a*DOT.opt.hete.sigma + DOT.opt.hete.c(1), ...
%    b*DOT.opt.hete.l + DOT.opt.hete.c(2),...
%    c*DOT.opt.hete.sigma + DOT.opt.hete.c(3))
%set(gca,'zdir','reverse'),axis equal,
%xlim([DOT.grid.x1 DOT.grid.x2]),...
%    ylim([DOT.grid.y1 DOT.grid.y2]),...
%    zlim([DOT.grid.z1 DOT.grid.z2])
% cylinder 2
%[c,a,b]= cylinder(ones(1000,1),1000);
%surf(a*DOT.opt.hete2.l + DOT.opt.hete2.c(1), ...
%    b*DOT.opt.hete2.sigma + DOT.opt.hete2.c(2),...
%    c*DOT.opt.hete2.sigma + DOT.opt.hete2.c(3))
%set(gca,'zdir','reverse'),axis equal,
%xlim([DOT.grid.x1 DOT.grid.x2]),...
%    ylim([DOT.grid.y1 DOT.grid.y2]),...
%    zlim([DOT.grid.z1 DOT.grid.z2])