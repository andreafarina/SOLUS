function ImageMeasure(DOT,im)%% plot sources
figure(485)
scatter(DOT.Source.Pos(:,1),DOT.Source.Pos(:,2),'k','filled'),
xlabel('x'),ylabel('y'),zlabel('z'),hold on
xlim([DOT.grid.x1 DOT.grid.x2]),...
    ylim([DOT.grid.y1 DOT.grid.y2]),...
    zlim([DOT.grid.z1 DOT.grid.z2])
[rows,cols]=find(DOT.dmask==1);
RES_FACT = 5;
data = abs(DOT.opt.Mua(:,:,:,1) - DOT.opt.muaB(1));
data = imresizen(data,RES_FACT);
xx = imresizen(DOT.grid.x',RES_FACT);
yy = imresizen(DOT.grid.y',RES_FACT);
zz = imresizen(DOT.grid.z',RES_FACT);
[x,y,z] = ind2sub(size(data),find(data>0.5*max(data(:))));
plot3(xx(x), yy(y), zz(z), '.k','MarkerSize',15);
scatter(DOT.Source.Pos(cols(im),1),DOT.Source.Pos(cols(im),2),'r','filled')
scatter(DOT.Detector.Pos(rows(im),1),DOT.Detector.Pos(rows(im),2),'r','filled')
title(['Measure #', num2str(im)])

