function roi = SelectROI(data,irf)
% plot time-domain data y, allow user to select a ROI. 
% If the irf is provided it will be plotted as weel (in green)

figure;
semilogy(data),hold on
if nargin > 1
    semilogy(irf/max(irf(:))*max(data(:)),'--g')
end
ylim([max(data(:))/10000 max(data(:))])%,legend(['data','irf']),
title('select ROI');
%r = drawrectangle;
% x(1) = r.Position(1);
% x(2) = r.Position(3) + x(1);
[x,~] = ginput(2);
roi = round(x);
disp(['Used ROI: ',num2str(roi(1)),'--' ,num2str(roi(2))]);
close