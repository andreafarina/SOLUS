function roi = SelectROI(data,irf)
% plot time-domain data y, allow user to select a ROI. 
% If the irf is provided it will be plotted as weel (in green)

figure;
semilogy(data),hold on
if nargin > 1
    semilogy(irf/max(irf)*max(data(:)),'g')
end
ylim([max(data(:))/10000 max(data(:))])%,legend(['data','irf']),
title('select ROI');
[x,~] = ginput(2);
roi = round(x);
disp(['Used ROI: ',num2str(roi')]);
close