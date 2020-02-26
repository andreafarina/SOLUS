function roi = SelectROI(data,irf)
% plot time-domain data y, allow user to select a ROI. 
% If the irf is provided it will be plotted as weel (in green)

<<<<<<< HEAD
figure();
=======
figure;
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
semilogy(data),hold on
if nargin > 1
    semilogy(irf/max(irf)*max(data(:)),'g')
end
ylim([max(data(:))/10000 max(data(:))])%,legend(['data','irf']),
<<<<<<< HEAD
title('select ROI');
=======
title('select ROI')
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
[x,~] = ginput(2);
roi = round(x);
disp(['Used ROI: ',num2str(roi')]);
close