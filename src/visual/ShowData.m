function ShowData(dmask,data,nwin)
% SHOWDATA Plots CW or TR data in the form of the permutation matrix dmask.
% optionally you can select the number of time-windows to show
% Copyright A. Farina 

y = zeros(size(dmask));
nt = size(data,1);
if nargin > 2
    delta = round(nt/nwin);
else
    delta = 1;
    nwin = nt;
end
%Suplot grid
hsub = ceil(sqrt(nwin));
vsub = ceil(nwin/hsub);

cmin = full(min(data(:)));
cmax = full(max(data(:)));
clim = [cmin cmax];
k = 1;
for i = 1:delta:nt
    y (dmask) = data(i,:);
    subplot(hsub,vsub,k);
    if clim(2)>clim(1)
        imagesc(y,clim),
    else
        imagesc(y),
    end
    axis equal tight;colorbar;
    title(['gate ',num2str(i)]),
    xlabel('Sources');ylabel('Meas');
    k = k + 1;
end

