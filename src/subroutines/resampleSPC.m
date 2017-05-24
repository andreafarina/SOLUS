function [data,area,peak,baric] = resampleSPC(y,t,dt,flag)
% Resample experimental data to the time-scale defined in DOT.time
% The final scale will be t(1):dt:t(end)
% INPUT
% y:    temporal curves n_time x n_meas
% t:    time-axis of data (ps)
% dt:   new sampling dt
% flag: 'norm' to get area-normalized curves
% OUTPTUT
% data: data resampled
% area: area before normalization
% peak: structure containing position,value and time of the peak
% baric: structure containing position and time of the baricenter

verbosity = 0;
%nmeas = size(y,2);
%dty = ty(2)-ty(1);
dty = t(2)-t(1);
bin_ratio = dt./dty;       % assume that y and irf time steps are in a integer ratio


%h = y.data;
% convert time scale of the IRF to fwd
%tirff = tirf;% - dtirf/2;% + IRF.baric.time;
%% Better to interp if bin is not integer
%tirf = tirff(1):0.1:tirff(end);
%h = interp1(tirff,h,tirf);
ti = t(1):dt:t(end);

% a = 1:bin:(numel(tirf)-1-bin);
% b = a + (bin-1);
% twin = [a',b'];
% irfii = WindowTPSF(h,twin);

% it create problems for boundary values of the bins
% [~,bin] = histc(t,ti);
% if ~all(bin)
%     y(bin==0) = [];
%     t(bin==0) = [];
%     bin(bin==0) = [];
% end

bin = meshgrid([1:numel(ti)],[1:uint8(bin_ratio)]);
bin = bin(1:numel(y));
bin = bin(:);
yi = accumarray(bin,y);
yi(yi<0) = 0;
nyi=numel(yi);
nti=numel(ti);
if nyi < nti
    ti(nyi+1:end) = [];
elseif nyi > nti
    yi(nti+1:end) = [];
end
%yi = yi./sum(yi);
if verbosity == 1
    figure;semilogy(ti,yi./sum(yi),t,y./sum(y));
end
% % resample y to the irf scale (expected more fine)
% 
% yi = interp1(ty,y,tirf + IRF.baric.time);
% yi(isnan(yi)) = 0;
%% Calculate area
area = sum(yi);

%% normalize data
if nargin==4
    data = yi./sum(yi);
else
    data = yi;
end
%% Calculate peak
[peak.value,peak.pos] = max(data);
peak.time = ti(peak.pos);
%% Calculate baricenter
x = 1:numel(ti);
baric.pos = round(x*data./sum(data));
baric.time = ti(baric.pos);

%% Convolution in time
% for i=1:nmeas
%     z1(:,i) = conv(y(:,i),irfi);
% end
% z = z1(1:numel(ty),:);
% %% resample in the forward timescale
% %z = interp1(tirf,z1,ty);
% if verbosity == 1
%     figure;
%     semilogy(z1);%(ty,z);%[y,z]);
%     figure;
%     semilogy(ty,[z(:,1),y(:,1)],tirfi,irfi),legend('conv','fwd','irf');
% end


