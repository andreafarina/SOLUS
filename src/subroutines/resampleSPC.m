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
%%
verbosity = 0;
dty = t(2)-t(1);
bin_ratio = round(dt./dty,6);       
%% Better to first oversample if bin_ratio is not integer
if floor(bin_ratio)>(bin_ratio+1e-6) || floor(bin_ratio)<(bin_ratio-1e-6)
    dt_ov = 100e-3;    % 100 fs
    tii = t(1):dt_ov:t(end);
    yii = interp1(t,y,tii);
    area = sum(y);
    yii = yii./sum(yii)*area;
    bin_ratio = round(dt./dt_ov);
else
    yii = y;
    bin_ratio = floor(bin_ratio);
end 

ti = t(1):dt:(t(end)+0.5*dt);
bin = meshgrid(1:numel(ti),1:bin_ratio);
bin = bin(1:numel(yii));
bin = bin(:);
yi = accumarray(bin,yii);
yi(yi<0) = 0;
nyi=numel(yi);
nti=numel(ti);
if nyi < nti
    ti(nyi+1:end) = [];
elseif nyi > nti
    yi(nti+1:end) = [];
end
if verbosity == 1
    figure;semilogy(ti,yi./sum(yi),t,y./sum(y)),
    legend('resampled','original')
end
%% Calculate area
area = sum(yi);

%% normalize data
if sum(yi(:))==0 || isnan(sum(yi(:))) 
    yi(isnan(yi))=0;
    data = yi;
    area = 0; 
    peak =[];
    baric =[];
    return
end
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