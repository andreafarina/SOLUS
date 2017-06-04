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

ti=t(1):dt:t(end); ti=ti';
yi = interp1(t,y,ti);

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

