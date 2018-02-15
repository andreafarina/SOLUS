function [] = ShowTimeWindows(y,twin,dt)
% Plot Time windows overlayed to time-resolved data y depending on twin.
% if the time step dt is provided the plot has time on abscissa
figure
if nargin < 3
    dt = 1;
    strx = 'channel';
else
    strx = 'time (ps)';
end
t = (1:size(y,1))*dt;
semilogy(t,y),ylim([max(y(:))./1e3 max(y(:))]),
xlabel(strx),
ylabel('counts')%,set(gca,'FontSize',20);
for i = 1:size(twin,1)
    rectangle('Position',[double(dt*twin(i,1)),min(y(y(:)>0)),...
        double(dt*twin(i,2)-dt*twin(i,1)+1),double(max(y(:)))]);
    
end
end