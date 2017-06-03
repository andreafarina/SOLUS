function winSize = tileFigs(monitor, winRatio)
% TILEFIGS Tiles figures across monitor number "monitor" of screen if no monitor specified  
%   @imput:
%   monitor: monitor number (e.g. if multiple monitors are used with one machine), 0: use entire screen size, 
%            as specified by get(0, 'ScreenSize')
%   winRatio: ratio of width to height of a window, default [4 3]
%   @output:
%   winSize: size of the tiled windows [width, height]

% Copyright 2011 Marta M. Betcke

if nargin < 1
  monitor = 0; %use all monitors, otherwise use monitor with the number "monitor"
end

%Get root object screen
screen = get(0);
screenSize = [screen.ScreenSize(3) screen.ScreenSize(4)];
menuOffset = [0, -2*27];
menuOffset = [0, 30];
monitors = get(0, 'MonitorPositions');

if monitor == 0
  monitorSize = [sum(monitors(:,3)) min(monitors(:, 4))];
  monitorOffset = [0 0];
else
  monitorSize = monitors(monitor, 3:4);
  monitorOffset = monitors(monitor, 1:2);
end
monitorSize = monitorSize + menuOffset;  

% archstr = computer('arch');
% switch archstr
%   case {'maci64', 'maci'}
%     screenOffset = [0, 0]; %modify for the top menue
%     screenSize = screenSize - [0,20];
%   case {'glnxa64', 'glnx86', 'pcwin', 'pcwin64'}
%     screenOffset = [0, 20]; %modify for the top menue
%     screenSize = screenSize - [20,0];
% end

if nargin < 2
  winRatio = [4 3];
end

%Find optimal tiling
nFigs = length(screen.Children);

swRatio = monitorSize./winRatio;
figRatio = swRatio(1)/swRatio(2);

n = sqrt(nFigs); m = n; 
if figRatio >= 1
  n = ceil(n);
  m = floor(m);
else
  n = floor(n);
  m = ceil(m);
end
while nFigs > n*m
  tn = (n+1)/m;
  tm = (m+1)/n;
  if abs(tn - swRatio) <= abs(tm - swRatio) 
    n = n+1;
  else
    m = m+1;
  end
end
winSize = floor([(monitorSize(1) - menuOffset(1))/n (monitorSize(2) - menuOffset(2))/m]);
winSize = winSize - menuOffset./[n m];
%Sort figures by ascending number
%set(0, 'Children', sort(get(0, 'Children')))
listChildren = sort(get(0, 'Children'));
%Distribute figures
 for jm = 1:m
   for jn = 1:n
     %j = (jm-1)*n+jn;
     j = (m-jm)*n+jn; 
    if j <= nFigs
      set(listChildren(j), 'OuterPosition', [-menuOffset(1)+monitorOffset(1)+(jn-1)*winSize(1)+1, -menuOffset(2)+monitorSize(2)-(m-jm+1)*winSize(2)+1, winSize(1), winSize(2)]);
      %set(listChildren(j), 'OuterPosition', [monitorOffset(1)+(jn-1)*winSize(1)+1, monitorSize(2)-(m-jm+1)*winSize(2)+1, winSize(1), winSize(2)]);
    else
      break
    end
  end
end
  