function ddisp(str, debug, PREFIX)
% DDISP Debug display, 'str' is displayed if 'debug' > 0
%   ddisp(str, debug)

%@author: Marta M. Betcke

if nargin < 3
  PREFIX = '$$';
end

if debug > 0
  disp([PREFIX ': ' str]);
end