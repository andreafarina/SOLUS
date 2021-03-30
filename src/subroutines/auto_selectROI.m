function out = auto_selectROI(x, p1,p2)
% retunrs the interval between two specified values of the curve x

if nargin < 2
   p1 = 0.90;
   p2 = 0.10;
end
    peak = max(x(:));
    out(1) = find(x>p1*peak, 1,'first');
    out(2) = find(x>p2*peak, 1,'last');
    if out(1)>out(2)
       error('Something wrong with the values') 
    end
    return
end