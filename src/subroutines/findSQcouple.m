function s_out = findSQcouple(dmask, idx, FLAG)
% retrieves the index on the mask for a given measurements
% complemtary of findMeasIndex


%dmask = dmask(:);

I = find(dmask(:) ~=0);
s_out = I(idx);
if nargin < 3
    return
end
if strcmpi(FLAG, 'sub')
    [s_out(1),s_out(2),s_out(3)] = ind2sub(size(dmask), s_out );
end


end