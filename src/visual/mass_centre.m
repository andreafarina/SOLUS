function [pos]=mass_centre(x)

% [pos]=mass_centre(x)
% gets the indexes of the centre of mass "pos" of a N-dimensional matrix "x"

x = squeeze(x);
pos = ndims(x);
    for i = 1:ndims(x)
        ksum = 1:ndims(x);
        ksum(i) = [];
        sumax = x;
        for is = ksum
            sumax = sum(sumax, is);
        end
        ax = 1:size(x,i);
        pos(i) = sum(ax(:).*sumax(:)) / sum(sumax(:));
    end
end


