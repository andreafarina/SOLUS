function out = auto_selectROI(x, p1,p2)
% returns the interval between two specified values of the curve x

% x = linspace(1,0,360) .* sin(linspace(0,2*pi,360)).^(2);

if nargin < 2
%    p1 = 0.50;
%    p2 = 0.0005;
   p1 = 0.1;
   p2 = 0.01;%005;

end
    peak = max(x(:));
    out(1) = find((x)>p1*peak, 1,'first');
    out(2) = find((x)>p2*peak, 1,'last');
    return
    % handle local maxima
    [pks,locs] = findpeaks(x/max(x(:)),'MinPeakProminence',p2);%2*std(x/max(x(:))));
    [pks, I] = sort( pks,'descend');
    locs = locs(I);
    % tail

    d = 1-pks; 
    
    if numel(d) ~= 1
        d = d(1:2);
        locs = locs(1:2);
        if any(d(2:end)<(1-p2) & (locs(2:end)>= locs(1)) & (locs(2:end)<=out(2)) )
            out(2) = out(2)-1;            
        end

        if any(d(2:end)<(1-p1) & locs(2:end)<= locs(1) & locs(2:end)>=out(1))
            out(1) = out(1)+1;
        end
        Imask = 0*x;
        Imask(out(1):out(2))=1;
        out = auto_selectROI(x.*Imask, p1, p2);
    else
        return
    end
    if out(1)>out(2)
       error('Something wrong with the values') 
    end
%     rect = 0*x;
%     rect(out(1):out(2)) = peak;
%     figure(1),plot(x), hold on, plot(rect)
%     drawnow
%     pause(0.5)
    return
end