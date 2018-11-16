
function out_curve = Lee_spline(points, npoints)
    % takes points as input and returns the spline interpolant using Eugene
    % Lee's centripetal method
    
    
    x = points(1,:);
    y = points(2,:);
    
    
    if points(:,1)==points(:,end)  
        endconds = 'periodic';
    else
        endconds = 'variational';
    end

    if length(x)==1 
        dt = 0;    
    else
        dt = sum((diff(points.').^2).'); 
    end

    t = cumsum([0,dt.^(1/4)]);
    
    sample_points = linspace(t(1), t(end), npoints);
    splinex = computeSpline(t,x, sample_points, endconds);
    spliney = computeSpline(t,y, sample_points, endconds);
    
    out_curve = [splinex';spliney'];
    
    
    

end
