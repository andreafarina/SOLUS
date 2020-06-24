function spline = computeSpline(t, y,  sample_points, endconds)
% Computes the spline function passing through the coordinates (t,y) and
% outputs the resulting function over the domain defined in sample_points
% endconds can be: 
% - "variational" for intepolating a natural spline
% - "periodic" for using periodic boundary conditions


    endcons = lower(endconds);
    h = t(2:end) - t(1:end -1); % dt
    b = (1./h) .* (y(2:end) - y(1:end - 1) ); % first derivative between points
    v0 = 2 * (h(end) + h(1));
    v = 2 * (h(1:end -1) + h(2:end));
    
    u0 = 6 * (b(1) - b(end)); 
    u = 6 * (b(2:end) - b(1:end - 1));
    
    if length(t) > 2   
        if strcmp(endcons,'periodic') == 1
            v = [v0,v,v0];
            u = [u0,u,u0];
            SIS = diag(h(1:end ), -1) + diag(v,0) + diag(h(1:end ),1);
            z = SIS \ u';
            z = [z;0];
        elseif strcmp(endcons,'variational') == 1
            SIS = diag(h(2:end - 1), -1) + diag(v,0) + diag(h(2:end - 1),1);
            z = SIS \ u';
            z = [0;z;0];
        else
            print('No method for spline selected. Using variational')
            SIS = diag(h(2:end - 1), -1) + diag(v,0) + diag(h(2:end - 1),1);
            z = SIS \ u;
            z = [0;z;0];
        end
    else
        z = [0,0]';
    end
    
    S = zeros(size(sample_points));
    %compute spline           
    for i = 1:length(h)
        idx = find(sample_points>= t(i) & sample_points <= t(i+1));
        S(idx) = (z(i+1)/(6*h(i)))*(sample_points(idx) - t(i)).^3 + ... 
                    (z(i)/(6*h(i)))*(t(i+1) - sample_points(idx)  ).^3 + ...
                    ( (y(i+1)/h(i)) - ( (z(i+1)/6)*h(i) ) ) * (sample_points(idx) - t(i)) + ...
                    ( (y(i)/h(i)) - ( (z(i)/6)*h(i) ) ) * (t(i+1) - sample_points(idx) );
    end

    spline = S(:);
    
end

