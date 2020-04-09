function sgm = snake_fitting(im, cor)
% This function takes as input a 2d rgb US-image and a set of user generated
% points and tries to fit a better contour to the inclusion


%load to_load4
im = double(rgb2gray(im));
im = im/max(im(:));
cor = cor(:, 1:end);
cor0 = cor;
npoints = 400;
nparam = 2*(numel(cor)-1);
points0 = Lee_spline(cor, npoints);
param = Lee_spline(cor, nparam+1);
param = param(:,1:end-1);
param0 = param;

order = 1;
%options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
%param = fminunc(@Loss2, param0, options);
dStep = 0.15; 

sigma_p = linspace(0.6,0.15,5);
i_s = 1;
dstruct = setDomain(im, param0, sigma_p, points0);
dstruct0 = dstruct;
dstruct.dist = dstruct0.dist(:,:,i_s);
dstruct.points0 = points0;
k = 1;
max_k = 10;
while k <= max_k
 %   drawfitting(param, npoints, im);
    
    L = Loss(forward(param,npoints), dstruct,param, param0);
    LL(k) = L;
    [dL, d2L, dUp] = computeGradient2(param, dStep, npoints, dstruct, param0,L, order);
    alpha = (10* dStep);%/max(dUp(:));
    ik = 1;
    Ln = L + 1;% Loss(forward((param - alpha*sign(dUp)),npoints), dstruct,param, param0);
    while Ln >= L && ik <=300 
        alpha = alpha * 0.8;
        Ln = Loss(forward((param - alpha*sign(dUp)),npoints), dstruct,param, param0);
        ik = ik + 1;
    end
    if Ln < L
        param = param - alpha*sign(dUp);
    else
        k = max_k + 1;
        disp('breaking cycle');
    end
    
    drawfitting(param, npoints, im);
    k = k + 1;
    if k>=max_k && i_s < numel(sigma_p)
        i_s = i_s + 1;
        dstruct.dist = dstruct0.dist(:,:,i_s); %setDomain(im, param0, sigma_p(i_s), points0);
        dstruct.points0 = points0;
        k = 1;
        disp('New sigma')
    end
end


%L = Loss(forward(param,npoints), dstruct,param, param0);
points=forward(param, npoints)';
sgm=roipoly(im,points(:,1),points(:,2));

disp('end')
% 
% function [L, dL] = Loss2(param)
%     lpoints = forward(param, npoints);
%     d = dstruct;
%     idx = sub2ind(size(d.dist),round(lpoints(1,:)),round(lpoints(2,:)));
%     eim = d.dist(idx);
%     Eim = sum(eim(:));  
% 
%     regu = param - param0;
%     Regu = sum(regu(1,:).^2 + regu(2,:).^2 );
%     
%     Eint = 0;
%     Eext = 0;
%     
%     L = Eint + Eext + Eim + 0 * Regu;
%     dStep = 0.15;
%     order = 1;
%     [dL,~, ~] = computeGradient2(param, dStep, npoints, dstruct, param0, L, order);
%     dL = dL/dStep;
% end
% 
% 
 end
% 

function drawfitting(param, npoints, im)
    points = forward(param, npoints);
    figure(1), imagesc(im);hold on
    plot(points(1,:),points(2,:),'r','LineWidth',2);
    plot(cat(2,param(1,:), param(1,1)),cat(2,param(2,:), param(2,1)),'y.');
    drawnow;
    
end

function out = forward(spline_points, num)

    corf = cat(2, spline_points, spline_points(:,1));
    %cor = spline_points;
    out = Lee_spline(corf, num);

end



function d = setDomain(im,cor, sigma_p, pp)
    [gx,gy] = imgradientxy(im);
    [gx2,~] = imgradientxy(gx);
    [~, gy2] = imgradientxy(gy);
    lap = gx2 + gy2;
    
    minx = min(cor(1,:));
    maxx = max(cor(1,:));
    miny = min(cor(2,:));
    maxy = max(cor(2,:));
    semix = (maxx-minx)/2;
    semiy = (maxy-miny)/2;
    minx = floor(minx - 0.2 * min(semix,semiy));
    maxx = ceil(maxx + 0.2 * min(semix,semiy));
    miny = floor(miny - 0.2 * min(semix,semiy));
    maxy = ceil(maxy + 0.2 * min(semix,semiy));

    [xx, yy] = meshgrid(1:size(lap',1), 1:size(lap',2));
    roi_idx = zeros(size(lap));
    roi_idx(xx>=minx & xx<=maxx & yy>= miny & yy <= maxy) = 1;
    roi_idx = logical(roi_idx(:));

    
    [sigma, add_dist] = findRadius(pp,im);
    d.dist = zeros(size(im,1), size(im,2), numel(sigma_p));
    for i_sigma = 1: numel(sigma_p)
        smlap = zeros(size(lap));    
        G = mygaussian(sigma_p(i_sigma)*sigma, size(im));
        smlap_ = conv2(lap, G, 'same');
        smlap = smlap+ smlap_;%*sum(abs(smlap_(:))); 
        smlap = abs(smlap);
        prc1 = 0;
        prc2 = 10;
        dist = zeros(size(smlap));
        i = prctile(smlap(roi_idx), prc1);
        j = prctile(smlap(roi_idx), prc2);
        dist = dist + double( bwdist(smlap >= i & smlap <= j ));%+ bwdist(smlap <= i).^2 ); 
        d.smlap = smlap;
        dist = dist - min(dist(:));
        dist = dist/ max(dist(:));
        add_dist = conv2(add_dist,mygaussian(1, size(add_dist)), 'same');
        d.dist(:,:,i_sigma) = dist.^2+ (add_dist*(max(dist(roi_idx))/max(add_dist(:)))).^2;
    end
    
    d.lap = lap;
%    imagesc(sum(d.dist(:,:,3)));

end

function g = mygaussian(sigma, siz)

    x = linspace(-0.5*siz(1),0.5*siz(1)  ,siz(1));
    y = linspace(-0.5*siz(2),0.5*siz(2)  ,siz(2));
    [xx, yy] = meshgrid(x,y);
    g = exp( - 0.5 * ((xx.^2 +yy.^2)/sigma.^2));
    g = g/ sum(g(:));
    
end

function [rmax,out_dist] = findRadius(pp,ima)
    points=pp';
    sgm=roipoly(ima,points(:,1),points(:,2))';
    sgm=sgm';

    dist_transform = double(bwdist(logical(1 - sgm)));
    rmax = max(dist_transform(dist_transform >0));
    out_dist = double(dist_transform)/max(double(dist_transform(:)));

end

function dL = computeGradient(param, dStep, npoints, domain, cor0, L)

%        L = computeLoss(points, domain,param, cor0);
        dL = 0 * param;
        for i = 1:size(param,1)
            for j = 1:size(param,2)
                dparam = 0*param;
                dparam(i,j) = dStep;
                dparam = param + dparam;
                dpoints = forward(dparam, npoints);
                dL(i,j) = Loss(dpoints, domain,dparam, cor0);
            end
        end
        dL = dL - L;
        dL = dL(:)/dStep;
end




function  [dL,d2L, dUp] = computeGradient2(param, dStep, npoints, domain, cor0, L, order)

%        L = computeLoss(points, domain,param, cor0);
        L = L(:);
        dL = 0 * param(:);
        d2L = zeros(numel(dL),numel(dL));
        for i = 1:numel(dL)
            dparam = 0*param(:);
            dparam(i) = dStep;
            dparam = param(:) + dparam;
            dpoints = forward(reshape(dparam, size(param)), npoints);
            dL(i) = Loss(dpoints, domain,reshape(dparam, size(param)), cor0) - L;
            dL(i) = dL(i)/dStep;
            if order ==2
                for j = 1:numel(dL)
                    ddparam = 0*param(:);
                    ddparam(j) = dStep;
                    ddparam = param(:) + dparam + ddparam;
                    ddpoints = forward(reshape(ddparam, size(param)), npoints);
                    d2L(i,j) = 0.5 * (Loss(ddpoints, domain,reshape(ddparam, size(param)), cor0) - dL(i));
                    d2L(i,j) = d2L(i,j)/dStep; 
                end
            end
        end
        if order ==2
            dUp = reshape(d2L'\dL, size(param)); 
            return;
        else 
            dUp = reshape(dL, size(param));
            d2L = 0;
            return;
        end
end


function L = Loss(points, d,param, cor0)
    Eim = 0;
    idx = sub2ind(size(d.dist),round(points(2,:)),round(points(1,:)));
    idx0 = sub2ind(size(d.dist),round(d.points0(2,:)),round(d.points0(1,:)));
    for i_s = 1:size(d.dist,3)
        dist_ = d.dist(:,:,i_s);
        eim = dist_(idx);
        Eim = Eim + sum(eim(:));  
    end
    Eim = Eim / size(d.dist,3);
    regu = param - cor0;
    Regu = sum(regu(1,:).^2 + regu(2,:).^2 );
    
    a = 0;
    b = 2;
    ab = a + b;
    a = a / ab;
    b = b / ab;
    eint0 = circshift(d.points0,1,2) - d.points0;
    eint0 = sum( a *( eint0(1,:).^2 + eint0(2,:).^2)/numel(idx0) );
    eint1 = (circshift(points,1,2) - points) /numel(idx); 
    eint1 = a*(eint1(1,:).^2 + eint1(2,:).^2);
    eint2 = 0.5 * (circshift(points,1,2) - 2 *  points + circshift(points,-1,2));
    eint2 =b * ( eint1(1,:).^2 + eint2(2,:).^2);      
    Eint = (sum(-eint1 + eint2)) - 0*eint0;
    
    
    centre = mean(cor0(:,:),2);
    dist2 = sqrt((cor0(1,:)  - centre(1)).^2 + (cor0(2,:)  - centre(2)).^2);
    Eext = 1/sum(dist2);
    
    L = 0.05*Eint + 0 * Eext + Eim + 0 * Regu;
  
end

