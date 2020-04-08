function sgm = snake_fitting(im, cor)
% This function takes as input a 2d rgb US-image and a set of user generated
% points and tries to fit a better contour to the inclusion


%load to_load
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

sigma_p = [0.5];
max_i_sigma = numel(sigma_p);
i_sigma = 1;
dstruct = setDomain(im, param0, sigma_p(i_sigma));
dstruct.points0 = points0;
k = 1;
max_k = 30;
while k <= max_k
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
    
    if k >=max_k && i_sigma < max_i_sigma
        k = 1;
        i_sigma = i_sigma + 1;
        dstruct = setDomain(im, param0, sigma_p(i_sigma));
        dstruct.points0 = points0;
        disp('refining sigma')
    end
        %param = param - dStep* (dUp/max(abs(dUp(:))));
    drawfitting(param, npoints, im);
    k = k + 1;
end

points=points';
if i>2
    sgm=roipoly(im,points(:,1),points(:,2));
else
    sgm=logical(zeros(size(im,1),size(im,2)));
end
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



function d = setDomain(im,cor, sigma_p)
    [gx,gy] = imgradientxy(im);
    [gx2,~] = imgradientxy(gx);
    [~, gy2] = imgradientxy(gy);
    lap = gx2 + gy2;
    
    
    if nargin == 1
        smlap = zeros(size(lap));
        for sigma = 10
            G = mygaussian(sigma, size(im));
            smlap = smlap+ conv2(lap, G, 'same');
        end
        smlap = (smlap)/max(abs(smlap(:)));
        prc1 = 40;
        prc2 = 60;
        roi_idx = 1:1:numel(smlap);
    elseif nargin ==3
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

        smlap = zeros(size(lap));
        for sigma = [sigma_p] * min(semix, semiy) 
            G = mygaussian(sigma, size(im));
            smlap_ = conv2(lap, G, 'same');
            smlap = smlap+ smlap_/sum(abs(smlap_(:))); 
        end
        smlap = (smlap)/max(abs(smlap(:)));
        
        [xx, yy] = meshgrid(1:size(smlap',1), 1:size(smlap',2));
        roi_idx = zeros(size(smlap));
        roi_idx(xx>=minx & xx<=maxx & yy>= miny & yy <= maxy) = 1;
        roi_idx = logical(roi_idx(:));
        prc1 = 30;
        prc2 = 70;
    end
        
    dist = zeros(size(smlap));
    for i=prctile(smlap(roi_idx), prc1):0.001:prctile(smlap(roi_idx), prc2)
        dist = dist + double( bwdist(smlap >= i).^2+ bwdist(smlap <= i).^2 ); 
    end
    d.smlap = smlap;
    dist = dist - min(dist(:));
    dist = dist/ max(dist(:));
    d.dist = dist;

end

function g = mygaussian(sigma, siz)

    x = linspace(-0.5*siz(1),0.5*siz(1)  ,siz(1));
    y = linspace(-0.5*siz(2),0.5*siz(2)  ,siz(2));
    [xx, yy] = meshgrid(x,y);
    g = exp( - 0.5 * ((xx.^2 +yy.^2)/sigma.^2));
    
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
    d.dist = d.dist;
    idx = sub2ind(size(d.dist),round(points(2,:)),round(points(1,:)));
    idx0 = sub2ind(size(d.dist),round(d.points0(2,:)),round(d.points0(1,:)));
    eim = d.dist(idx);
    Eim = sum(eim(:));  

    regu = param - cor0;
    Regu = sum(regu(1,:).^2 + regu(2,:).^2 );
    
    a = 1;
    b = 0.5;
    ab = a + b;
    a = a / ab;
    b = b / ab;
    eint0 = circshift(d.points0,1,2) - d.points0;
    eint0 = sum( a *( eint0(1,:).^2 + eint0(2,:).^2)/numel(idx0) );
    eint1 = (circshift(points,1,2) - points) /numel(idx); 
    eint1 = a*(eint1(1,:).^2 + eint1(2,:).^2);
    eint2 = 0.5 * (circshift(points,1,2) - 2 *  points + circshift(points,-1,2))/ numel(idx)^2;
    eint2 =b * ( eint1(1,:).^2 + eint2(2,:).^2);      
    Eint = (sum(eint1 + eint2)) - eint0;
    
    
    centre = mean(cor0(:,:),2);
    dist2 = sqrt((cor0(1,:)  - centre(1)).^2 + (cor0(2,:)  - centre(2)).^2);
    Eext = 1/sum(dist2);
    
    L = 0.0001*Eint + 0 * Eext + Eim + 0 * Regu;
  
end

