function sgm = snake_fitting(im, cor)
% This function takes as input a 2d rgb US-image and a set of user generated
% points and tries to fit a better contour to the inclusion


%load to_load
im = double(rgb2gray(im));
im = im/max(im(:));
cor = cor(:, 1:end);
cor0 = cor;
npoints = 400;
nparam = 40;
points = Lee_spline(cor, npoints);
param = Lee_spline(cor, nparam+1);
param = param(:,1:end-1);
param0 = param;
dstruct = setDomain(im);
order = 1;
%options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
%param = fminunc(@Loss2, param0, options);
dStep = 1; 
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



function d = setDomain(im)
    [gx,gy] = imgradientxy(im);
    [gx2,~] = imgradientxy(gx);
    [~, gy2] = imgradientxy(gy);
    lap = gx2 + gy2;
    smlap = zeros(size(lap));
    for sigma = 20
        G = mygaussian(sigma, size(im));
        smlap = smlap+ conv2(lap, G, 'same');
    end
    smlap = (smlap)/max(abs(smlap(:)));
    
    dist = zeros(size(smlap));
    for i=prctile(smlap(:), 40):0.001:prctile(smlap(:), 60)
        dist = dist + double( bwdist(smlap >= i).^2+ bwdist(smlap <= i).^2 ); 
    end
    d.smlap = smlap;
    %dist = 1./(1+abs(smlap));
    d.dist = (dist /  max(dist(:))).^(1/2); %double(bwdist(smlap.* (smlap>prctile(smlap(:), 20))));% bwdist(abs(smlap).* (abs(smlap)>=1e-3)).^2; %dist - min(dist(:))) / (max(double(dist(:))) - min(dist(:)));    
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
    eim = d.dist(idx);
    Eim = sum(eim(:));  

    regu = param - cor0;
    Regu = sum(regu(1,:).^2 + regu(2,:).^2 );
    
    a = 1;
    b = 2;
    ab = a + b;
    a = a / ab;
    b = b / ab;
    eint1 = circshift(points,1,2) - points; 
    eint1 = (eint1(1,:).^2 + eint1(2,:).^2)/numel(idx);
    eint2 = 0.5 * (circshift(points,1,2) - 2 *  points + circshift(points,-1,2))/ numel(idx)^2;
    eint2 =b * ( eint1(1,:).^2 + eint2(2,:).^2);      
    Eint = sum(eint1 + eint2);
    
    
    centre = mean(cor0(:,:),2);
    dist2 = sqrt((cor0(1,:)  - centre(1)).^2 + (cor0(2,:)  - centre(2)).^2);
    Eext = 1/sum(dist2);
    
    L = Eint + 0 * Eext + Eim + 0 * Regu;
  
end

