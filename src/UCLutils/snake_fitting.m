function [sgm, param] = snake_fitting(im, cor)
%
% [SGM, PARAM] = snake_fitting(IM, COR)
% This function takes as input a 2d rgb US-image IM and a set of user generated
% points COR and tries to fit a better contour to the inclusion.
% It returns the resulted segmented image SGM and an updated set of points
% PARAM to generate the spline if needed.


%load to_load3 % load for safety test
% convert image to double and normailise
im = double(rgb2gray(im));
im = im/max(im(:));
% keep user selected points
cor0 = cor;
% number of points after spline
npoints = 400;
% number of parameters though which the updated spline will be run
nparam = 2*(numel(cor)-1);
points0 = Lee_spline(cor, npoints);
param = Lee_spline(cor, nparam+1);
param = param(:,1:end-1);
param0 = param;

order = 1;
dStep = 0.1; 
%drawfitting(param, npoints, im);
nsigma = 15;
minsigma = 0.15;
maxsigma = 0.7;
sigma_p = linspace(maxsigma,minsigma,nsigma);
i_s = 1;
% calculate image features
dstruct = setDomain(im, param0, sigma_p, points0);
dstruct0 = dstruct;
dstruct0dist = dstruct0.dist;
dstruct.dist = dstruct0.dist(:,:,i_s);
dstruct.points0 = points0;

% initialise cycle
k = 1;
max_k = 16;
max_ik = 100;
alpha0 = 30*dStep;
% compute initial loss
Ln = Loss(forward(param,npoints), dstruct,param, param0);
while k <= max_k
    % update loss
    L = Ln;
    % compute gradient
    [~, ~, dUp] = computeGradient2(param, dStep, npoints, dstruct, param0, L, order);
    dUp = dUp/max(abs(dUp(:)));
    % update parameters
    alpha = alpha0;
    ik = 1;
    % compute loss
    Ln = Loss(forward((param - alpha*(dUp)),npoints), dstruct,param, param0);
    while Ln > L && ik <= max_ik 
        alpha = alpha * 0.9;
        %dUp = dUp;% .* (abs(dUp/L) > 0 );
        Ln = Loss(forward((param - alpha*(dUp)),npoints), dstruct,param, param0);
        ik = ik + 1;
    end
    if Ln < L % check if line search went well
        param = param - alpha*(dUp); 
    else % end minimisation for current smoothing
        k = max_k +1;
    end
    if  k>=max_k && i_s < numel(sigma_p) 
        % move to narrower smoothing
        i_s = i_s + 1;
        dstruct.dist = dstruct0dist(:,:,i_s);
        dstruct.points0 = points0;
        Ln = Loss(forward(param,npoints), dstruct,param, param0);
        k = 0;
    end
    k = k+1;
%    drawfitting(param, npoints, im);
end


drawfitting(param, npoints, im);
points=forward(param, npoints)';
sgm=roipoly(im,points(:,1),points(:,2));

end
function drawfitting(param, npoints, im)
    % draws 
    points = forward(param, npoints);
    figure(1), imagesc(im);colormap('gray'),hold on
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
    add_dist = conv2fft(add_dist,mygaussian(2, size(add_dist)));
    G = mygaussian(sigma_p.*sigma, size(im));
    d.dist = zeros(size(im,1), size(im,2), numel(sigma_p));
    prc1 = 0;
    prc2 = 10;
    % calculate smoothed images via fft2
    smlap3 = conv2fft(lap, G);
    for i_sigma = 1:numel(sigma_p)      
        %smlapc = abs(conv2(lap, G(:,:,i_sigma), 'same'));
        smlap = smlap3(:,:,i_sigma);
        i = prctile(smlap(roi_idx), prc1);
        j = prctile(smlap(roi_idx), prc2);
        dist = double( bwdist(smlap >= i & smlap <= j ));%+ bwdist(smlap <= i).^2 ); 
        dist = dist - min(dist(:));
        dist = dist/ max(dist(:));        
        d.dist(:,:,i_sigma) = dist.^2+ (add_dist*(max(dist(roi_idx))/max(add_dist(:)))).^2;
    end    
    d.lap = lap;
%    imagesc(sum(d.dist(:,:,3)));

end

function smlap3 = conv2fft(lap, G)
    
    padsize = ceil([size(lap,1), size(lap,2)]);
    lap_pad = padarray(repmat(lap, [1,1,size(G,3)]),padsize);    
    smlappad = abs(ifftshift(ifft2(fftshift(fft2(lap_pad)).*fftshift(fft2(padarray(G(:,:,:), padsize))))));
    smlap3 = smlappad(padsize(1):end - padsize(1)-1,padsize(2):end - padsize(2)-1,:);

end

function g = mygaussian(sigma, siz)
    
    g = zeros(siz(1),siz(2), numel(sigma));
    x = linspace(-0.5*siz(2),0.5*siz(2)  ,siz(2));
    y = linspace(-0.5*siz(1),0.5*siz(1)  ,siz(1));
    [xx, yy] = meshgrid(x,y);
    eucl = repmat(xx.^2 +yy.^2, [1,1, numel(sigma)]);
    sigma_ = zeros(1,1,numel(sigma));
    sigma_(1,1,:) = sigma;
    sigma = repmat(sigma_, [siz(1), siz(2), 1]);
    
    g = exp( - 0.5 * (eucl./sigma.^2));
    g = g ./ sum(sum(g, 1),2);
    
end

function [rmax,out_dist] = findRadius(pp,ima)
    points=pp';
    sgm=roipoly(ima,points(:,1),points(:,2))';
    sgm=sgm';
    dist_transform = double(bwdist(logical(1 - sgm)));
    rmax = max(dist_transform(dist_transform >0));
    out_dist = double(dist_transform)/max(double(dist_transform(:)));
end


function  [dL,d2L, dUp] = computeGradient2(param, dStep, npoints, domain, cor0, L, order)

        L = L(:);
        dL = 0 * param(:);
        if order ==2
                d2L = zeros(numel(dL),numel(dL));
        end
        
        sizparam = size(param);
        for i = 1:numel(dL)
            dparam = 0*param(:);
            dparam(i) = dStep;
            dparam = param(:) + dparam;
            dpoints = forward(reshape(dparam,sizparam), npoints);
            dL(i) = Loss(dpoints, domain,reshape(dparam,sizparam), cor0) - L;
            %dL(i) = dL(i)/dStep;
            if order ==2
                for j = 1:numel(dL)
                    ddparam = 0*param(:);
                    ddparam(j) = dStep;
                    ddparam = dparam + ddparam;
                    ddpoints = forward(reshape(ddparam, size(param)), npoints);
                    d2L(i,j) = 0.5 * (Loss(ddpoints, domain,reshape(ddparam, sizparam), cor0) - dL(i));
                    d2L(i,j) = d2L(i,j)/dStep; 
                end
            end
        end
        if order ==2
            dUp = reshape(d2L'\dL, size(param)); 
            return;
        else 
            dUp = reshape(dL/dStep, sizparam);
            d2L = 0;
            return;
        end
end


function L = Loss(points, d,~,~)%param, cor0)

    Eim = 0;
    idx = sub2ind(size(d.dist),round(points(2,:)),round(points(1,:)));
    %idx0 = sub2ind(size(d.dist),round(d.points0(2,:)),round(d.points0(1,:)));
    distim = d.dist;
    dsize3 = size(distim,3);
    for i_s = 1:dsize3
        dist_ = distim(:,:,i_s);
        eim = dist_(idx);
        Eim = Eim + sum(eim(:));  
    end
    Eim = Eim / dsize3;
%     regu = param - cor0;
%     Regu = sum(regu(1,:).^2 + regu(2,:).^2 );
    
    a = 0;
    b = 2;
    ab = a + b;
    a = a / ab;
    b = b / ab;
%     eint0 = circshift(d.points0,1,2) - d.points0;
%     eint0 = sum( a *( eint0(1,:).^2 + eint0(2,:).^2)/numel(idx0) );
    eint1 = (circshift(points,1,2) - points) /numel(idx); 
    eint1 = a*(eint1(1,:).^2 + eint1(2,:).^2);
    eint2 = 0.5 * (circshift(points,1,2) - 2 *  points + circshift(points,-1,2));
    eint2 =b * ( eint1(1,:).^2 + eint2(2,:).^2);      
    Eint = (sum(-eint1 + eint2)); % - 0*eint0;
        
%     centre = mean(cor0(:,:),2);
%     dist2 = sqrt((cor0(1,:)  - centre(1)).^2 + (cor0(2,:)  - centre(2)).^2);
%     Eext = 1/sum(dist2);    
    L = 0.05*Eint + Eim ;% 0 * Eext +0 * Regu;
  
end

