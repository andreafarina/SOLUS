function [mu,cov_matrix,area] = Gaussian3D(x,y,z,dx,dy,dz,data,V)
% =========================================================================
%                            Gaussian 3D
% Find center-of-mass (mu), covariance matrix(cov_matrix), and area of 
% a 3-fold variate gaussian distribution.
% size(data) = [numel(x),numel(y),numel(z)];
% =========================================================================
debug = 0;
%% check dimensions consistency
size_space = [numel(x),numel(y),numel(z)];
if sum(size_space-size(data))~=0
    disp('Attention! Check the dimensions of x,y,z and data!');
    return
end

%% mask again
mxm = max(data(:));
th = 0.01;
%V(data<mxm*th) = 0;
%data = data.*V;
%% Prepare space variables
[X,Y,Z] = ndgrid(x,y,z);
V = logical(V);
%V = logical(ones(size(V)));
%X = X(V);Y = Y(V);Z = Z(V);
space = [X(:), Y(:), Z(:)];

%data = data(V);
dV = dx*dy*dz;
%% Center of mass
data = data(:)';
mu = ((V(:).*data(:))' * space )./ sum(V(:).*data(:));
%posmu = round(mu./[dx,dy,dz]);

%% Covariance matrix
data3dd = repmat((V(:).*data(:))',[3 1]);
off_space = bsxfun(@minus,space,mu);
cov_matrix = ((data3dd.*off_space')*off_space)./sum(V(:).*data(:));
%% Area
area = sum(data(:).*V(:))*dV;
%% try a fit of the amplitude
%  G = @(A,b)(b+(exp(-1/2*sum(off_space/(cov_matrix).*off_space,2))*...
%      1./(sqrt((2*pi)^3*det(cov_matrix)))*A));
G = @(A,b,mu,cov_matrix)(b+(exp(-1/2*sum((space-mu)/(cov_matrix).*(space-mu),2))*...
   1./(sqrt((2*pi)^3*det(cov_matrix)))*A));
%ff = G(1,0);
%area = (ff(:)' * data(:))/(ff(:)'*ff(:));
% 
cov_matrix = eye(3,3);
%mu = [mean(x),mean(y),mean(z)];
%mu = [rand*max(x),rand*max(y),rand*max(z)];
 data = reshape(data,size_space);
opts = optimoptions('lsqcurvefit',...
    'Jacobian','off',...
    ...'Algorithm','levenberg-marquardt',...
    'DerivativeCheck','off',...
    'MaxIter',100,...
    ...'Display','iter-detailed',...
    'TolFun',1e-10);%,'TolX',1e-10)
 Gaussian = @(p,space)Gfun(p,space,size_space);
 p0 = [area,0,mu,cov_matrix(1,1),cov_matrix(1,2),cov_matrix(1,3),...
     cov_matrix(2,2),cov_matrix(2,3),cov_matrix(3,3)];

 p = lsqcurvefit(Gaussian,p0,space,data,[],[],opts);
 area = p(1);
 b = p(2);
 mu = [p(3),p(4),p(5)];
 cov_matrix = [p(6),p(7),p(8);
                p(7),p(9),p(10);
                p(8),p(10),p(11)];
 area = area + p(2)*dV*sum(V(:));           
if debug == 1
    Gfit = Gaussian(p,space);
    % slice on the inclusion
    iyslice = find(y == round(mu(2)));
    for i = 1:numel(z)
        data_slice = squeeze(data(:,iyslice,i));
        fit_slice = squeeze(Gfit(:,iyslice,i));
        figure(10000),
        plot([data_slice,fit_slice]),
        ylim([0 max(data(:))*1.2]), legend('data','fit');
        pause(0.2);
    end
    disp(['background = ',num2str(p(2))]);
end
            
end