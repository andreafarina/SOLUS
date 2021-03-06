function [fbmua,fbmus,fbConc,fbA,fbbB]=FitVoxel(bmua,bmus,spe)
nL = spe.nLambda;
nV = size(bmua,1);
opts = optimoptions('lsqcurvefit',...
    'Jacobian','off',...
    ...'Algorithm','levenberg-marquardt',...
    'DerivativeCheck','off',...
    'MaxIter',100,'Display','none','FinDiffRelStep',repmat(1e-5,2,1));%,'TolFun',1e-10,'TolX',1e-10)

fbmua  = zeros(nV,nL); fbmus = zeros(nV,nL);
fbConc  = zeros(nV,spe.nCromo); fbA = zeros(nV,1); fbbB = zeros(nV,1); 
for iv = 1:nV
fbConc(iv,:) = spe.ext_coeff0\bmua(iv,:)';    
x = lsqcurvefit(@forward,ones(2,1),[],bmus(iv,:),[],[],opts);% x0 = x;
fbA(iv) = x(1);fbbB(iv) = x(2);
fbmus(iv,:) = x(1).*(spe.lambda./spe.lambda0).^(-x(2));
fbmua(iv,:) = spe.ext_coeff0*fbConc(iv,:)';
if rem(iv,1000)==0, display(['Voxel : ' num2str(iv) '/' num2str(nV)]); end
end
err_mus = (fbmus-bmus)./bmus *100;
err_mua = (fbmua-bmua)./bmua *100;
display(mean(err_mus,1));
display(mean(err_mua,1));
function [mus] = forward(x,~)
    
    mus = x(1).*(spe.lambda./spe.lambda0).^(-x(2));
%     mus = mus';
%     figure(345)
%     subplot(1,2,1),plot(spe.lambda,[data(1:nL) mua]), legend('data_mua','fit_mua')
%     subplot(1,2,2),plot(spe.lambda,[data(nL+1:end) mus]), legend('data_mus','fit_mus')
%     display(x')
end
end