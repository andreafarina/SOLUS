function Q = QuantifyDOT(REC,ref_true)
% =========================================================================
%%                            Quantify DOT 
% =========================================================================
%% Prepare variables
x = REC.grid.x;
y = REC.grid.y;
z = REC.grid.z;
dx = REC.grid.dx;
dy = REC.grid.dy;
dz = REC.grid.dz;
[X,Y,Z] = ndgrid(x,y,z);
%% subtract background
dmua_rec3d = reshape(REC.opt.bmua-REC.opt.muaB,REC.grid.dim); 
%dmua_rec3d = reshape(REC.opt.bmua,REC.grid.dim); 

%% Prepare Mask for gaussian fit
radius_max = 20; %max radius of region of interest;
if ref_true == 0
    temp = dmua_rec3d;
    temp(:,:,1:10) = 0;
    [mxm,mind] = max(temp(:));
%mind = find(dmua_rec3d == mxm);
[Xrec,Yrec,Zrec] = ind2sub(size(dmua_rec3d),mind);
Xrec = x(Xrec);
Yrec = y(Yrec);
Zrec = z(Zrec);
else
    Xrec = REC.opt.hete1.c(1);
    Yrec = REC.opt.hete1.c(2);
    Zrec = REC.opt.hete1.c(3);
    temp = dmua_rec3d;
    
end
V = zeros(size(dmua_rec3d));
nzm = sqrt((X-Xrec).^2+(Y-Yrec).^2+(Z-Zrec).^2)<radius_max;
V(nzm) = 1;
%V = 1;
dmua_rec3d = dmua_rec3d.*V;
if ~exist('mxm','var')
    mxm = max(dmua_rec3d(:));
end
%% Gaussian 3d estimation
[COM_REC,COV_REC,TOT_VOL_REC] = Gaussian3D(x,y,z,dx,dy,dz,dmua_rec3d,V);

%% compare with reference
if ref_true == 1
    dmua_true3d = REC.opt.Mua-REC.opt.muaB; %% subtract background
    dmua_true3d = dmua_true3d.*V;
    [COM_TRUE,COV_TRUE,TOT_VOL_TRUE] = Gaussian3D(x,y,z,dx,dy,dz,dmua_true3d,V);
end

%% Print report 

for iinc = 1:1
disp(['------- REPORT for inclusion ',num2str(iinc),' ----------']);
if ref_true == 1
    disp(['REF_Tot_Volume: ',num2str(TOT_VOL_TRUE(iinc)),' mm^2']);
end
disp(['REC_Tot_Volume: ',num2str(TOT_VOL_REC(iinc)),' mm^2']);
disp('--------------------------------------------');
if ref_true == 1
    disp(['REF_COM(x,y,z): ',num2str(COM_TRUE(iinc,:))]);
end
disp(['REC_COM(x,y,z): ',num2str(COM_REC(iinc,:))]);
if ref_true == 1
    disp(['X_COM_error: ',num2str(COM_REC(iinc,1)-COM_TRUE(iinc,1))]);
    disp(['Y_COM_error: ',num2str(COM_REC(iinc,2)-COM_TRUE(iinc,2))]);
    disp(['Z_COM_error: ',num2str(COM_REC(iinc,3)-COM_TRUE(iinc,3))]);
    disp(['total COM_error: ',num2str(norm(COM_TRUE(iinc,:)-COM_REC(iinc,:)))]);
end
disp('--------------------------------------------');
if ref_true == 1
    disp('REF_COVARIANCE: ');
%disp(num2str(COV_TRUE(:,:,iinc),3));
    COV_TRUE(:,:,iinc)
end
disp('REC_COVARIANCE: ');
%disp(num2str(COV_REC(:,:,iinc)));
COV_REC(:,:,iinc)
%disp('------------------------------------------');
end


%% Centre of Mass
% COM_TRUE = zeros(3,1);
% vv = sum(dmua_true3d(:));
% COM_TRUE(1) = sum(dmua_true3d(:).*X(:))/vv;
% COM_TRUE(2) = sum(dmua_true3d(:).*Y(:))/vv;
% COM_TRUE(3) = sum(dmua_true3d(:).*Z(:))/vv;
% disp(['TRUE centre of mass [',num2str(COM_TRUE(1)),',',num2str(COM_TRUE(2)),',',num2str(COM_TRUE(3)),']']);
%  
% COM_REC = zeros(3,1);
% vr = sum(dmua_rec3d(:));
% COM_REC(1) = sum(dmua_rec3d(:).*X(:))/vr;
% COM_REC(2) = sum(dmua_rec3d(:).*Y(:))/vr;
% COM_REC(3) = sum(dmua_rec3d(:).*Z(:))/vr;
% disp(['REC centre of mass [',num2str(COM_REC(1)),',',num2str(COM_REC(2)),',',num2str(COM_REC(3)),']']);
%  
% COM_err = norm(COM_REC-COM_TRUE);
% disp(['Difference in centre of mass [',num2str(COM_REC(1)-COM_TRUE(1)),',',num2str(COM_REC(2)-COM_TRUE(2)),',',num2str(COM_REC(3)-COM_TRUE(3)),']']);
% disp(['Error in centre of mass ',num2str(COM_err)]);
% %F(errXP).Value(iL)=COM_REC(1)-COM_TRUE(1);
% %F(errYP).Value(iL)=COM_REC(2)-COM_TRUE(2);
% %F(errZP).Value(iL)=COM_REC(3)-COM_TRUE(3);
% 
% COM_REC = zeros(3,1);
%  vr = sum(V(:).*dmua_rec3d(:));
%  COM_REC(1) = sum(V(:).*dmua_rec3d(:).*X(:))/vr;
%  COM_REC(2) = sum(V(:).*dmua_rec3d(:).*Y(:))/vr;
%  COM_REC(3) = sum(V(:).*dmua_rec3d(:).*Z(:))/vr;
% COM_err = norm(COM_REC-COM_TRUE);
% disp(['Difference in centre of mass in region [',num2str(COM_REC(1)-COM_TRUE(1)),',',num2str(COM_REC(2)-COM_TRUE(2)),',',num2str(COM_REC(3)-COM_TRUE(3)),']']);
% disp(['Error in centre of mass in region ',num2str(COM_err)]);
% 
% %% covariance of mass
% VAR_TRUE = zeros(3,1);
% VAR_TRUE(1) = sum(dmua_true3d(:).*((X(:)-COM_TRUE(1)).^2 ))/vv;
% VAR_TRUE(2) = sum(dmua_true3d(:).*((Y(:)-COM_TRUE(2)).^2 ))/vv;
% VAR_TRUE(3) = sum(dmua_true3d(:).*((Z(:)-COM_TRUE(3)).^2 ))/vv;
% 
% VAR_REC = zeros(3,1);
% vr = sum(dmua_rec3d(:));
% VAR_REC(1) = sum(dmua_rec3d(:).*((X(:)-COM_REC(1)).^2 ))/vr;
% VAR_REC(2) = sum(dmua_rec3d(:).*((Y(:)-COM_REC(2)).^2 ))/vr;
% VAR_REC(3) = sum(dmua_rec3d(:).*((Z(:)-COM_REC(3)).^2 ))/vr;
% 
% disp(['Sqrt of covariance of mass [',num2str(sqrt(VAR_REC(1))),',',num2str(sqrt(VAR_REC(2))),',',num2str(sqrt(VAR_REC(3))),']']);
% 
% vr = sum(V(:).*dmua_rec3d(:));
% VAR_REC(1) = sum(V(:).*dmua_rec3d(:).*((X(:)-COM_REC(1)).^2 ))/vr;
% VAR_REC(2) = sum(V(:).*dmua_rec3d(:).*((Y(:)-COM_REC(2)).^2 ))/vr;
% VAR_REC(3) = sum(V(:).*dmua_rec3d(:).*((Z(:)-COM_REC(3)).^2 ))/vr;
% disp(['Sqrt of covariance of mass in region [',num2str(sqrt(VAR_REC(1))),',',num2str(sqrt(VAR_REC(2))),',',num2str(sqrt(VAR_REC(3))),']']);

%% Sigma and CNR
sigma=std(dmua_rec3d(:));
disp(['Standard deviation ',num2str(sigma)]);

av=max(dmua_rec3d(:));
disp(['Max ',num2str(av)]);
disp(['CNR ',num2str(av./sigma)]);

% av=max(V(:).*dmua_rec3d(:));
% disp(['Max in region ',num2str(av)]);
% disp(['CNR in region ',num2str(av./sigma)]);
% 
% 
% %% Quantitation
%  factor=REC.opt.hete1.sigma^3/(REC.grid.dx*REC.grid.dy*REC.grid.dz)*(4/3*pi);
%  DmuaVolTrue=sum(dmua_true3d(:))/factor;
%  DmuaVolRec=sum(dmua_rec3d(:))/factor;
%  disp(['DmuaVolTrue ',num2str(DmuaVolTrue)]);
%  disp(['DmuaVolRec ',num2str(DmuaVolRec)]);
%  disp(['error ',num2str(DmuaVolRec/DmuaVolTrue-1)]);
% %F(relMUAP).Value(iL)=(DmuaVolRec-DmuaVolTrue)/DmuaVolTrue;
%% prepare mask for CNR (calculated as Arridge, Ducros..)
W = zeros(size(dmua_rec3d));
W(temp>mxm/2) = 1;
Wback = 1 - W; 
dmuaROI = W(:).*dmua_rec3d(:);
dmuaROI(W==0) = [];
muROI = mean(dmuaROI);
dmuaBack = Wback(:).*dmua_rec3d(:);
dmuaBack(Wback==0) = [];
muBack = mean(dmuaBack);

CNR2 =(muROI - muBack)./...
    sqrt(sum(W(:))./numel(W).*(std(dmuaROI(:)).^2) + ...
     sum(Wback(:))./numel(W).*(std(dmuaBack(:)).^2));

Q.COM.rec = COM_REC;
Q.COV.rec = COV_REC;
Q.max.rec = mxm;
Q.volume.rec = TOT_VOL_REC;
Q.cnr = max(dmua_rec3d(:))./std(dmua_rec3d(:));
Q.cnr2 = CNR2;
disp(['CNR2 = ',num2str(CNR2)]);
if ref_true == 1
    Q.COM.true = COM_TRUE;
    Q.COM.error = COM_REC - COM_TRUE;
    Q.COV.true = COV_TRUE;
    
    Q.max.true = max(dmua_true3d(:));
    Q.max.error = mxm -Q.max.true;
    Q.max.rel_error = Q.max.error./Q.max.true;
    
    disp(['Relative max error = ',num2str(Q.max.rel_error)]);
    
    Q.volume.true = TOT_VOL_TRUE;
    Q.volume.rel_error = (TOT_VOL_REC - TOT_VOL_TRUE)./TOT_VOL_TRUE;
    disp(['Relative volume error = ',num2str(Q.volume.rel_error)]);
end