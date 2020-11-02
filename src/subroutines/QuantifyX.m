function Q = QuantifyX(grid,xB,xrec,ref_true,c_true,xtrue,FullPrint)
%% Prepare variables
x = grid.x;
y = grid.y;
z = grid.z;
dx = grid.dx;
dy = grid.dy;
dz = grid.dz;
dV = dx*dy*dz;
[X,Y,Z] = ndgrid(x,y,z);
%% subtract background
dx_rec3d = reshape(xrec-xB,grid.dim);
%dmua_rec3d = reshape(xrec,REC.grid.dim);

%% Prepare Mask for gaussian fit
radius_max = 20; %max radius of region of interest;
if ref_true == 0
    temp = dx_rec3d;
    %temp(:,:,1:10) = 0;
    [mxm,mind] = max(temp(:));
    %mind = find(dmua_rec3d == mxm);
    [Xrec,Yrec,Zrec] = ind2sub(size(dx_rec3d),mind);
    Xrec = x(Xrec);
    Yrec = y(Yrec);
    Zrec = z(Zrec);
else
    Xrec = c_true(1);
    Yrec = c_true(2);
    Zrec = c_true(3);
    temp = dx_rec3d;
    
end
V = zeros(size(dx_rec3d));
nzm = sqrt((X-Xrec).^2+(Y-Yrec).^2+(Z-Zrec).^2)<radius_max;
V(nzm) = 1;
%V = ones(size(dmua_rec3d));
dx_rec3d = dx_rec3d.*V;
if ~exist('mxm','var')
    mxm = max(dx_rec3d(:));
end
%% Gaussian 3d estimation
fit_type = 'normal'; % 'normal' or 'gaussian'
try
    [COM_REC,COV_REC,TOT_VOL_REC] = Gaussian3D(x,y,z,dx,dy,dz,dx_rec3d,V,fit_type);
    
    %% compare with reference
    if ref_true == 1
        dx_true3d = xtrue-xB; %% subtract background
        dx_true3d = dx_true3d.*V;
        [COM_TRUE,COV_TRUE,TOT_VOL_TRUE] = Gaussian3D(x,y,z,dx,dy,dz,dx_true3d,V,fit_type);
    end
catch ME
    if contains(ME.message,'lsqcurvefit')
        disp('---- Error in fitting procedure ----')
        Q = [];
        return
    else
        throw(ME);
    end
end

%% Print report
if FullPrint == true
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
    
    %% Sigma and CNR
    sigma=std(dx_rec3d(:));
    disp(['Standard deviation ',num2str(sigma)]);
    
    av=max(dx_rec3d(:));
    disp(['Max ',num2str(av)]);
    disp(['CNR ',num2str(av./sigma)]);
    
    av=max(V(:).*dx_rec3d(:));
    disp(['Max in region ',num2str(av)]);
    disp(['CNR in region ',num2str(av./sigma)]);
    
end
%% Quantitation
%factor = REC.opt.hete1.sigma^3/(REC.grid.dx*REC.grid.dy*REC.grid.dz)*(4/3*pi);
if ref_true == 1
    DxVolTrue = sum(dx_true3d(:)*dV);%/factor;
end
DxVolRec = sum(dx_rec3d(:)*dV);%/factor;
%disp(['DmuaVolTrue ',num2str(DmuaVolTrue)]);
%disp(['DmuaVolRec ',num2str(DmuaVolRec)]);
%disp(['error ',num2str(DmuaVolRec/DmuaVolTrue-1)]);
%F(relMUAP).Value(iL)=(DmuaVolRec-DmuaVolTrue)/DmuaVolTrue;
%% prepare mask for CNR (calculated as Arridge, Ducros..)
W = zeros(size(dx_rec3d));
W(temp>mxm/2) = 1;
Wback = 1 - W;
dxROI = W(:).*dx_rec3d(:);
dxROI(W==0) = [];
muROI = mean(dxROI);
dxBack = Wback(:).*dx_rec3d(:);
dxBack(Wback==0) = [];
muBack = mean(dxBack);

CNR2 =(muROI - muBack)./...
    sqrt(sum(W(:))./numel(W).*(std(dxROI(:)).^2) + ...
    sum(Wback(:))./numel(W).*(std(dxBack(:)).^2));

CNR3 =(muROI - muBack)./...
    sqrt(sum(Wback(:))./numel(W).*(std(dxBack(:)).^2));

Q.COM.rec = COM_REC;
Q.COV.rec = COV_REC;
Q.max.rec = mxm;
Q.volumeG.rec = TOT_VOL_REC;
Q.volume.rec = DxVolRec;
Q.cnr = max(dx_rec3d(:))./std(dx_rec3d(:));
Q.cnr2 = CNR2;
Q.cnr3 = CNR3;
if FullPrint == true
    disp(['CNR2 = ',num2str(CNR2)]);
end
if ref_true == 1
    Q.COM.true = COM_TRUE;
    Q.COM.error = COM_REC - COM_TRUE;
    Q.COV.true = COV_TRUE;
    
    Q.max.true = max(dx_true3d(:));
    Q.max.error = mxm -Q.max.true;
    Q.max.rel_error = Q.max.error./Q.max.true;
    
    
    if strcmpi(fit_type,'gaussian')
        Q.volumeG.true = TOT_VOL_TRUE;
        Q.volumeG.rel_error = (TOT_VOL_REC - TOT_VOL_TRUE)./TOT_VOL_TRUE;
        if FullPrint == true
            disp(['Relative volume (gaussian fit) error = ',num2str(Q.volumeG.rel_error)]);
        end
    end
    Q.volume.true = DxVolTrue;
    Q.volume.rel_error = (DxVolRec-DxVolTrue)/DxVolTrue;
    Q.total.rel_error = norm(dx_true3d(:)-dx_rec3d(:))./norm(xrec+dx_true3d(:));
    if FullPrint == true
        disp(['Relative max error = ',num2str(Q.max.rel_error)]);
        disp(['Relative volume error = ',num2str(Q.volume.rel_error)]);
        disp(['Integral relative error = ',num2str(Q.total.rel_error)]);
    end
end