function J = toastJacobianTimedomain_INPROGRESS(hMesh,hBasis,qvec,mvec, dmask, mua, mus, ref,dt,nstep,twin,irf, GRADIENT)
% Time domain Jacobian for absorption c and K in  basis space
% A. Farina CNR - Polimi 09/2016
% M. Betcke
% S. Arridge
% possibility to use time-domain convolution or FFT
% CALC =
%   'conv'      : use time-domain convolution (less memory, slow)
%   'fft'       : use fft for convolution (much memory,  fast)
% possibility to calculate the gradient in 3 ways
% GRADIENT =
%   'toast'     : use toastImageGradient
%   'operator'  : use gradientOperator function by UCL
%   'matlab'    : use gradient function in Matlab (only 3D)
%   'none'      : avoid calcualation of pdmf_K
%
%
%% FLAG DEFINITION
CALC = 'fft';% 'fft';
USEGPU = 1;

if nargin < 13
    GRADIENT = 'operator';
end

if strcmpi(CALC,'fft')
    disp('TD Jacobian is calculated using fft');
else
    disp('TD Jacobian is calculated using direct convolution');
end
if strcmpi(GRADIENT,'none')
    disp('J_K part is not calculated')
end
        
%% Basis or mesh space
if hBasis == 0
    %    JMESH = 1;
    nsol = hMesh.NodeCount;
else
    %    JMESH = 0;
    solmask = (find(hBasis.GridElref>0));
    nsol = numel(solmask);
end

%%  vector lengths
nwin = size(twin,1);
nQM = sum(dmask(:));
NDIM = hMesh.Dimension;     %   2D or 3D
Nmesh = hMesh.NodeCount;

BB = hMesh.BoundingBox();
bmin = BB(1,:); bmax = BB(2,:);
bSize = hBasis.Dims;bSize = bSize';
vscale = (bmax-bmin)./(bSize-1);
blen = prod(bSize);
%% Internal fields calculation
theta = 1;%0.5;
disp('Computing direct and adjoint fields')
param = assignFields([], 'hMesh', hMesh, 'mua', mua, 'mus', mus, 'ref', ref);
timeSolver_('compute', [], dt, nstep, theta, param, USEGPU);
tphi = timeSolver_('apply', qvec);
taphi = timeSolver_('apply', mvec);

%% Convolution with IRF if present
% if numel(irf) > 1
%     %for j = 1:size(tphi)
%     tphi = convn(irf,tphi);
%     tphi(nstep + 1:end,:,:) = [];
% end
% Elegant but memory consuming
% if numel(irf)>1
%     lf = numel(irf) + 2*nstep -1;
%     Firf = repmat(fft(irf,lf),[1 size(tphi,2) size(tphi,3)]);
%     tphi = fft(tphi,lf,1).*Firf;
%     taphi = fft(taphi,lf,1);
%     clear Firf
% end

%% Allocate memory space for Jacobian
if strcmpi(GRADIENT,'none')
    J = zeros(nQM * nwin, nsol);
else
    J = zeros(nQM * nwin,2*nsol);
end
% Preallocate variables
if ~strcmpi(GRADIENT,'none')
    dim = zeros(nstep,blen);
    aim = zeros(nstep,blen);
    gd = zeros(nstep,NDIM,nsol);
    ga = zeros(size(gd));
    rhost = zeros(2*nstep-1,nsol);
end
if strcmpi(CALC,'conv')
    da = zeros(2*nstep-1,Nmesh);
else
    da = zeros(2*nstep-1,nsol);
end
rhoat = zeros(nstep,nsol);

switch lower(CALC)
%% CALCULATION USING DIRECT CONVOLUTION IN TIME
    case 'conv'
        row_off = 0;
        for q = 1:size(qvec,2)
            ind_m = find(dmask(:,q));
            disp(['Building Jacobian for source: ', num2str(q)]);
            
            for im = 1:length(ind_m)
                m = ind_m(im);
                
                disp(['Building Jacobian for source / measurement: ', num2str(q) '/' num2str(m)]);
                if ~strcmpi(GRADIENT,'none')
                    for it = 1:nstep
                        % Map internal fields to basis
                        dim(it,:) = hBasis.Map('M->B',tphi(it,:,q));
                        aim(it,:) = hBasis.Map('M->B',taphi(it,:,m));
                        % Map gradient field
                        dd = hBasis.ImageGradient(dim(it,:)); % is a vector!
                        gd(it,1:NDIM,:) = dd(:,solmask);
                        dd = hBasis.ImageGradient(aim(it,:));
                        ga(it,1:NDIM,:) = dd(:,solmask);
                    end
                    % loop on solution elements - Scattering
                    for iel = 1:nsol
                        for idim=1:NDIM
                            rhost(:,iel) = rhost(:,iel) + squeeze(conv(gd(:,idim,iel),ga(:,idim,iel)));
                        end
                    end
                    dm = findIndexfromSQ(dmask,q,m);
                    J(row_off + (1:nwin), nsol+(1:nsol)) = -WindowTPSF(ConvIRF(rhost,irf(:,dm)),twin(:,:,dm));
                end
                % loop on the mesh elements - Absorption
                for iel=1:size(qvec,1)
                    da(:,iel) = conv(tphi(:,iel,q),taphi(:,iel,m));
                end
                % loop on time steps - Absorption
                for it=1:nstep
                    rhoat(it,:) = hBasis.Map('M->S',da(it,:));
                end
                dm = findIndexfromSQ(dmask,q,m);
                J(row_off + (1:nwin), 1:nsol) = -WindowTPSF(ConvIRF(rhoat,irf(:,dm)),twin(:,:,dm));
                row_off = row_off + nwin;
                toc
            end
        end
        %% CALCULATION USING FFT IN THE TIME-DOMAIN FOR SPEED UP THE CONVOLUTION
    case 'fft'
        %Overwrite the tphi, aphi for memory efficiency
        tphi = fft(tphi,2*nstep-1,1);
        taphi = fft(taphi,2*nstep-1,1);
        row_off = 0;
        switch lower(GRADIENT)
            %% GRADIENT
            case 'none'
                for q = 1:size(qvec,2)
                    ind_m = find(dmask(:,q));
                    disp(['Building Jacobian for source: ', num2str(q)]);
                    for im = 1:length(ind_m)
                        m = ind_m(im);
                        %disp(['Building Jacobian for source / measurement: ', num2str(q) '/' num2str(m)]);
                        for k = 1:size(tphi,1)
                            rhoat(k,:) = - hBasis.Map('M->S',tphi(k,:,q).*taphi(k,:,m)); %absorbtion
                        end
                        rhoat = ifft(rhoat,[],1);
                        dm = findIndexfromSQ(dmask,q,m);
                        J(row_off + (1:nwin), 1:nsol) = WindowTPSF(ConvIRF(rhoat,irf(:,dm)),twin(:,:,dm));
                        row_off = row_off + nwin;
                    end
                end
                
            case {'toast','matlab'}
                for q = 1:size(qvec,2)
                    ind_m = find(dmask(:,q));
                    disp(['Building Jacobian for source: ', num2str(q)]);
                    for im = 1:length(ind_m)
                        m = ind_m(im);
                        %disp(['Building Jacobian for source / measurement: ', num2str(q) '/' num2str(m)]);
                        for k = 1:size(tphi,1)
                            dim = hBasis.Map('M->B',tphi(k,:,q));
                            aim = hBasis.Map('M->B',taphi(k,:,m));
                            if strcmpi(GRADIENT,'toast')
                                gr = - sum(hBasis.ImageGradient(dim).*...
                                    hBasis.ImageGradient(aim),1);
                            else
                                [Gadx,Gady,Gadz] = gradient(reshape(aim,bSize),...
                                    vscale(1),vscale(2),vscale(3));
                                [Gdx,Gdy,Gdz] = gradient(reshape(dim,bSize),...
                                    vscale(1),vscale(2),vscale(3));
                                gr = -(Gdx.*Gadx + Gdy.*Gady + Gdz.*Gadz);
                            end
                            rhost(k,:) = gr(solmask);
                            rhoat(k,:) = - hBasis.Map('M->S',tphi(k,:,q).*taphi(k,:,m)); %absorbtion
                        end
                        rhoat = ifft(rhoat,[],1);
                        rhost = ifft(rhost,[],1);
                        dm = findIndexfromSQ(dmask,q,m);
                        J(row_off + (1:nwin), 1:nsol) = WindowTPSF(ConvIRF(rhoat,irf(:,dm)),twin(:,:,dm));
                        J(row_off + (1:nwin), nsol+(1:nsol)) = WindowTPSF(ConvIRF(rhost,irf(:,dm)),twin(:,:,dm));
                        row_off = row_off + nwin;
                    end
                end
            case 'operator'
                % direct field map
                tphi = permute(tphi,[2 1 3]); %permute to have space x frequency(in time) x source
                dim = zeros([prod(bSize), size(tphi,2), size(tphi,3)]);
                %dim(:,:) = M2B*tphi(:,:);  no more possible without the explicit matrix
                for i = 1:size(tphi,2)*size(tphi,3)
                    dim(:,i) = hBasis.Map('M->B',tphi(:,i));
                end
                
                % adjoint field map
                taphi = permute(taphi,[2 1 3]); %permute to have space x frequency(in time) x detector
                aim = zeros([prod(bSize), size(taphi,2), size(taphi,3)]);
                %aim(:,:) = M2B*taphi(:,:);
                for i = 1:size(taphi,2)*size(taphi,3)
                    aim(:,i) = hBasis.Map('M->B',taphi(:,i));
                end
                % allocate matrix for gradient operator
                tmp = zeros(bSize); tmp(solmask) = 1;
                [Dy, Dx, Dz] = gradientOperator(bSize,vscale,ones(prod(bSize),1),'none',tmp);
                row_off = 0;
                for q = 1:size(qvec,2)
                    %ind_m = find(dmask(q,:));
                    ind_m = find(dmask(:,q));	% A. Farina: consistent with dmask order
                    
                    disp(['Building Jacobian for source: ', num2str(q)]);
                    for im = 1:length(ind_m)
                        m = ind_m(im);
                        
                        rhoat = - dim(solmask,:,q).*aim(solmask,:,m); %absorbtion
                        gr = - ( (Dy*dim(:,:,q)) .* (Dy*aim(:,:,m)) ...
                            +(Dx*dim(:,:,q)) .* (Dx*aim(:,:,m)) ...
                            +(Dz*dim(:,:,q)) .* (Dz*aim(:,:,m)) );
                        rhost = gr(solmask,:);
                        
                        rhoat = permute(rhoat,[2 1]);
                        rhoat = ifft(rhoat,[],1); % time domain PMDF for absorption
                        rhost = permute(rhost,[2 1]);
                        rhost = ifft(rhost,[],1); % time domain PMDF for diffusion
                        dm = findIndexfromSQ(dmask,q,m);
                        J(row_off + (1:nwin), 1:nsol)        = WindowTPSF(ConvIRF(rhoat,irf(:,dm)),twin(:,:,dm)); %toastMapMeshToSol(hBasis,rhoaw(w,:));
                        J(row_off + (1:nwin), nsol+(1:nsol)) = WindowTPSF(ConvIRF(rhost,irf(:,dm)),twin(:,:,dm)); %toastMapMeshToSol(hBasis,rhosw(w,:));
                        row_off = row_off + nwin;
                    end
                end
        end
end















