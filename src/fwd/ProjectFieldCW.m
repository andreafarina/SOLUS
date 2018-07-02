function [proj, proj_exc] = ProjectFieldCW(hMesh,qvec,mvec,dmask,...
            mua,musp,conc,n, freq,fwd,verbosity)
% calculate CW projection based on qvec, mvec, and dmask
% if fwd='diff' the ouptut is only diffusive fluence
% if fwd='fluo' the output is also the excitation

%spacer = ' ------- ';
%disp([spacer mfilename ': setup' spacer]);

tic;
N = hMesh.NodeCount();
nQ = size(dmask,2);       % effective sources
%nM = size(dmask,1);      % effective detectors
%nQM = sum(dmask(:));     % effective measurements
%==========================================================================
%%                             CW fluence (excitation)
%==========================================================================

%disp([spacer,'Excitation fluence',spacer]);
tic;

%S.phi = toastFields (S.mesh.hMesh, 0, S.qvec, S.mvec, S.opt.mua,...
%                    S.opt.musp, S.opt.n, S.freq);
Smat = dotSysmat2_noC(hMesh, mua, musp, n);
phi=Smat\qvec;

if verbosity==1
    figure;
    for i=1:nQ
        subplot(ceil(sqrt(nQ/1.5)),ceil(1.5*sqrt(nQ/1.5)),i);
        hMesh.Display(log(abs(phi(:,i))),'showcolorbar',0);
        %c = colorbar;
        %c.delete;
        %toastShowMesh(DOT.mesh.hMesh,toastMapSolToMesh(DOT.grid.hBasis,log(abs(DOT.phi))));
    end
    suptitle('Excitation fluence log(abs)');
    drawnow;
%     figure(4444)
%     subplot(1,3,1)
%     toastShowMesh(hMesh,log(abs(phi(:,1))));
%     subplot(1,3,2)
%     toastShowMesh(hMesh,log(abs(phi(:,3))));
%     subplot(1,3,3)
%     toastShowMesh(hMesh,log(abs(phi(:,7))));
    %saveas(gca,'/media/psf/Dropbox/UCL_Polimi/Final_Figures/pattern_fluence.pdf')% -transparent
end

%==========================================================================
%%                             CW fluence (emission)
%==========================================================================

if strcmpi(fwd,'fluo')
    phi_exc = phi;
    disp([spacer,'Fluo multiplication',spacer]); tic;
    sourceF = zeros(N, nQ);
    for i = 1:nQ
        sourceF(:,i) = hMesh.IntFG(phi_exc(:,i), conc);
    end
    phi = Smat\sourceF;
    
%     if verbosity == 1
%         figure;
%         for i=1:nQ
%             subplot(ceil(sqrt(nQ/1.5)),ceil(1.5*sqrt(nQ/1.5)),i);
%             toastShowMesh(hMesh,log(abs(phi(:,i))));
%             %toastShowMesh(DOT.mesh.hMesh,toastMapSolToMesh(DOT.grid.hBasis,log(abs(DOT.phi))));
%         end
%         suptitle('Emission fluence log(abs)');
%         drawnow;
%     end
end

tmp = mvec.' * phi;
proj = tmp(dmask(:));
%proj = proj';

if verbosity==1
    figure;
    imagesc(tmp.*dmask);title('CW emission');colorbar;
    ylabel('Meas');xlabel('Sources');
end


if strcmpi(fwd,'fluo')
    tmp = mvec' * phi_exc;
    proj_exc = tmp(dmask(:));
    %proj_exc = proj_exc';
    if verbosity == 1
    figure;
    imagesc(tmp.*dmask);title('Excitation Sinogram');axis square;colorbar
    ylabel('Meas');xlabel('Sources');
    % figure;
%     plot(proj_exc./proj);title('Born Sinogram');axis square;%colorbar
    end
end



% phi=reshape(phi,S.mesh.N,S.CameraS.Ns,S.rotation.Ntheta);
% for k = 1:S.rotation.Ntheta
%      proj(:,:,k) = S.M(:,:,k)'*phi(:,:,k);
% end

% if S.datalog ==1
%     proj=reshape(log(proj),[],1);
% else
%     proj=reshape(proj,[],1);
end