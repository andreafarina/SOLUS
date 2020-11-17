%==========================================================================
% This function contains a solver for fitting optical properties of
% homogeneous phantom using routines in Matlab Optimization Toolbox
%
% Andrea Farina 10/15
% Andrea Farina 11/20 fixed error arising by fitting with fraction
%==========================================================================

function [bmua,bmus] = FitMuaMus_TD(~,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, type_fwd)
geom = 'semi-inf';
weight_type = 'rect';%'none'; % 'none','rect'
fract_first = 0.5; fract_last = 0.01;
data = data;%ref;%data;%ref;
data2 = data;
nQM = sum(dmask(:));
nwin = size(twin,1);
Jacobian = @(mua, mus) JacobianTD (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom, 'muaD', type_fwd, self_norm,0);
%% Inverse solver
[proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
                [],[], A, dt, nstep, self_norm,...
                geom,'linear',irf);
    proj = WindowTPSF(proj,twin);
    switch lower(weight_type)
        case 'none'
            if self_norm == true
                proj = NormalizeTPSF(proj);
            end
        case 'rect'
            [ROI,ROI_d] = ThresholdTPSF(data,fract_first,fract_last);
            data = NormalizeTPSF(data.*ROI_d);
            data = data(ROI);
            sd = sd(ROI);
            ref = ref(ROI);
            
            if self_norm == true
                proj = NormalizeTPSF(proj.*ROI_d);
            end
            proj = proj(ROI);
    end

proj = proj(:);
data = data(:);
ref = ref(:);
%factor = proj./ref;

%data = data .* factor;
%ref = ref .* factor;
%% data scaling
%sd = sd(:).*(factor);
%% mask for excluding zeros
mask = ((ref(:).*data(:)) == 0) | ...
    (isnan(ref(:))) | (isnan(data(:)));
%mask = false(size(mask));


if ref == 0
    ref = proj(:);
end

ref(mask) = []; 
data(mask) = [];
proj(mask) = [];
figure(1002);semilogy([proj,data]),legend('proj','ref')
sd = sqrt(ref);

data = data./sd;


%% solution vector
x = [mua0;mus0;-0];     % [mua0,mus0,t0]


%% Setup optimization for lsqcurvefit
opts = optimoptions('lsqcurvefit',...
    'Jacobian','off',...
    ...'Algorithm','levenberg-marquardt',...
    'DerivativeCheck','off',...
    'MaxIter',100,'Display','final-detailed','FinDiffRelStep',[1e-3,1e-2,2]);%,'TolFun',1e-10,'TolX',1e-10)
%% Setup optimization for lsqnonlin
% opts = optimoptions('lsqnonlin',...
%     'Jacobian','off',...
%     ...'Algorithm','levenberg-marquardt',...
%     'DerivativeCheck','off',...
%     'MaxIter',20,'Display','iter-detailed')
%% Setup optimization for fminunc
%opts = optimoptions('fminunc','GradObj','on','Algorithm','quasi-newton','MaxIter',2,...
%   'Display','iter')
% opts = optimoptions('fminunc',...
%     'Algorithm','quasi-newton',...
%     'Display','iter-detailed',...
%      'GradObj','off',...
%     'MaxIter',20)


%% Solve
%x = fminunc(@objective,x0,opts);
%x = fminsearch(@objective,x0);
x0 = x;
%if strcmpi(weight_type,'none')
    x = lsqcurvefit(@forward,x0,[],data,[],[],opts);
%else    
%    x = fminunc(@Loss_func,x0);
%end
%x = lsqnonlin(@objective2,x0,[],[],opts);

%% Map parameters back to mesh

bmua = x(1)*ones(grid.N,1);
bmus = x(2)*ones(grid.N,1);


%% display fit result
% AF: occhio che mua0 e mus0 sono i valori da cui parte il fit e non quelli
% simulati. Non ha senso calcolare l'errore su quei dati.
fprintf(['<strong>mua = ',num2str(bmua(1)),'</strong>\n']);
fprintf(['<strong>musp = ',num2str(bmus(1)),'</strong>\n']);
fprintf(['<strong>t0 = ',num2str(x(3)),'</strong>\n']);
% display(['MuaErr= ',num2str(bmua(1)-mua0)])
% display(['MusErr= ',num2str(bmus(1)-mus0)])
% display(['MuaErr%= ',num2str(((bmua(1)-mua0)./mua0).*100)])
% display(['MusErr%= ',num2str(((bmus(1)-mus0)./mus0).*100)])
%% extract the amplitude area
self_norm = 0;
mask = true(nwin*nQM,1);
[~, Aproj_fit] = ForwardTD(grid,Spos, Dpos, dmask, x(1), x(2), n, ...
                [],[], A, dt, nstep, self_norm,...
                geom, 'linear',irf);
%proj_fit = circshift(proj_fit,round(x(3)/dt));

%proj_fit = WindowTPSF(proj_fit,twin);
A_data = sum(data2);
factor = Aproj_fit./A_data;
save('factor_ref.mat','factor');
%% ===================== OBJECTIVE FUNCTIONS=============================

%% Callback function for objective evaluation
    function [p,g] = objective(x,~) %#ok<DEFNU>
        verbosity = 1;
        xx = [x(1)*ones(nsol,1);x(2)*ones(nsol,1)];
        [mua,mus] = toastDotXToMuaMus(hBasis,xx,refind);
        mua(notroi) = mua0;
        mus(notroi) = mus0;
        %mua(mua<0) = 1e-4;
        %mus(mus<0) = 0.2;
%         for j = 1:length(mua) % ensure positivity
%             mua(j) = max(1e-4,mua(j));
%             mus(j) = max(0.2,mus(j));
%         end
%         
        [proj,Aproj] = ProjectFieldTD(hMesh,qvec,mvec,dmask,...
            mua,mus,conc,tau,n,dt,nstep, 0,self_norm,'diff',0);
        if numel(irf)>1
            for i = 1:nQM
                z(:,i) = conv(full(proj(:,i)),irf);
            end
            nmax = max(nstep,numel(irf));
            proj = z(1:nmax,:);
            clear nmax
            if self_norm == true
                proj = proj * spdiags(1./sum(proj,'omitnan')',0,nQM,nQM);
            end
            clear z
        end
        
        proj = WindowTPSF(proj,twin);
        proj = proj(:);
        
        %proj = privProject (hMesh, hBasis, x, ref, freq, qvec, mvec);
        [p, p_data, p_prior] = privObjective (proj, data, sd);
        if verbosity > 0
            fprintf (1, '    [LH: %f, PR: %f]\n', p_data, p_prior);
        end
        nwin = size(twin,1);
        if nargout>1
            JJ = Jacobian (mua, mus, qvec, mvec);
            
            % Normalized measruements
            if self_norm == true
                for i=1:nQM
                    sJ = sum(JJ((1:nwin)+(i-1)*nwin,:),'omitnan');
                    sJ = repmat(sJ,nwin,1);
                    sJ = spdiags(proj((1:nwin)+(i-1)*nwin),0,nwin,nwin) * sJ;
                    JJ((1:nwin)+(i-1)*nwin,:) = (JJ((1:nwin)+(i-1)*nwin,:) - sJ)./Aproj(i);
                end
            end
            J(:,1) = sum(JJ(:,1:nsol),2);
            J(:,2) = sum(JJ(:,nsol + (1:nsol)),2);
            g = - 2 * J' * ((data-proj)./sd);
        end
    end

%% Callback function for objective evaluation
    function [p,J] = objective2(x,~) %#ok<DEFNU>
        xx = [x(1)*ones(nsol,1);x(2)*ones(nsol,1)];
        [mua,mus] = toastDotXToMuaMus(hBasis,xx,refind);
        mua(notroi) = mua0;
        mus(notroi) = mus0;
        %mua(mua<0) = 1e-4;
        %mus(mus<0) = 0.2;
%         for j = 1:length(mua) % ensure positivity
%             mua(j) = max(1e-4,mua(j));
%             mus(j) = max(0.2,mus(j));
%         end
%         
        [proj,Aproj] = ProjectFieldTD(hMesh,qvec,mvec,dmask,...
            mua,mus,conc,tau,n,dt,nstep, 0,self_norm,'diff',0);
        if numel(irf)>1
            for i = 1:nQM
                z(:,i) = conv(full(proj(:,i)),irf);
            end
            proj = z(1:nstep,:);
            if self_norm == true
                proj = proj * spdiags(1./sum(proj,'omitnan')',0,nQM,nQM);
            end
            clear z
        end
        proj = WindowTPSF(proj,twin);
        proj = proj(:);
        
        p = (data - proj)./sd;
        
        nwin = size(twin,1);
        if nargout>1
            JJ = Jacobian (mua, mus, qvec, mvec);
            JJ = spdiags(1./sd,0,nQM*nwin,nQM*nwin) * JJ;
            % Normalized measruements
            
            J(:,1) = - sum(JJ(:,1:nsol),2);
            J(:,2) = - sum(JJ(:,nsol + (1:nsol)),2);
       end
    end
%% Callback function for forward problem
function [proj,J] = forward(x,~)
    %xx = [x(1)*ones(nsol,1);x(2)*ones(nsol,1)];
    t0 = x(3);
    [proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, x(1), x(2), n, ...
                [],[], A, dt, nstep, self_norm,...
                geom, 'linear',irf);
   
    proj = circshift(proj,round(t0/dt));
    proj = WindowTPSF(proj,twin);
            
    switch lower(weight_type)
        case 'none'
            if self_norm == true
                proj = NormalizeTPSF(proj);
            end
        case 'rect'
            if self_norm == true
                proj = NormalizeTPSF(proj.*ROI_d);
            end
            proj = proj(ROI);
            
    end
        
    
    proj(mask) = [];
    proj = proj(:)./sd;
    
% plot forward
    t = (1:numel(data)) * dt;
    figure(1003);
    semilogy(t,proj,'-',t,data,'.'),ylim([1e-3 1])
    drawnow;
    nwin = size(twin,1);
    if nargout>1
        JJ = Jacobian (x(1), x(2));
        %save('J1','JJ');
        njac = size(JJ,2)/2;
        % Normalized measruements
        
        J(:,1) = sum(JJ(:,1:njac),2);% * 0.3;
        J(:,2) = sum(JJ(:,njac + (1:njac)),2);% * 0.3;
       % J = spdiags(1./proj,0,nQM*nwin,nQM*nwin) * J;
    end
end
function L =  Loss_func(coeff,~) %#ok<DEFNU>
    fwd = forward(coeff);
    err_square = (fwd-data).^2;
    data_=reshape(data,numel(fwd)/nQM,nQM);
    switch lower(weight_type)
        case 'rect'
            weight = zeros(numel(fwd)/nQM,nQM);
            for im = 1:nQM
                idx1 = find(data_(:,im)>max(data_(:,im))*fract_first,1,'first');
                idx2 = find(data_(:,im)>max(data_(:,im))*fract_last,1,'last');
                weight(idx1:idx2,im) = 1;
                interval(1,im) = idx1+(im-1)*numel(fwd)/nQM; interval(2,im)= idx2+(im-1)*numel(fwd)/nQM;
            end
            L=sum(weight(:).*err_square);
            vline(dt*interval(:),'r');
    end
end
end
