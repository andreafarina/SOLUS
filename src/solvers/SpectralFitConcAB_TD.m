
%==========================================================================
% This function contains a solver for fitting optical properties of
% homogeneous phantom using routines in Matlab Optimization Toolbox
%
% Andrea Farina 10/15
%==========================================================================

function [bMua, bMusp, bConc,bB, bA] = SpectralFitConcAB_TD(solver,grid, conc0,a0,b0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd,TYPE_FWD,radiometry,spe,geometry)
BULK = 1;% set to 1 if wanting to fit the bulk.
INCL = 1;% set to 1 if wanting to fit the inclusion defined by refimage
        % set to INCL = 0 and BULK = 1 to make a homogenous fit
geom = geometry;
refimage = solver.prior.refimage;
self_norm = true;
if strcmpi(TYPE_FWD, 'linear')
    warning('TYPE_FWD set to linear but it should be fem instead. TYPE_FWD set to fem')
    TYPE_FWD = 'fem';    
end
weight_type = 'none';%'none'
min_func  = 'lsqrfit';%'lsqrcurvefit'
first_lim = 0.1; last_lim = 0.1;
ForceConstitSolution = spe.ForceConstitSolution;

nQM = sum(dmask(:));
nwin = size(twin,1);

%% Inverse solver
mua0 = (spe.ext_coeff0*conc0)';          
mus0 = (a0 .* (spe.lambda./spe.lambda0).^(-b0));

[proj, Aproj] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom,TYPE_FWD,radiometry);
if numel(irf)>1
    for inl = 1:radiometry.nL
        meas_set = (1:nQM)+(inl-1)*nQM;
        proj_single = proj(:,meas_set);
        z = convn(proj_single,irf(:,inl));
        nmax = max(nstep,numel(irf(:,inl)));
        proj_single = z(1:nmax,:);
        clear nmax
        
        if self_norm == true
            proj_single = proj_single * spdiags(1./sum(proj_single,'omitnan')',0,nQM,nQM);
        end
        clear z
        if inl == 1
            proj(1:size(proj_single),meas_set) = proj_single;
            proj(size(proj_single,1)+1:end,:) = [];
        else
            proj(1:size(proj_single,1),meas_set) = proj_single;
        end
    end
end
if self_norm == true
    for inl = 1:radiometry.nL
        meas_set = (1:nQM)+(inl-1)*nQM;
        data_single = data(:,meas_set);
        data_single = data_single * spdiags(1./sum(data_single,'omitnan')',0,nQM,nQM);
        ref_single = ref(:,meas_set);
        ref_single = ref_single * spdiags(1./sum(ref_single,'omitnan')',0,nQM,nQM);
        data(:,meas_set) = data_single;
        ref(:,meas_set) = ref_single;
    end
end
dummy_proj = zeros(size(twin,1),nQM*radiometry.nL);
for inl = 1:radiometry.nL
    meas_set =(1:nQM)+(inl-1)*nQM; twin_set = (1:2)+(inl-1)*2;
    proj_single = proj(:,meas_set);
    proj_single = WindowTPSF(proj_single,twin(:,twin_set));
    if self_norm == true
        proj_single = proj_single * spdiags(1./sum(proj_single,'omitnan')',0,nQM,nQM);
    end
    dummy_proj(:,meas_set) = proj_single;
end
weight = ones(nwin,nQM*radiometry.nL);
if ~strcmpi(weight_type,'none')
    weight = zeros(nwin,nQM*radiometry.nL);
    interval = zeros(2,nQM*radiometry.nL);
    for im = 1:nQM*radiometry.nL
        idx1 = find(ref(:,im)>max(ref(:,im)*first_lim),1,'first');
        idx2 = find(ref(:,im)>max(ref(:,im)*last_lim),1,'last');
        weight(idx1:idx2,im) = 1;
        interval(1,im) = idx1+(im-1)*nwin;
        interval(2,im) = idx2+(im-1)*nwin;
    end
    weight = weight(:);
    figure, semilogy(ref(:)), vline(interval(:),'r');
    ref =ref(:).*weight;
end
proj = dummy_proj;
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
%weight(mask) = [];
figure(1002);semilogy([proj,data]),legend('proj','ref')
sd = sqrt(ref);%%ones(size(proj));%proj(:);
%sd = ones(size(data));
data = data./sd;


%% solution vector
% x = [mua0;mus0;-0];     % [mua0,mus0,t0]


%% Setup optimization for lsqcurvefit
if spe.SPECTRA == 0 && ForceConstitSolution == false
    FinDiffRelStep = 10*repmat(1e-3,spe.nLambda*2,INCL+BULK);
else
    FinDiffRelStep = 10*repmat(1e-3,spe.nCromo+2,INCL+BULK);
end

%% Setup optimization for lsqnonlin
% opts = optimoptions('lsqnonlin',...
%     'Jacobian','off',...
%     ...'Algorithm','levenberg-marquardt',...
%     'DerivativeCheck','off',...
%     'MaxIter',20,'Display','iter-detailed')
%% Setup optimization for fminunc
% opts = optimoptions('fminunc','GradObj','on','Algorithm','quasi-newton','MaxIter',2,...
%   'Display','iter')


%% initialise solution vector
if (BULK == 1) && (INCL == 1)  
    x0 = [[conc0', a0, b0];[conc0', a0, b0]];   
elseif ( (BULK == 0) && (INCL == 1) )|| ((BULK == 1) && (INCL == 0))    
    x0 = [conc0', a0, b0];
    xB = x0;
else
    disp('No region to fit has been selected')
    return;
end
    
%% Solve

%x = fminsearch(@objective,x0);
if spe.SPECTRA == 0 && ForceConstitSolution == false
 %   x0 = {mua0 mus0};
    low_bound = zeros(INCL + BULK,numel([x0{:}]));
else
    if ForceConstitSolution
        spe.opt.conc0 = ones(spe.nCromo,1).*spe.active_cromo';
        spe.opt.a0 = 1; spe.opt.b0 = 1;
    end
%    x0 = {spe.opt.conc0',[spe.opt.a0 spe.opt.b0]};
    low_bound = zeros(INCL + BULK,spe.nCromo+2);
end

 
if strcmpi(min_func,'lsqrfit')
    opts = optimoptions('lsqcurvefit',...
    'Jacobian','off',...
    'Algorithm','trust-region-reflective',...
    'DerivativeCheck','off',...
    'MaxIter',200*radiometry.nL,'Display','iter-detailed',...
    'FinDiffRelStep', FinDiffRelStep,'TolFun',1e-10,'TolX',1e-10);
   
     [x,res] = lsqcurvefit(@comp_forward,x0,[],data,low_bound,[],opts);
elseif strcmpi(min_func,'fminunc')
    opts = optimoptions('fminunc',...
    'Algorithm','quasi-newton',...
    'Display','iter-detailed',...
     'GradObj','off',...
    'MaxIter',200 * radiometry.nL);    
    [x, res] = fminunc(@Loss_func,x0, opts);
end
%x0 = [t0 x0{:}];
%[x,res] = lsqcurvefit(@forward,x0,[],data,low_bound,[],opts);

%x = lsqnonlin(@objective2,x0,[],[],opts);

%% display fit result
disp(['Residual: ' num2str(res)])
% chromofores, b & a
if ((BULK == 1) && (INCL == 1))
    bConc = x(1,1:end-2) .* ~refimage(:) + x(2,1:end-2) .* refimage(:);   
    bA = x(1,end-1) * ~refimage(:) + x(2,end-1) * refimage(:);   
    bB = x(1,end) * ~refimage(:) + x(2,end) * refimage(:); 

elseif ((BULK == 1) && (INCL == 0))
    bConc = x(1,1:end-2)' .* ones(size(refimage(:)));   
    bA = x(1,end-1) * ones(size(refimage(:))) ;   
    bB = x(1,end) * ones(size(refimage(:))); 
elseif (BULK == 0) && (INCL == 1 )
    bConc = xB(1,1:end-2) .* ~refimage(:) + x(2,1:end-2) .* refimage(:);   
    bA = xB(1,end-1) * ~refimage(:) + x(2,end-1) * refimage(:);   
    bB = xB(1,end) * ~refimage(:) + x(2,end) * refimage(:); 
end

% Optical Coefficients from chromofores, b & a
bMua = (spe.ext_coeff0 * bConc(:,:)')' ;
bMusp = bA.*(spe.lambda./spe.lambda0).^(-bB);



return

%% extract the amplitude area
self_norm = 0;
mask = true(nwin*nQM,1);
[proj_fit, Aproj_fit] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, bmua, bmus, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom, TYPE_FWD,radiometry);
% proj_fit = circshift(proj_fit,round(x(3)/dt));
if numel(irf)>1
    z = convn(proj_fit,irf);
    nmax = max(nstep,numel(irf));
    proj_fit = z(1:nmax,:);
    clear nmax
end
proj_fit = WindowTPSF(proj_fit,twin);
Aproj_fit = sum(proj_fit);
A_data = sum(data2);
factor = Aproj_fit./A_data;
save('factor_ref.mat','factor');

%% Callback function OBJECTIVE for forward problem
    function [proj,J] =comp_forward(x,~)
        %xx = [x(1)*ones(nsol,1);x(2)*ones(nsol,1)];
        if spe.SPECTRA == 0 && ForceConstitSolution == false
            mua = x(2:spe.nLambda+1);
            mus = x(spe.nLambda+2:end);
        else
            if (BULK == 1) && (INCL == 1)
                b_out = x(1,end); a_out = x(1,end-1); conc_out = x(1,1:end-2).*spe.active_cromo;
                b_in = x(2,end); a_in = x(2,end-1); conc_in = x(2,1:end-2).*spe.active_cromo;
                mua_out = (spe.ext_coeff0*conc_out');
                mus_out = a_out.*(spe.lambda./spe.lambda0).^(-b_out);
                mua_in = (spe.ext_coeff0*conc_in');
                mus_in = (a_in.*(spe.lambda./spe.lambda0).^(-b_in));
                mua = zeros([size(refimage), radiometry.nL]);
                mus = zeros([size(refimage), radiometry.nL]);
                for inL = 1: radiometry.nL
                    mua(:,:,:, inL) = mua_out(inL) * (1 -refimage) + mua_in(inL) * refimage;
                    mus(:,:,:, inL) = mus_out(inL) * (1 -refimage) + mus_in(inL) * refimage;
                end
            elseif (BULK == 1) && (INCL == 0)
                b_out = x(1,end); a_out = x(1,end-1); conc_out = x(1,1:end-2).*spe.active_cromo;
                mua_out = (spe.ext_coeff0*conc_out');
                mus_out = a_out.*(spe.lambda./spe.lambda0).^(-b_out);
                mua = zeros([size(refimage), radiometry.nL]);
                mus = zeros([size(refimage), radiometry.nL]);
                for inL = 1: radiometry.nL
                    mua(:,:,:, inL) = mua_out(inL) * ones(size(refimage));
                    mus(:,:,:, inL) = mus_out(inL) * ones(size(refimage));
                end
            elseif (BULK == 0) && (INCL == 1 )
                b_out = xB(1,end); a_out = xB(1,end-1); conc_out = xB(1,1:end-2).*spe.active_cromo;
                b_in = x(1,end); a_in = x(1,end-1); conc_in = x(1,1:end-2).*spe.active_cromo;
                mua_out = (spe.ext_coeff0*conc_out');
                mus_out = a_out.*(spe.lambda./spe.lambda0).^(-b_out);
                mua_in = (spe.ext_coeff0*conc_in');
                mus_in = (a_in.*(spe.lambda./spe.lambda0).^(-b_in));
                mua = zeros([size(refimage), radiometry.nL]);
                mus = zeros([size(refimage), radiometry.nL]);
                for inL = 1: radiometry.nL
                    mua(:,:,:, inL) = mua_out(inL) * ~(refimage) + mua_in(inL) * refimage;
                    mus(:,:,:, inL) = mus_out(inL) * ~(refimage) + mus_in(inL) * refimage;
                end 
            end

        end
        [proj, Aproj] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, mua_out, mus_out, n, ...
            mua,mus, A, dt, nstep, self_norm,...
            geom,TYPE_FWD,radiometry);
        if spe.SPECTRA == 0 && ForceConstitSolution == false
            display(['mua = ',num2str(mua)]);
            display(['musp = ',num2str(mus)]);
        else
            if BULK == 1
            disp('bulk:');disp(spe.cromo_label); disp(conc_out); disp({'a','b'}); disp([a_out b_out]);
            end
            if INCL == 1
            disp('inclusion:');disp(spe.cromo_label); disp(conc_in); disp({'a','b'}); disp([a_in b_in]);
            end
        end
        % [~,proj] = Contini1997(0,(1:nstep)*dt/1000,20,mua(1),mus(1),1,n(1),'slab','Dmus',200);
        % proj = proj';%./sum(proj);
        
        if numel(irf)>1
            for inl_ = 1:radiometry.nL
                meas_set_ = (1:nQM)+(inl_-1)*nQM;
                proj_single_ = proj(:,meas_set_);
                z = convn(proj_single_,irf(:,inl_));
                nmax = max(nstep,numel(irf(:,inl_)));
                proj_single_ = z(1:nmax,:);
                clear nmax
                if self_norm == true
                    proj_single_ = proj_single_ * spdiags(1./sum(proj_single_,'omitnan')',0,nQM,nQM);
                end
                clear z
                if inl_ == 1
                    proj(1:size(proj_single_,1),meas_set_) = proj_single_;
                    proj(size(proj_single_,1)+1:end,:) = [];
                else
                    proj(1:size(proj_single_,1),meas_set_) = proj_single_;
                end
            end
        end
        clear meas_set_
        dummy_proj_ = zeros(size(twin,1),sum(dmask(:))*radiometry.nL);
        for inl_ = 1:radiometry.nL
            meas_set_ = (1:nQM)+(inl_-1)*nQM;twin_set_ = (1:2)+(inl_-1)*2;
            proj_single_ = proj(:,meas_set_);
            proj_single_ = WindowTPSF(proj_single_,twin(:,twin_set_));
            if self_norm == true
                proj_single_ = proj_single_ * spdiags(1./sum(proj_single_,'omitnan')',0,nQM,nQM);
            end
            dummy_proj_(:,meas_set_) = proj_single_;
        end
        proj = dummy_proj_(:);
        if ~strcmpi(weight_type,'none')
            proj = proj.*weight;
        end
        proj(mask) = [];
        proj = proj(:)./sd;
        % plot forward
        t = (1:numel(data)) * dt;
        figure(1003);
        
        semilogy(t,proj,'-',t,data,'.'),ylim([1e-3 1])
        %         semilogy([proj data]),ylim([1e-3 1])
        %         vline(sum(dmask(:))*size(twin,1)*(1:radiometry.nL),repmat({'r'},radiometry.nL,1))
        drawnow;
        nwin = size(twin,1);

    end
 

    function L = Loss_func(x,~)
        out = comp_forward(x);
        L = sum((out(:) - data(:)).^2);
    end


end