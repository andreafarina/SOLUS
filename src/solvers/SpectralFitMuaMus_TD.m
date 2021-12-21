%==========================================================================
% This function contains a solver for fitting optical properties of
% homogeneous phantom using routines in Matlab Optimization Toolbox
%
% Andrea Farina 10/15
%==========================================================================

function [bMua,bMus,out_conc,out_a,out_b] = SpectralFitMuaMus_TD(~,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd,~,radiometry,spe)
geom = 'semi-inf';
self_norm = 1;
weight_type = 'none'; %'rect';
first_lim = 0.0; last_lim = 0.0;
ForceConstitSolution = spe.ForceConstitSolution;
PLOT = 0;
% mua0 = 0.01;
% mus0 = 1.0;

%spe.ext_coeff0(:,5) = 0.001;
%spe.ext_coeffB = spe.ext_coeff0;

nQM = sum(dmask(:));
nwin = size(twin,1);
%% Inverse solver
[proj, Aproj] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom,'linear',radiometry,irf);

if self_norm == true
    tmpA = sum(data,1);
    data = NormalizeTPSF(data);
    ref = NormalizeTPSF(ref);
    sd = 1+0*sqrt(data).*sqrt(1./tmpA);
end

dummy_proj = zeros(size(twin,1),nQM*radiometry.nL);

proj_single = WindowTPSF(proj,twin);
if self_norm == true
    proj_single = NormalizeTPSF(proj_single);
end

dummy_proj = proj_single;
weight = ones(nwin,nQM);
if ~strcmpi(weight_type,'none')
    weight = zeros(nwin,nQM);
    interval = zeros(2,nQM);
    for im = 1:nQM
        idx1 = find(ref(:,im)>max(ref(:,im)*first_lim),1,'first');
        idx2 = find(ref(:,im)>max(ref(:,im)*last_lim),1,'last');
        weight(idx1:idx2,im) = 1;
        interval(1,im) = idx1+(im-1)*nwin;
        interval(2,im) = idx2+(im-1)*nwin;
    end
    weight = weight(:);
    if PLOT
        figure, semilogy(ref(:)), vline(interval(:),'r');
    end
    ref =ref(:).*weight;
end
proj = dummy_proj;
proj = proj(:);
data = data(:);
ref = ref(:);
sd = sd(:);
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
sd(mask) = [];
proj(mask) = [];
%weight(mask) = [];
if PLOT
    figure(1002);semilogy([proj,data]),legend('proj','ref')
end

if self_norm
end
%sd = sqrt(ref);%%ones(size(proj));%proj(:);
%sd = ones(size(data));
data = data./sd;


%% solution vector
% x = [mua0;mus0;-0];     % [mua0,mus0,t0]


%% Setup optimization for lsqcurvefit
if spe.SPECTRA == 0 && ForceConstitSolution == false
    FinDiffRelStep = [2; repmat(1e-3,spe.nLambda*2,1)];
else
    FinDiffRelStep = [1,1e-1,1e-1,1e3,1e3,1e3,1e-1,1e-1 ]*0.00005;%repmat(1e-1,spe.nCromo+2,1)];
    TypX = [1e-2,1e0,1e0,1e3,1e3,1e2,1e0,1e0];
end
opts = optimoptions('lsqcurvefit',...
    'Jacobian','off',...
    'Algorithm','levenberg-marquardt','InitDamping',1e-10,...%'TypicalX',TypX,...
    'DerivativeCheck','off','FiniteDifferenceType','forward',...
    'MaxIter',3200,'Display','iter-detailed','FinDiffRelStep',FinDiffRelStep,'MaxFunctionEvaluation',15000,...
    'TolFun',1e-14,'TolX',1e-14)
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
if spe.SPECTRA == 0 && ForceConstitSolution == false
    x0 = {mua0 mus0};
    low_bound = [-20 zeros(1,numel([x0{:}]))];
else
    if ForceConstitSolution
        spe.opt.conc0 = [1,1,1000,1000,200]'.*ones(spe.nCromo,1).*spe.active_cromo';
        spe.opt.a0 = 1; spe.opt.b0 = 1;
    end
    x0 = {spe.opt.conc0,[spe.opt.a0 spe.opt.b0]};
    low_bound = [-100 zeros(1,spe.nCromo+2)];
end
t0 = 0;
if ~isrow(x0{1}), x0{1} = x0{1}';end
    
x0 = [t0 [x0{1} x0{2} ]];
[x,res] = lsqcurvefit(@forward,x0,[],data,low_bound,[],opts);

%x = lsqnonlin(@objective2,x0,[],[],opts);

%% display fit result
disp(['Residual: ' num2str(res)])
if spe.SPECTRA == 0 && ForceConstitSolution == false
    bmua = x(2:spe.nLambda+1);
    bmus = x(spe.nLambda+2:end);
    out_conc = 0;
    out_a = 0;
    out_b = 0;
else
    if isrow(x), x = x'; end
    bmua = spe.ext_coeffB*x(2:end-2);
    bmus = x(end-1).*(spe.lambda/spe.lambda0).^(-x(end));
    if iscolumn(bmua), bmua = bmua'; end
    if iscolumn(bmus), bmus = bmus'; end
    out_conc= (x(2:end-2));
    out_a=(x(end-1));
    out_b=(x(end));
    display(['a = '  num2str(x(end-1))]); display(['b =' num2str(x(end))]);
    display([char(spe.cromo_label) repmat('=',spe.nCromo,1) num2str(x(2:end-2))]);
    disp([char(spe.cromo_label) repmat('=',spe.nCromo,1) num2str(x(2:end-2)) repmat(' ',spe.nCromo,1) char(spe.cromo_units)]);
end
fprintf(['<strong>mua = ',num2str(bmua),'</strong>\n']);
fprintf(['<strong>musp = ',num2str(bmus),'</strong>\n']);
display(['t0 = ',num2str(x(1))]);
% display(['MuaErr= ',num2str(bmua-mua0)])
% display(['MusErr= ',num2str(bmus-mus0)])
% display(['MuaErr%= ',num2str(((bmua-mua0)./mua0).*100)])
% display(['MusErr%= ',num2str(((bmus-mus0)./mus0).*100)])

%% Map parameters back to mesh

bMua = bmua.*ones(grid.N,1);
bMus = bmus.*ones(grid.N,1);
out_conc = out_conc'.*ones(grid.N,numel(out_conc));
out_a = out_a.*ones(grid.N,1);
out_b = out_b.*ones(grid.N,1);
return

%% extract the amplitude area
self_norm = 0;
mask = true(nwin*nQM,1);
[proj_fit, Aproj_fit] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, bmua, bmus, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom, 'linear',radiometry);
% proj_fit = circshift(proj_fit,round(x(3)/dt));

proj_fit = ConvIRF(proj_fit, irf);
proj_fit = WindowTPSF(proj_fit,twin);
Aproj_fit = sum(proj_fit);
A_data = sum(data2);
factor = Aproj_fit./A_data;
self_norm = 1;
save('factor_ref.mat','factor');

%% Callback function for forward problem
    function [proj,J] = forward(x,~)
        %xx = [x(1)*ones(nsol,1);x(2)*ones(nsol,1)];
        disp(x)
        t0_ = x(1);
        if spe.SPECTRA == 0 && ForceConstitSolution == false
            mua = x(2:spe.nLambda+1);
            mus = x(spe.nLambda+2:end);
        else
            b_ = x(end); a_ = x(end-1); conc_ = x(2:end-2).*spe.active_cromo;
            if isrow(conc_), conc_ = conc_'; end
            mua = (spe.ext_coeff0*conc_)';
            mus = a_.*(spe.lambda./spe.lambda0).^(-b_);
        end
        [proj, Aproj] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, mua, mus, n, ...
            [],[], A, dt, nstep, self_norm,...
            geom,'linear',radiometry,irf);
%         if spe.SPECTRA == 0 && ForceConstitSolution == false
%             display(['mua = ',num2str(mua)]);
%             display(['musp = ',num2str(mus)]);
%         else
%             disp(spe.cromo_label); disp(conc_'); disp({'a','b'}); disp([a_ b_]);
%         end
        % [~,proj] = Contini1997(0,(1:nstep)*dt/1000,20,mua(1),mus(1),1,n(1),'slab','Dmus',200);
        % proj = proj';%./sum(proj);
        

        clear meas_set_
        proj_single_ = WindowTPSF(proj,twin);
        if self_norm == true
            proj_single_ = NormalizeTPSF(proj_single_);
        end
        dummy_proj_ = proj_single_;
        proj = dummy_proj_(:);
        if ~strcmpi(weight_type,'none')
            proj = proj.*weight;
        end
        proj(mask) = [];
        proj = proj(:)./sd(:);
        
        if PLOT
            % plot forward
            t = (1:numel(data)) * dt;
            figure(1003);
            semilogy(t,proj,'-',t,data,'.'),ylim([1e-3 1])
            %         semilogy([proj data]),ylim([1e-3 1])
            %         vline(sum(dmask(:))*size(twin,1)*(1:radiometry.nL),repmat({'r'},radiometry.nL,1))
            drawnow;
            
        end
        nwin = size(twin,1);
        if nargout>1
            JJ = Jacobian (mua, mus, qvec, mvec);
            %save('J1','JJ');
            njac = size(JJ,2)/2;
            % Normalized measruements
            if self_norm == true
                for i=1:nQM
                    sJ = sum(JJ((1:nwin)+(i-1)*nwin,:),'omitnan');
                    sJ = repmat(sJ,nwin,1);
                    sJ = spdiags(proj((1:nwin)+(i-1)*nwin),0,nwin,nwin) * sJ;
                    JJ((1:nwin)+(i-1)*nwin,:) = (JJ((1:nwin)+(i-1)*nwin,:) - sJ)./Aproj(i);
                end
            end
            J(:,1) = sum(JJ(:,1:njac),2);% * 0.3;
            J(:,2) = sum(JJ(:,njac + (1:njac)),2);% * 0.3;
            % J = spdiags(1./proj,0,nQM*nwin,nQM*nwin) * J;
        end
    end
end