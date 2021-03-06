%-------------------------------------------------------------------------%
% Add a spherical inhomogeneity
%
% c         -- [1x3]    -- position of the center of the sphere
% var       -- string   -- the field you wanna mofify
% back      -- [1x1]    -- background value of 'var'
% DISTRIB   -- string   -- indicates wheter you wanna distribute the
%                          intensity within the inhomogeneity ('ON' or
%                          'OFF')
% INTENSITY -- [1x1]    -- maximum value of the inhomogeneous 'var'
% SHAPE     -- string   -- 'gaussian' or 'step'
% SIGMA     -- [1x1]    -- Width of the inhomogeneity
%
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 29/11/10
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 16/12/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 22/12/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 14/01/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 02/02/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 17/05/09
%-------------------------------------------------------------------------%


function [DOT,hete] = sphere3D(DOT, hete, solver)
%%--
SPE = 0;
SPE_CONC = 0;
SPE_SCA = 0;
sp = 1;
if contains(solver,'spectral')
    sp = 2;
end
c    = hete.c;
for is = 1:sp
    for itype = 1:size(hete.type,2)
        var  = hete.type{itype}; %change the inclusion type Mua in absorption and Musp in scattering
        if is >= 2
            SPE = 1;
            if strcmpi(var,'Mua')
                SPE_CONC = 1;
                SPE_SCA = 0;
            elseif strcmpi(var,'Musp')
                SPE_SCA = 1;
                SPE_CONC = 0;
            end
        end
        back = eval(['DOT.opt.',lower(hete.type{itype}),'B']); %muaB (and muspB in the next cycle) for background calculated from concentrations in TraslateXVar
        sigma= hete.sigma;
        intensity = hete.val((1:DOT.radiometry.nL)+(itype-1)*DOT.radiometry.nL); %val has the 8 values of absorption and 8 of scattering for the inclusion calculated in TranslateXVar
        shape     = hete.profile;
        distrib   = hete.distrib;
        if SPE
            if SPE_CONC
                back = eval(['DOT.opt.concB']);
                if iscolumn(back),back=back';end
                intensity = hete.conc;
                if iscolumn(intensity),intensity=intensity';end
            end
            if SPE_SCA
                back = [eval('DOT.opt.aB') eval('DOT.opt.bB')];
                intensity = [hete.a hete.b];
            end
        end
        %% ----------------------- distances au carree mesh centre ---------------%
        %-- Unstructured mesh --%
        if isfield(DOT,'grid')
            [X,Y,Z] = ndgrid(DOT.grid.x, DOT.grid.y, DOT.grid.z);
            X = reshape(X,DOT.grid.N,[]);
            Y = reshape(Y,DOT.grid.N,[]);
            Z = reshape(Z,DOT.grid.N,[]);
            M = [X-c(1) Y-c(2) Z-c(3)];
            add = zeros(DOT.grid.N,numel(intensity));
            Nx = DOT.grid.Nx;
            Ny = DOT.grid.Ny;
            Nz = DOT.grid.Nz;
            
            %-- Regular mesh --%
        else
            M = [DOT.mesh.pos(:,1)-c(1) DOT.mesh.pos(:,2)-c(2) DOT.mesh.pos(:,3)-c(3)];
            add = zeros(DOT.mesh.N,numel(intensity));
            Nx = DOT.mesh.Nx;
            Ny = DOT.mesh.Ny;
            Nz = DOT.mesh.Nz;
        end
        
        dist2 = sum(M.^2,2);
        if SPE
            if SPE_CONC
                var = 'Conc';
            end
            if SPE_SCA
                var = 'AB';
                DOT.opt.AB = cat(4,DOT.opt.A,DOT.opt.B);
            end
        end
        switch upper(shape)
            %% ---------------------- profil de concentration gaussien ---------------%
            case 'GAUSSIAN'
                %-- selection des indices --%
                indx = find(dist2 < 9*sigma*sigma); %in absolute value all the indexes of dist2 numbers different from zero less than 9*sigma*sigma
                %-- update concentration --%
                param = getfield(DOT.opt, var);
                switch upper(distrib)
                    case 'ON'
                        add(indx,:) = add(indx,:) + exp(-dist2(indx,1)/sigma/sigma);
                        add = add.*intensity./sum(add,1);
                    case 'OFF'
                        add(indx,:) = add(indx,:) + (intensity-back).*exp(-dist2(indx,1)/sigma/sigma);
                        %here the inclusion is created: for every value of absorption and
                        %scattering i subtract the background values
                end
                
                %% --------------------- profil de concentration creneau -----------------%
            case 'STEP'
                %-- selection des indices --%
                indx = find(dist2 <= sigma*sigma);
                %-- update concentration --%
                param = getfield(DOT.opt, var);
                switch upper(distrib)
                    case 'ON',      add = intensity./length(indx);
                    case 'OFF',
                        param = reshape(param,DOT.grid.N,numel(intensity));
                        param(indx,:) = param(indx,:) + (intensity-back);
                        param = reshape(param,Nx,Ny,Nz,numel(intensity));
                end
        end
        
        %% ---------------------------- Update -----------------------------------%
        add = squeeze(reshape(add,Nx,Ny,Nz,numel(intensity)));
        indx_dummy = indx;
        for inl = 1:numel(intensity)-1
            indx = [indx;indx_dummy+inl*DOT.grid.N];
        end
        param(indx) = param(indx) + add(indx);
        
        if strcmpi(var,'AB')
            hete    = setfield(hete, ['d','A'], add(:,:,:,1));
            hete    = setfield(hete, ['d','B'], add(:,:,:,2));
            DOT.opt = rmfield(DOT.opt,'AB');
            DOT.opt = setfield(DOT.opt, var(1), param(:,:,:,1));
            DOT.opt = setfield(DOT.opt, var(2), param(:,:,:,2));
        else
            hete    = setfield(hete, ['d',var], add);
            DOT.opt = setfield(DOT.opt, var, param);
        end
    end
end

