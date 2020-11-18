%-------------------------------------------------------------------------%
% Add a cylindrical inhomogeneity
%
% c         -- [1x3]    -- a point on the axis of the cylinder
% d         -- [1x3]    -- the direction of the axis of the cylinder
% var       -- string   -- the field you wanna mofify
% back      -- [1x1]    -- background value of 'var'
% INTENSITY -- [1x1]    -- maximum value of the inhomogeneous 'var'
% SHAPE     -- string   -- 'gaussian' or 'step'
% SIGMA     -- [1x1]    -- radius of the cylinder
%
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 08/02/10
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 06/07/10
% A. Farina - CNR-IFN - Dip. Fisica  - Politecnico di Milano - 10/04/15
%-------------------------------------------------------------------------%


function [DOT,hete] = cylinder3D(DOT,hete,solver,varargin)
%--
SPE = 0;
SPE_CONC = 0;
SPE_SCA = 0;
sp = 1;

if contains(solver,'spectral')
    sp = 2;
end

d = hete.d;     % direction vector
d = d./norm(d); % unitary direction vector
c = hete.c;    % first point on the axis
l = hete.l;    % length of the cylinder
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
        back = eval(['DOT.opt.',lower(hete.type{itype}),'B']);
        sigma = hete.sigma;
        intensity = hete.val((1:DOT.radiometry.nL)+(itype-1)*DOT.radiometry.nL);
        shape = hete.profile;
        distrib = hete.distrib;
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

%-- --%
        if isfield(DOT,'grid'),
            %-- distance to the axis of the cylinder --%
            D = [d(1)*ones(DOT.grid.N,1),d(2)*ones(DOT.grid.N,1),d(3)*ones(DOT.grid.N,1)];
            [X,Y,Z] = ndgrid(DOT.grid.x,DOT.grid.y,DOT.grid.z);
            X = reshape(X,DOT.grid.N,[]);
            Y = reshape(Y,DOT.grid.N,[]);
            Z = reshape(Z,DOT.grid.N,[]);
            M = [X-c(1) Y-c(2) Z-c(3)];
            dist = sum(cross(D,M).^2,2).^0.5./norm(d);
            add = zeros(DOT.grid.N,numel(intensity));
            Nx = DOT.grid.Nx;
            Ny = DOT.grid.Ny;
            Nz = DOT.grid.Nz;
            else
            %-- distance to the axis of the cylinder --%
            D = [d(1)*ones(DOT.mesh.N,1),d(2)*ones(DOT.mesh.N,1),d(3)*ones(DOT.mesh.N,1)];
            M = [DOT.mesh.pos(:,1)-c(1) DOT.mesh.pos(:,2)-c(2) DOT.mesh.pos(:,3)-c(3)];
            dist = sum(cross(D,M).^2,2).^0.5./norm(d);
            add = zeros(DOT.mesh.N,numel(intensity));
            Nx = DOT.grid.Nx;
            Ny = DOT.grid.Ny;
            Nz = DOT.grid.Nz;
            end
            %-- axial distance --%
            L = M*d';
            ind1 = find(L>=0);
            ind2 = find(L<=l);
            ind = intersect(ind1,ind2);
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
    
        %-- profil de concentration gaussien --%       
            case 'GAUSSIAN'
    
            %-- selection des indices --%
            indx = find(dist < 3*sigma);
            indx = intersect(indx,ind);

            %-- update concentration --%
           if isfield(DOT,'grid'),
            param = getfield(DOT.opt, var);
             else
                 a=lower(hete.type{itype});
                param=getfield(DOT.opt,a);
            end
    
            switch upper(distrib)
                case 'ON'
                    add(indx,:) = add(indx,:) + exp(-dist(indx,1).^2/sigma/sigma);
                    add = add.*intensity./sum(add,1);
                    add = squeeze(reshape(add,Nx,Ny,Nz,numel(intensity)));
                    indx_dummy = indx;
                    for inl = 1:numel(intensity)-1
                        indx = [indx;indx_dummy+inl*DOT.grid.N];
                    end
                    param(indx) = param(indx) + add(indx);    
    
                case 'OFF'  
                    add(indx,:) = add(indx,:) + (intensity-back).*exp(-dist(indx,1).^2/sigma/sigma);
                    add = squeeze(reshape(add,Nx,Ny,Nz,numel(intensity)));
                    indx_dummy = indx;
                    for inl = 1:numel(intensity)-1
                        indx = [indx;indx_dummy+inl*DOT.grid.N];
                    end
                    param(indx) = param(indx) + add(indx);

            end
    
            if strcmpi(var,'AB')
                hete    = setfield(hete, ['d','A'], add(:,:,:,1));
                hete    = setfield(hete, ['d','B'], add(:,:,:,2));
                DOT.opt = rmfield(DOT.opt,'AB');
                DOT.opt = setfield(DOT.opt, var(1), param(:,:,:,1));
                DOT.opt = setfield(DOT.opt, var(2), param(:,:,:,2));
            else
                hete    = setfield(hete, ['d',var], add);
                if isfield(DOT,'grid'),
                   DOT.opt = setfield(DOT.opt, var, param);
                else
                   DOT.opt = setfield(DOT.opt, lower(var), param);
                end
            end
    
%-- profil de concentration creneau --%       
            case 'STEP'
    
                %add = zeros(MOL.mesh.N,1);
    
                %-- selection des indices --%
                indx = find(dist <= sigma);
                indx = intersect(indx,ind);
    
                %-- update concentration --%
                if isfield(DOT,'grid'),
                    param = getfield(DOT.opt, var);
                else
                    a=lower(hete.type{itype});
                    param=getfield(DOT.opt,a);
                end
                
                switch upper(distrib)
                    case 'ON'
                        add(indx,:) = repmat(intensity./length(indx),numel(indx),1);
                        add = squeeze(reshape(add,Nx,Ny,Nz,numel(intensity)));
                        indx_dummy = indx;
                        for inl = 1:numel(intensity)-1
                            indx = [indx;indx_dummy+inl*DOT.grid.N];
                        end
                        param(indx) = param(indx) + add(indx);

                    case 'OFF'
                        add(indx,:) = repmat((intensity-back),numel(indx),1);
                        add = squeeze(reshape(add,Nx,Ny,Nz,numel(intensity)));
                        indx_dummy = indx;
                        for inl = 1:numel(intensity)-1
                            indx = [indx;indx_dummy+inl*DOT.grid.N];
                        end
                        param(indx) = param(indx) + add(indx);
                end
                
                hete    = setfield(hete, ['d',hete.type{itype}], add);
                if strcmpi(var,'AB')
                    hete    = setfield(hete, ['d','A'], add(:,:,:,1));
                    hete    = setfield(hete, ['d','B'], add(:,:,:,2));
                    DOT.opt = rmfield(DOT.opt,'AB');
                    DOT.opt = setfield(DOT.opt, var(1), param(:,:,:,1));
                    DOT.opt = setfield(DOT.opt, var(2), param(:,:,:,2));
                else
                    hete    = setfield(hete, ['d',var], add);
                    if isfield(DOT,'grid'),
                       DOT.opt = setfield(DOT.opt, var, param);
                    else
                       DOT.opt = setfield(DOT.opt, lower(var), param);
                    end
                end
        end
    end
end