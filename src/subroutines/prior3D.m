%-------------------------------------------------------------------------%
% Add a generic inhomogeneity
%
% var       -- string   -- the field you wanna mofify
% back      -- [1x1]    -- background value of 'var'
% INTENSITY -- [1x1]    -- maximum value of the inhomogeneous 'var'
% 
% No profiles are implemented. We get a binary tridimensional mask.
%-------------------------------------------------------------------------%





function [DOT, hete] = prior3D(DOT,hete,solver)

SPE = 0;
SPE_CONC = 0;
SPE_SCA = 0;
sp = 1;
if contains(solver,'spectral')
    sp = 2;
end

for is = 1:sp
    for itype = 1:size(hete.type,2)
        var  = hete.type{itype};
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
        intensity = hete.val((1:DOT.radiometry.nL)+(itype-1)*DOT.radiometry.nL);
        %distrib   = hete.distrib;
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


        if isfield(DOT,'grid')
            nx = DOT.grid.Nx;
            ny = DOT.grid.Ny;
            nz = DOT.grid.Nz;
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

        %% here!
        smask = load(hete.path);
        fn = fieldnames(smask);
        mask = smask.(fn{1});
        delta = smask.(fn{2});
        %% swap fields... in case
        if ~isvector(delta)
            dd = mask;
            mask = delta;
            delta = dd;
        end

        prior = priormask3D(hete.path,DOT.grid);
        param = getfield(DOT.opt, var);
        for inl = 1:numel(intensity)
            param(:,:,:,inl) = param(:,:,:,inl) + double(prior).*(intensity(inl)-back(inl));
        end
        
        DOT.opt = setfield(DOT.opt, var, param);
        
    end
end
return;
end



