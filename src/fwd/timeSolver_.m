function tPhi = timeSolver_(mode, phi, dt, nstep, theta, param, USEGPU)

%theta = 1; %0.5;
tol = 1e-8;

persistent timeSolver

switch mode
  case 'compute'
    c = 0.299./param.ref(1);
    timeSolver = assignFields(timeSolver, 'dt', dt, 'nstep', nstep, 'theta', theta, 'USEGPU', USEGPU);

    %smat = real(dotSysmat (param.hMesh, param.mua, param.mus, param.ref, 0));
    smat = dotSysmat2_noC(param.hMesh, param.mua, param.mus, param.ref);
    mmat = param.hMesh.Massmat; % the bmat is only the boundary part!

    K0 = -(smat * (1-theta) - mmat * 1/(c*dt));            % backward difference matrix
    K1 = smat * theta + mmat * 1/(c*dt);                   % forward difference matrix

    
    % Compute inverse / factorization
    if ~timeSolver.USEGPU
      [timeSolver.L, timeSolver.U, timeSolver.P, timeSolver.Q, timeSolver.R] = lu(K1);
      timeSolver.K0 = K0;   %AF added because of an error if USEGPU=0
      %[L,U] = lu(K1);
    else
      %Select GPU (this should happen at higher level), 
      %print info abour GPU devices, this resets all memory too.
      timeSolver.gpu = gpuDevice;%(timeSolver.USEGPU);
      timeSolver.K0 = gpuArray((K0));
      timeSolver.K1 = gpuArray(K1);
    end
    
    tPhi = [];
    
    
  case 'apply'
       
    nh = size(phi,1);
    ns = size(phi,2);
    tPhi = zeros(nh,ns,timeSolver.nstep);

    % initial condition
    %tPhi(1,:,:) = phi;
    phi = phi/timeSolver.dt * 1000;
    tPhi(:,:,1) = phi;
    
    % AF initial conditions
    if ~timeSolver.USEGPU
        %tPhi(:,:,1) = timeSolver.Q*(timeSolver.U\(timeSolver.L\(timeSolver.P*(timeSolver.R\phi))));
        
        spacer = ' ------- ';
        disp([spacer mfilename ': time stepping' spacer]);
        dstep = ceil(timeSolver.nstep/4);
        
        tstart = cputime;
        
        for i = 2:timeSolver.nstep
            q = timeSolver.K0 * tPhi(:,:,i-1);
            %q = timeSolver.K0 * tPhi(:,:,i);
            
            % % standard LU
            %tPhi(:,:,i) = U\(L\q);
            % LU umfpack: P*(R\A)*Q = L*U
            tPhi(:,:,i) = timeSolver.Q*(timeSolver.U\(timeSolver.L\(timeSolver.P*(timeSolver.R\q))));
            
            ddisp(['Step: ' num2str(i) '. Elapsed time: ' num2str(cputime - tstart)],...
                mod(i,dstep)==0, mfilename);
            %    phi = max(phi,0);
        end
        %tPhi(:,:,1) = [];
        
    else
        spacer = ' ------- ';
        disp([spacer mfilename ': time stepping' spacer]);
        dstep = ceil(timeSolver.nstep/4);
        
        tstart = cputime;
        
        for is = 1:ns
            phi_g = gpuArray(phi(:,is));
            [phiG(:,is),~] = pcg(timeSolver.K1,phi_g,tol,100);
        end
       tPhi(:,:,1) = gather(phiG);
        for i = 2:timeSolver.nstep
            q = timeSolver.K0 * tPhi(:,:,i-1);
            for is = 1:ns
                phi_g = gpuArray(q(:,is));
                [phiG(:,is),~] = pcg(timeSolver.K1,phi_g,tol,100);
            end
            tPhi(:,:,i) = gather(phiG);
            
            ddisp(['Step: ' num2str(i) '. Elapsed time: ' num2str(cputime - tstart)], ...
                mod(i,dstep)==0, mfilename);
        end
        %tPhi(:,:,1) = [];
    end
    tPhi = permute(tPhi, [3 1 2])/1000;
    if (sum(tPhi(:)<0) > 0)
        warning([mfilename,spacer,num2str(sum(tPhi(:)<0)),' elements < 0 set to 0']);
        tPhi(tPhi<0) = 0;
    end
        
    
        
        
        
        
        
% % %     %%% here if GPU was in the the loop
% % %     
% % %     
% % %     if ~timeSolver.USEGPU
% % %         tPhi(:,:,1) = timeSolver.Q*(timeSolver.U\(timeSolver.L\(timeSolver.P*(timeSolver.R\phi))));
% % %     else
% % %         for is = 1:ns
% % %             phi_g = gpuArray(phi(:,is));
% % %             [phiG(:,is),flag] = pcg(timeSolver.K1,phi_g,tol,100);
% % %         end
% % %     end
% % %     tPhi(:,:,1) = gather(phiG);
% % %     %Put the right hand side and K0 on GPU too
% % % %     if timeSolver.USEGPU
% % % %       phi = gpuArray(full(phi));
% % % %        q = zeros(size(phi));
% % % %       q = gpuArray(zeros(size(phi)));
% % % %     end
% % %       
% % % 
% % %     t = [1:timeSolver.nstep] * timeSolver.dt;
% % % 
% % %     spacer = ' ------- ';
% % %     disp([spacer mfilename ': time stepping' spacer]);
% % %     dstep = ceil(timeSolver.nstep/4);
% % % 
% % % 
% % %     % loop over time steps
% % %     tstart = cputime;
% % %         
% % %     %for i = 2:timeSolver.nstep+1
% % %     for i = 2:timeSolver.nstep
% % %       
% % %       if ~timeSolver.USEGPU
% % %         q = timeSolver.K0 * tPhi(:,:,i-1);
% % %         %q = timeSolver.K0 * tPhi(:,:,i);
% % %         
% % %         % % standard LU
% % %         %tPhi(:,:,i) = U\(L\q);
% % %         % LU umfpack: P*(R\A)*Q = L*U
% % %         tPhi(:,:,i) = timeSolver.Q*(timeSolver.U\(timeSolver.L\(timeSolver.P*(timeSolver.R\q))));
% % %         
% % %       else
% % %         q = timeSolver.K0 * tPhi(:,:,i-1);
% % %         for is = 1:ns
% % %             phi_g = gpuArray(q(:,is));
% % %             [phiG(:,is),flag] = pcg(timeSolver.K1,phi_g,tol,100);
% % %         end
% % %         tPhi(:,:,i) = gather(phiG);
% % %         %phi = timeSolver.iK1*q;
% % %          
% % %       end
% % %       
% % %       ddisp(['Step: ' num2str(i) '. Elapsed time: ' num2str(cputime - tstart)], mod(i,dstep)==0, mfilename);
% % %       %    phi = max(phi,0);
% % %     end
% % %     %tPhi(:,:,i) = gather(phi);
% % %        
% % %     %tPhi(:,:,1) = 0;
% % %     %tPhi(:,:,1) = [];
% % %     tPhi = permute(tPhi, [3 1 2]);
% % %     
% % %    % if timeSolver.USEGPU    %AF
% % %    %     q = gather(q);
% % %    % end
% % %    % clear q;
% % %       
    
  case 'finalize'
       
    disp('Finalizing timeSolver')
    if timeSolver.USEGPU
      %phi = gather(phi);
      %q = gather(q); 
      K0 = gather(timeSolver.K0);
     % iK1 = gather(timeSolver.iK1);
      
      %Reset GPU device, because of memory leak
      %reset(timeSolver.gpu);
      gpu.delete;
    end
      
    clear timeSolver
end  


end