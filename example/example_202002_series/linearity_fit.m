%% analyisis
%clear

direc =  { 'ex_bulk10_0.1_incl10_0.05',...
    'ex_bulk10_0.1_incl10_0.4',...
   'ex_bulk10_0.2_incl10_0.2',...
    'ex_bulk10_0.1_incl10_0.1',...
    'ex_bulk10_0.1_incl10_0.6',...
    'ex_bulk10_0.2_incl10_0.4',...   
    'ex_bulk10_0.1_incl05_0.1',...
    'ex_bulk10_0.1_incl10_0.2',...   
    'ex_bulk10_0.1_incl15_0.1' };
direc =  { 'ex_bulk10_0.1_incl10_0.05',...
    'ex_bulk10_0.1_incl10_0.4',...
    'ex_bulk10_0.1_incl10_0.1',...
    'ex_bulk10_0.1_incl10_0.6',...
    'ex_bulk10_0.1_incl10_0.2' };


% direc =  {'ex_bulk10_0.1_incl10_0.2_p0.0001',...
%     'ex_bulk10_0.1_incl10_0.2_p0.00025',...     
% 'ex_bulk10_0.1_incl10_0.2_p0.0005',...
% 'ex_bulk10_0.1_incl10_0.2_p0.001',...
% 'ex_bulk10_0.1_incl10_0.2_p0.0025',...     
% 'ex_bulk10_0.1_incl10_0.2_p0.005',...
% 'ex_bulk10_0.1_incl10_0.2_p0.01',...    
% 'ex_bulk10_0.1_incl10_0.2_p0.025',...               
% 'ex_bulk10_0.1_incl10_0.2_p0.05',...    
% 'ex_bulk10_0.1_incl10_0.2_p0.1'...
% 'ex_bulk10_0.1_incl10_0.2_p0.2'};
% 
% %load Test_Standard_REC.mat

iii = 1;
i_fig = 1001;
for i = 1:numel(direc)
    %['dtprior/',direc{i},'/Test_Standard_REC.mat']
    load( ['./dtprior_nonLinFit/',direc{i},'/Test_Standard_REC.mat'])
    for lambda = 1:8
        [bx, by, bz] = ndgrid(REC.grid.x, REC.grid.y, REC.grid.z);
        if lambda >0 %== 3  
        abs_real(iii) = REC.opt.hete1.val(lambda);  
        mua_recon = REC.opt.bmua(:,lambda);
        if numel(size(REC.solver.prior.refimage)) > 3
            REC.solver.prior.refimage = REC.solver.prior.refimage(:,:,:,1); 
            
        end
        refim = REC.solver.prior.refimage;
        abs_recon(iii) = mean(REC.opt.bmua(abs(bx(:)) <= 1 & abs(by(:)) <= 1 & abs(bz(:)) >= 9 & abs(bz(:)) <= 11, lambda));%mean(mua_recon(REC.solver.prior.refimage(:) >  mean(REC.solver.prior.refimage(:)))); %mean(REC.opt.bmua(abs(bx(:)) <= 1 & abs(by(:)) <= 1 & abs(bz(:)) > 6 & abs(bz(:)) < 10, lambda));
        abs_recon(iii) = mean(REC.opt.bmua(refim(:), lambda)); 
        %mean(mua_recon > 0.5*mean(mua_recon)); %
        abs_bulk(iii) = REC.opt.mua0(lambda);
        sca_real(iii) = REC.opt.hete1.val(lambda+ 8);
        mus_recon = REC.opt.bmusp(:,lambda);
        sca_recon(iii) = mean(REC.opt.bmusp(abs(bx(:)) <= 1 & abs(by(:)) <= 1 &  abs(bz(:)) >= 9 & abs(bz(:)) <= 11, lambda));%mean(mus_recon > 0.5*mean(mus_recon));%
        sca_bulk(iii) = REC.opt.musp0(lambda);
        tau(iii) = REC.solver.tau;
        iii = iii +1;
        if lambda == 0
            PlotMua = reshape(REC.opt.bmua,[REC.grid.dim REC.radiometry.nL]);
            PlotMus = reshape(REC.opt.bmusp,[REC.grid.dim REC.radiometry.nL]);
            PlotMua = PlotMua(:,:,:,lambda);
            PlotMus = PlotMus(:,:,:,lambda);
            figure(i_fig);ShowRecResults(REC.grid,PlotMua,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto',0.00,0.05);
                pause(1)
                    suptitle(sprintf('Recon Mua,target incl = %g mm-1, bulk = %g mm-1, lambda = %g nm, tau = %g',...
                        round(abs_real(iii-1),4),...
                        round(abs_bulk(iii-1),4),...
                        830,...
                        0*tau(iii-1)))
                    pause(1)
%             figure(i_fig+1000);ShowRecResults(REC.grid,PlotMus,...
%                     REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto',0.00,0.05);
%                     suptitle(sprintf('Recon Mus, ,nominal incl = %g , bulk = %g, lambda= %g',sca_real(iii-1), sca_bulk(iii-1), lambda))    
%                 
                i_fig = i_fig + 1;
        end
        end
    end
end
delta_abs_real = (abs_real)./abs_bulk;
delta_abs_recon = (abs_recon)./abs_bulk;

delta_sca_real = (sca_real)./sca_bulk;
delta_sca_recon = (sca_recon)./sca_bulk;

err_abs = (abs_recon - abs_real)./abs_real; 
err_sca = (sca_recon - sca_real)./abs_real;

figure(2);scatter((delta_abs_real), delta_abs_recon),
xlabel('nominal contrast \mu_{a,in}/ \mu_{a,bulk}'),
ylabel('Reconstructed contrast \mu_{a,in}^{recon} / \mu_{a,bulk}'),
title('ABS: Recontructed vs nominal contrast with DT prior, \tau = 0.05' )

figure(3);scatter(delta_sca_real, delta_sca_recon),
xlabel('nominal contrast \mu_{s,in}/ \mu_{s,bulk}'),
ylabel('Reconstructed contrast \mu_{s,in}^{recon} / \mu_{s,bulk}'),
title('SCA: Recontructed vs nominal contrast with DT prior, \tau = 0.05' )

figure(4);scatter(delta_abs_real, err_abs),
xlabel('nominal contrast \mu_{a,in}/ \mu_{a,bulk}'),
ylabel('Reconstruction Error (\mu_{a,in}^{recon} - \mu_{a,in}) / \mu_{a,in}'),
title('ABS: ERROR vs nominal contrast with DT prior, \tau = 0.05' )

figure(5);scatter(delta_sca_real, err_sca),
xlabel('nominal contrast \mu_{s,in}/ \mu_{s,bulk}'),
ylabel('Reconstruction Error (\mu_{s,in}^{recon} - \mu_{s,in}) / \mu_{s,in}'),
title('SCA: ERROR vs nominal contrast with DT prior, \tau = 0.05' )

figure(6);scatter(abs_real, abs_recon),
xlabel('nominal abs \mu_{a,in}'),
ylabel('Reconstruction \mu_{a,in}^{recon} '),
title('ABS: abs vs nominal abs with DT prior, \tau = 0.05' )

figure(7);semilogx(tau, err_abs, 'r--+',...
    'MarkerEdgeColor','r','MarkerFaceColor','r', 'LineWidth',2, 'MarkerSize', 15),
xlabel('Regularisation Parameter'),
ylabel('Reconstruction Error (\mu_{a,in}^{recon} -\mu_{a,in})/\mu_{a,in} '),
title('Absorption Reconstruction Error vs Regularisation Parameter' )

xx = reshape(delta_abs_real,[8,5]);
yy = reshape(delta_abs_recon, [8,5]);
for i_s = 1:5
    a(i_s,:) = rand(1,3);
figure(8);plot(xx(:,i_s), yy(:,i_s), '+', 'MarkerSize', 8, 'color', a(i_s,:)), hold on
%plot(xx(:,i_s), yy_old(:,i_s), 'o', 'MarkerSize', 12, 'color', a(i_s,:)), hold on
end
pause(2);

grid on
xlabel('nominal contrast \mu_{a,in}/ \mu_{a,bulk}'),
ylabel('Reconstructed contrast \mu_{a,in}^{recon} / \mu_{a,bulk}'),
legend('\mu_{a}  = 0.005 mm^{-1}','\mu_{a}  = 0.01 mm^{-1}', '\mu_{a}  = 0.02 mm^{-1}', '\mu_{a}  = 0.04 mm^{-1}', '\mu_{a}  = 0.06 mm^{-1}')
title('ABS: Recontructed Contrast vs target contrast with DT prior, \tau = 0.01 (+), \mu_{s,bulk} \sim 1 mm^{-1}  \mu_{a,bulk} \sim 0.01 mm^{-1}' )
for i_s = 1:5
 %   a(i_s,:) = rand(1,3);
figure(8);%plot(xx(:,i_s), yy(:,i_s), '+', 'MarkerSize', 12, 'color', a(i_s,:)), hold on
plot(xx(:,i_s), yy_old(:,i_s), 'o', 'MarkerSize', 8, 'color', a(i_s,:)), hold on, legend off
end
legend('\mu_{a}  = 0.005 mm^{-1}','\mu_{a}  = 0.01 mm^{-1}', '\mu_{a}  = 0.02 mm^{-1}', '\mu_{a}  = 0.04 mm^{-1}', '\mu_{a}  = 0.06 mm^{-1}')

disp('END')
function [a,b,c] = fit_linearity(vec_real, vec_recon)

    vec_real = vec_real/ 0.01
    a = 1;
    b = 1;
    c = 1;

    f = fit(vec_real, vec_recon, 'poly2');
    return 

end
