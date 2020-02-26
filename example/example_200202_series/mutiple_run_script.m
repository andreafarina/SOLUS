%% run series of experiments


direc =  { 'ex_bulk10_0.1_incl10_0.05',...
    'ex_bulk10_0.1_incl10_0.4',...
    'ex_bulk10_0.2_incl10_0.2',...
    'ex_bulk10_0.1_incl10_0.1',...
    'ex_bulk10_0.1_incl10_0.6',...
    'ex_bulk10_0.2_incl10_0.4',...   
    'ex_bulk10_0.1_incl05_0.1',...
    'ex_bulk10_0.1_incl10_0.2',...   
    'ex_bulk10_0.1_incl15_0.1' };

% 
% for i = 1:numel(direc)
%     close all
%     old_pwd = pwd;
%     run_([old_pwd,'/free/',direc{i}], pwd);
%     system('rm ./precomputed_jacobians/J*')
%     cd ../../
% end

% % 
for i = 1:numel(direc)
    close all
    old_pwd = pwd;
    run_([old_pwd,'/cylprior/',direc{i},], pwd);
    system('rm ./precomputed_jacobians/J*')
    cd ../../
end


for i = 1:numel(direc)
    close all
    old_pwd = pwd;
    run_([old_pwd,'/dtprior/',direc{i}], pwd);
    system('rm ./precomputed_jacobians/J*')
    cd ../../
end

% % pp = {'0.0001','0.00025','0.0005','0.001','0.0025', '0.005', '0.01', '0.025', '0.05', '0.1', '0.2'};
% % for i = 1:numel(pp)
% %     close all
% %     old_pwd = pwd;
% %     run_([old_pwd,'/cylprior_priortest/','ex_bulk10_0.1_incl10_0.2_p',pp{i}], pwd);
% %     system('rm ./precomputed_jacobians/J*')
% %     cd ../../
% % end
% % 
% % for i = 1:numel(pp)
% %     close all
% %     old_pwd = pwd;
% %     run_([old_pwd,'/dtprior_priortest/','ex_bulk10_0.1_incl10_0.2_p',pp{i}], pwd);
% %     system('rm ./precomputed_jacobians/J*')
% %     cd ../../
% % end




function run_(cd_direc, varargin)
    
    cd(cd_direc);
    run_DOT;

end





