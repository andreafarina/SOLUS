function genRecSetts(config_file,output_folder,to_Update)
% genRecSetts(config_file,output_folder,to_Update)
% genRecSetts generate and initialization file RecSettings_DOT in the
% output folder output_folder starting from a template for RecSettings.
% config_file: is the path to a .mat filed containing two variables config
% and config_name. They are both monodimensional arrays of cells of the
% same dimensions. If the config_file is empty a predefined file is loaded. 
% If config_file is set to "help" a template config_name is shown.
% output_folder: output directory
% to_Update: an array of cells in the form {name1,value1,name2,values2,...}
% that is used to update the default values in config_file
% 
% % run genRecSetts('help',[],[]) to get help


% Load template
id_template = fopen('Template_RecSettings.txt','r');
str = fscanf(id_template,'%c');
fclose(id_template);

% load config
if ~isempty(config_file) && strcmpi(config_file,'help') == 0
    load(config_file)
elseif strcmpi(config_file,'help') == 1
    disp('Displaying suggested parameters')
    load('config_file_Rec')
    for i = 1:numel(CONFIG),disp([ CONFIG_NAME{i}]),disp('      '),disp(CONFIG{i}),end 
    return
else

     load('config_file_Rec')    
%     CONFIG = {'td';'linear';'auto';0.8;0.1;20;'{ ''mua'',''mus''}' ; 0.01;1;1.4;...
%         '[ 0.0091189   0.0082529   0.0074575    0.011498   0.0073212   0.0072694   0.0091525   0.0076381]';'[ 1.4267      1.3066     0.94913     0.78723     0.77073     0.70609     0.63522     0.59899]';...
%         0.78;1.05;'[9.50498	0.92704	0.15076	0.63606	1.11337]';0.1;'true';'false';'linear';'USprior';
%'';'false';
%'prejacobian/J';'spectral'};
%     CONFIG_NAME = {'domain','type_fwd','roi','range1','range2','NUM_TW','variables','mua0','musp0','nB',...
%         'mua','musp','a','b','conc','tau','fit_reference','fit_reference_far','fit_referecence_fwd','solver_type',...
%         'prior_path','prejacobian_load','prejacobian_path','fit_type'};
%      save('config_file_Rec','CONFIG','CONFIG_NAME')
end

%update config with userdefined parameters
CONFIG = update_config(CONFIG,CONFIG_NAME, to_Update);


wr_str = sprintf(str,CONFIG{:});

id_out = fopen([output_folder,filesep,'RecSettings_DOT.m'],'w');
fprintf(id_out,'%s',wr_str);
fclose(id_out);



end