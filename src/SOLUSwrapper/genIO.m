function genIO(config_file,output_folder,to_Update)
% genIO(config_file,output_folder,to_Update)
% genIO generates an initialization file SetIO_DOT in the
% output folder output_folder starting from a template for IO.
% config_file: is the path to a .mat filed containing two variables config
% and config_name. They are both monodimensional arrays of cells of the
% same dimensions. If the config_file is empty a predefined file is loaded. 
% If config_file is set to "help" a template config_name is shown.
% output_folder: output directory
% to_Update: an array of cells in the form {name1,value1,name2,values2,...}
% that is used to update the default values in config_file
% 
% run genIO('help',[],[]) to get help

% Load template
id_template = fopen('Template_SetIO.txt','r');
str = fscanf(id_template,'%c');
fclose(id_template);

% load config
if ~isempty(config_file) && strcmpi(config_file,'help') == 1
    load(config_file)
elseif strcmpi(config_file,'help') == 1
    disp('Displaying suggested parameters')
    load('config_file_QM')
    for i = 1:numel(CONFIG),disp([ CONFIG_NAME{i}]),disp('      '),disp(CONFIG{i}),end 
    return
else
    disp('loading template configuration')
    load('config_file_IO')  
    % CONFIG = {'Tomo.mat';'Test_Standard';'202109',;'./';'./';'spectra_polimi'};
    % CONFIG_NAME = {'exp_file','filename','session','exp_path','res_path','spectra_file'};
    % save('config_file_IO','CONFIG','CONFIG_NAME')
end




%update config with userdefined parameters
CONFIG = update_config(CONFIG,CONFIG_NAME, to_Update);


wr_str = sprintf(str,CONFIG{:});

id_out = fopen([output_folder,filesep,'SetIO_DOT.m'],'w');
fprintf(id_out,'%s',wr_str);
fclose(id_out);



end