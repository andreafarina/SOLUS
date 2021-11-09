function genQM(config_file,output_folder,to_Update)
% genQM(config_file,output_folder,to_Update)
% genQM generates an initialization file setQM_DOT in the
% output folder output_folder starting from a template for QM.
% config_file: is the path to a .mat filed containing two variables config
% and config_name. They are both monodimensional arrays of cells of the
% same dimensions. If the config_file is empty a predefined file is loaded. 
% If config_file is set to "help" a template config_name is shown.
% output_folder: output directory
% to_Update: an array of cells in the form {name1,value1,name2,values2,...}
% that is used to update the default values in config_file
% 

% Load template
id_template = fopen('Template_SetQM.txt','r');
str = fscanf(id_template,'%c');
fclose(id_template);

% load config
if ~isempty(config_file) &&  strcmpi(config_file,'help') == 0
    load(config_file)
elseif strcmpi(config_file,'help') == 1
    disp('Displaying suggested parameters')
    load('config_file_QM')
    for i = 1:numel(CONFIG),disp([ CONFIG_NAME{i},CONFIG{i}]),end 
else
    disp('loading template configuration')
    load('config_file_QM')    
%      CONFIG = {'[14 -14]';'[20.0000    6.6667   -6.6667  -20.0000]';'0';'0.5';...
%          '[14 -14]';'[20.0000    6.6667   -6.6667  -20.0000]';'0';'0.5';'[]';'[]'};
%      CONFIG_NAME = {'ys','xs','zs','SA','yd','xd','zd','DA','del_detector','del_source'};
%      save('config_file_QM','CONFIG','CONFIG_NAME')
end

%update config with userdefined parameters
CONFIG = update_config(CONFIG,CONFIG_NAME, to_Update);


wr_str = sprintf(str,CONFIG{:});

id_out = fopen([output_folder,filesep,'SetQM_DOT.m'],'w');
fprintf(id_out,'%s',wr_str);
fclose(id_out);



end