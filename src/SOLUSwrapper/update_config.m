function out = update_config(config, config_name,varargin)
% out = update_config(config, config_name,varargin)
% This function overwrites the configuration in the variable config with
% those funished in vararging in the form of {Name1, newvalue1,Name2, newvalue2,...}
% config: monodimensional array of cells
% config_name: array of strings of the same length of config and contaning
% the name descriptors of the file config
% 
% 
    for i = 1:numel(config_name)
       for j =1:2:(numel(varargin{1}))
           if strcmpi(config_name{i},varargin{1}{j})==1
               config{i} =varargin{1}{j+1};
           end
       end
    end
    out = config;
    return
end