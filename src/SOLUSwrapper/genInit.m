function genInit(config_file,output_folder,to_Update)
% genInit(config_file,output_folder,to_Update)
% genInit generates an initialization file Init_DOT in the
% output folder output_folder starting from a template for Init.
% config_file: is the path to a .mat filed containing two variables config
% and config_name. They are both monodimensional arrays of cells of the
% same dimensions. If the config_file is empty a predefined file is loaded. 
% If config_file is set to "help" a template config_name is shown.
% output_folder: output directory
% to_Update: an array of cells in the form {name1,value1,name2,values2,...}
% that is used to update the default values in config_file
% 
% % run genInit('help',[],[]) to get help

switch nargin 
    case 0
        config_file =[];
        output_folder = './';
    case 1
        output_folder = './';
end
% Load template
id_template = fopen('Template_Init_DOT.txt','r');
str = fscanf(id_template,'%c');
fclose(id_template);

% load config
if ~isempty(config_file) && strcmpi(config_file,'help') == 0
    load(config_file)
    
elseif strcmpi(config_file,'help') == 1
    disp('Displaying suggested parameters')
    load('config_file_Init')
    for i = 1:numel(CONFIG),disp([ CONFIG_NAME{i}]),disp('      '),disp(CONFIG{i}),end 
    return
else
    disp('loading template configuration')
     load('config_file_Init')    
% %     CONFIG = {1;1;1;'linear';'semi-inf';1;1;1;'all';1;1;...
% %         1;'Born';1;1;0;0;0;...
% %         '{''Hb'',''HbO2'',''Lipid'',''H2O'',''Collagen''}';'[1,1,1,1,1]';'[1,1,10*100*0.91,10*100, 10*100*0.196]';' {''microM'',''microM'',''mg/cm^3'',''mg/cm^3'',''mg/cm^3''}';0;...
% %         '1:8';'[ 0.0091189   0.0082529   0.0074575    0.011498   0.0073212   0.0072694   0.0091525   0.0076381]';'[ 1.4267      1.3066  0.94913  0.78723  0.77073  0.70609  0.63522  0.59899]'; 0.79; 1.06 ; '[9.50498 0.92704	0.15076	0.63606	1.11337]';...
% %         0.01;1.0;1.4;1;-31;31;2;-29;29;2;0;60;2;' {''Mua'',''Musp''} ';'cylinder';'[0,0,-1]';15;11;'[0, 0, 15+ DOT.opt.hete1.l]';'OFF';'Step';...
% %         '[ 0.0052469   0.0048377   0.0036964   0.0080475   0.0035186   0.0033081   0.0052744   0.0038201]';'[ 1.1469  1.055  0.7657  0.615  0.60022  0.55903 0.48579  0.46547]';1.54;1;'[0.4085 0.71424 0.48844	0.5868	0.37051]';...
% %         '[''VICTRE.mat'']';91.7105;128;'Poisson';1e-3;true;1e6;1;1;0.9; '[640 675 830 905 930 970 1020 1050 ]'  ;4;0.05;620};
% %     CONFIG_NAME = {'SHOWPLOTS','FORWARD','REF','TYPE_FWD','geom','RECONSTRUCTION','EXPERIMENTAL','EXP_IRF','EXP_DELTA','EXP_DATA',...
% %         'TD','sigma','type','RADIOMETRY','SAVE_FWD','LOAD_FWD_TEO','TOAST2DOT','SPECTRA',...
% %         'cromo_label','active_cromo','cromo_factor' ,'cromo_units','ForceConstitSolution',...
% %         'lamda_id','mua_','musp_','a_','b_','conc_','muaB','muspB','nB','nE','x1','x2','dx','y1','y2','dy','z1','z2','dz',...
% %         'hete1_type','hete1_geom','hete1_d','hete1_l','hete1_sigma','hete1_c','hete1_distrib','hete1_profile',...
% %        'muap_','muspp_','ap_','bp','concp_','hete1_path',...
% %        'dt','nstep','noise','tsigma','self_norm','tot_counts','power','acqtime','opteff','lambda','area','qeff','lambda0'};
% %         save('config_file_Init','CONFIG','CONFIG_NAME')
 %      for i = 1:numel(CONFIG),disp(CONFIG_NAME{i}),disp(CONFIG{i}),end 
end

%update config with userdefined parameters
CONFIG = update_config(CONFIG,CONFIG_NAME, to_Update);


wr_str = sprintf(str,CONFIG{:});

id_out = fopen([output_folder,filesep,'Init_DOT.m'],'w');
fprintf(id_out,'%s',wr_str);
fclose(id_out);



end

% 
% FORWARD = %g;
% REF = %g;   
% TYPE_FWD = '%s';
% geom = %s; 
% RECONSTRUCTION = %g;
% EXPERIMENTAL = %g; 
% EXP_IRF = %g;       
% EXP_DELTA = '%s';  
% EXP_DATA = %g;      
% DOT.TD = %g;     
% DOT.sigma = %g;  
% type = '%s';
% RADIOMETRY = %g;       
% SAVE_FWD = %g;           
%                         
% LOAD_FWD_TEO = %g;         
% TOAST2DOT = %g;         
% SPECTRA = %g; 
% DOT.spe.cromo_label = %s; 
% DOT.spe.active_cromo = %s; 
% DOT.spe.cromo_factor = %s; 
% DOT.spe.cromo_units = %s; 
% DOT.spe.ForceConstitSolution =%g; 
% lamda_id = %s; 
% if SPECTRA == 0 
%     mua_=%s; 
%     musp_=%s;  
%     mua_ = mua_(lamda_id); 
%     musp_ = musp_(lamda_id); 
%     Xd = {mua_,musp_}; 
% else 
%     a_ = %g; 	
%     b_ =%g; 
%     conc_ = %s; 
%     Xd = {conc_,[a_ b_]}; 
% end 
% 
% DOT.opt.muaB = %g;  
% DOT.opt.muspB = %g;      
% DOT.opt.nB = %g;      
% DOT.opt.nE = %g;     
% 
% DOT.grid.x1 = %g;
% DOT.grid.x2 = %g;
% DOT.grid.dx = %g;
% 
% DOT.grid.y1 = %g;
% DOT.grid.y2 = %g;           
% DOT.grid.dy = %g;
% 
% DOT.grid.z1 = %g;        
% DOT.grid.z2 = %g;         
% DOT.grid.dz = %g;
% 
% NUM_HETE = %g;
% 
% DOT.opt.hete1.type  = '%s';  
% DOT.opt.hete1.geometry = '%s';  
% DOT.opt.hete1.d     = %g;  
% DOT.opt.hete1.l     = %g; 
% DOT.opt.hete1.sigma = %g; 
% DOT.opt.hete1.c     = %s; 
% DOT.opt.hete1.distrib = %g; 
% DOT.opt.hete1.profile ='%s'; 
% DOT.opt.hete1.val   = %s; 
% if SPECTRA == 0  
%     muap_=%s; 
%     muspp_=%s; 
%     Xp = {muap_,muspp_}; 
% else 
%     a_ = %g;	b_ = %g;  
%     concp_ = %s; 
%     Xp = {concp_,[a_ b_]}; 
% end 
% DOT.opt.hete1.path =[%s] ;  
% 
% DOT.time.dt =%g  
% dt = %s;
% DOT.time.nstep = %g;             
% nstep = %s;
% DOT.time.noise = '%s';        
%                                  
% DOT.time.sigma = %g;              
% DOT.time.self_norm = %g;        
% DOT.time.TotCounts = %g;          
%                                    
% DOT.radiometry.power = [%g]; 
% DOT.radiometry.timebin = DOT.time.dt; 
% DOT.radiometry.acqtime = %g;    
% DOT.radiometry.opteff = %g;   
% DOT.radiometry.lambda = %s;
% DOT.radiometry.lambda = DOT.radiometry.lambda(lamda_id);
% DOT.radiometry.nL = numel(DOT.radiometry.lambda);
%                                 
% DOT.radiometry.area = %g;       
% DOT.radiometry.qeff =%g;    
% DOT.radiometry.lambda0 = %g;
% if size(DOT.time.TotCounts,2)<DOT.radiometry.nL
%     DOT.time.TotCounts = repmat(DOT.time.TotCounts,1,DOT.radiometry.nL);
% end
% CUT_COUNTS = 1;         
% NumDelays = 1;    