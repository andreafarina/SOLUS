function genExample(ExampleFolder,SOLUSfile,EXPsavename,priorname,flag_segm, varargin)
% genExample(ExampleFolder,SOLUSfile,EXPsavename,priorname,flag_segm,I,IO,QM,R,E) 
% generates initialisation files for the reconstruction of multimodal data
% from SOLUS
% - ExampleFolder: folder where to write the Initialisation files
% - SOLUSfile: SOLUS file used for reconstruction
% - EXPsavenane : name of output file containing the EXP structure
% - priorname : name of output file containing the prior
% - flag_segm: regulates the behaviour of the function with respect to the
% prior extraction. Run "help genEXP_fromSOLUS" for details
% - I : optional cell in form {'param1', param1,'param2', param2,...} that
% can be used to change the default of the file Init_DOT. Run
% "genInit('help',[],[])" for options on the default and on the possible
% parameters
% - IO : optional cell in form {'param1', param1,'param2', param2,...} that
% can be used to change the default of the file SetIO_DOT. Run
% "genIO('help',[],[])" for options on the default and on the possible
% parameters
% - QM : optional cell in form {'param1', param1,'param2', param2,...} that
% can be used to change the default of the file SetQM_DOT. Run
% "genQM('help',[],[])" for options on the default and on the possible
% parameters
% - R : optional cell in form {'param1', param1,'param2', param2,...} that
% can be used to change the default of the file RecSettings_DOT. Run
% "genRecSetts('help',[],[])" for options on the default and on the possible
% parameters
% - E : optional cell in form {homofile, hetefile, log1,log2, nom_vals{5}} that
% can be used to change the default of the experimental files EXP such as which optical data to reconstruct from. Run
% "help genEXP_fromSOLUS" for options






disp(['Creation of Example from file \n', SOLUSfile, ' in Folder\n',ExampleFolder])
% creating folder
disp('Folder Creation')
mkdir(ExampleFolder)

if isempty(varargin)
    for i = 1:5
        config{i} =[];
    end
    
end
for i = 1:numel(varargin)
   if ~isempty(varargin{i})
       config{i} = varargin{i}; 
   else
       config{i}=[];
   end
end
 

cfol = [getenv('DOTSRC'),filesep,'SOLUSwrapper',filesep];

% creating
disp('Generating: ')

disp('-EXP structure')
[dt,nstep]=genEXP_fromSOLUS([ExampleFolder,filesep,EXPsavename],[ExampleFolder,filesep,priorname],SOLUSfile,flag_segm,config{1});
disp('-Init')
genInit([cfol,'config_file_Init'], ExampleFolder,cat(2,{'dt',dt,'nstep',nstep},config{2})) 
disp('-IO')
genIO([cfol,'config_file_IO'], ExampleFolder,cat(2,{'Tomoname',EXPsavename},config{3}))
disp('-QM')
genQM([cfol,'config_file_QM'], ExampleFolder,config{4})
disp('-REC')
genRecSetts([cfol,'config_file_Rec'], ExampleFolder,cat(2,{'prior_path',priorname},config{5}))

disp('-Dummy Override')
id_out = fopen([ExampleFolder,filesep,'Override_MultiSim.m'],'w');
fprintf(id_out,'');
fclose(id_out);



end
