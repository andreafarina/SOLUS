function genExample(ExampleFolder,SOLUSfile,EXPsavename,priorname,flag_segm, varargin)

% 
%
%

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
genRecSetts([cfol,'config_file_Rec'], ExampleFolder,cat(2,'prior_path',priorname,config{5}))

disp('-Dummy Override')
id_out = fopen([ExampleFolder,filesep,'Override_MultiSim.m'],'w');
fprintf(id_out,'');
fclose(id_out);



end
