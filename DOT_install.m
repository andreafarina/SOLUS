% Add path and set environment variables.
StartVars = who;
disp('Setting environment variable DOTDIR');
setenv('DOTDIR',pwd);
disp('Setting environment variable DOTSRC');

setenv('DOTSRC',[getenv('DOTDIR'),filesep,'src']);
disp('Adding path');
addpath(getenv('DOTSRC'));
addpath([getenv('DOTSRC'),filesep,'util']);
d = dir([getenv('DOTSRC'),filesep,'util']);
names = {d.name};
isdir = {d.isdir};
num_names = size(names,2);
n_name = size(names,2);
for i=3:n_name
    f = names{i};
    if isdir{i}
        addpath([getenv('DOTSRC'),filesep,'util',...
            filesep,f]);
    end
end
addpath([getenv('DOTSRC'),filesep,'fwd']);
addpath([getenv('DOTSRC'),filesep,'subroutines']);
addpath([getenv('DOTSRC'),filesep,'solvers']);
addpath([getenv('DOTSRC'),filesep,'experimental']);
addpath([getenv('DOTSRC'),filesep,'visual']);
addpath([getenv('DOTSRC'),filesep,'multiSimProcedures']);
clearvars('-except',StartVars{:});



