%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a structure EXP containing all data of the experiment
% This structure is loaded in the reconstruction script
% Andrea Farina (CNR-Polimi)
% October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InitScript
addpath('../../src/experimental/');
folder_data = 'Dati phantom carne/VEAL+LARD';
addpath((['./Tomo structs/' folder_data]));
lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
%lambda = [620 670 740 800 910 1020 1050 1090];

for iw = 1:numel(lambda)
EXP_file = strcat('Tomo_wave_',num2str(lambda(iw)));
load(EXP_file)
% =========================================================================
%%                                Path
% =========================================================================
%EXP.path.data_folder = '/media/psf/Home/Documents/Politecnico/Work/Data/Compressive Imaging/';
%EXP.path.data_folder = '/Users/andreafarina/Documents/Politecnico/Work/Data/Compressive Imaging/';
output_folder = './';
output_suffix = '';%'_sum100';
% =========================================================================
%%                      Measurement info
% =========================================================================
% Measurement folder and name
EXP.path.day = '201807';
EXP.path.data_folder = T.path.data_folder;%['../../results/201710/dynamic phantom'];

EXP.path.file_name = T.path.file_name;
% Time resolved windows
EXP.spc.gain = T.gain;
EXP.spc.n_chan = T.Nchannel;
EXP.spc.factor = 50e3/T.Nchannel/T.gain;
%EXP.time.roi = [500,1500];
EXP.time.roi = T.roi;  % just for info
EXP.bkg.ch_start = T.bkg.ch_start;
EXP.bkg.ch_end = T.bkg.ch_stop;
% % time windows
%  a = 1:2001;
%  b = a + 0;
%  twin1 = [a', b'];
%  twin = twin1;
% EXP.time.twin = twin;
% irf file
EXP.path.irf_file = 'IRF.mat';
% irf shift
EXP.irf.t0 = T.irf.t0; %40;%-15;%-30;   % ps
% sum repetitions

% =========================================================================
%%                         Pre-process SPC and CCD
% =========================================================================
[homo,hete,t,EXP.irf] = PreProcessing_SPC(EXP);

% =========================================================================
%%                         Determine size of data
% =========================================================================
spc_size = size(hete);
ref_size = size(homo);


% =========================================================================
%%                            Save forward data
% =========================================================================
% EXP.data.spc = int32(reshape(hete,size(hete,1),[]));
% EXP.data.ref = int32(reshape(homo,size(homo,1),[]));
EXP.data.spc = (reshape(hete,size(hete,1),[]));
EXP.data.ref = (reshape(homo,size(homo,1),[]));
EXP.time.axis = t;
EXP.grid.dmask = T.dmask;
EXP.grid.xs = T.xs;
EXP.grid.ys = T.ys;
EXP.grid.zs = T.zs;
EXP.grid.SourcePos = T.sourcepos;
EXP.grid.DetPos = T.detpos;
EXP.optp.homo.abs = T.opt.homo.abs;
EXP.optp.homo.sca = T.opt.homo.sca;
EXP.optp.hete.abs = T.opt.hete.abs;
EXP.optp.hete.sca = T.opt.hete.sca;
EXP.lambda = T.lambda;
mkdir(fullfile(pwd,'EXP structs',folder_data))
str_file = [fullfile(pwd,'EXP structs',folder_data),filesep,...
    'EXP_',EXP_file,output_suffix];
save(str_file,'EXP');
%%assemble full exp struct
EXP_full.irf.area(iw) = EXP.irf.area;
EXP_full.irf.baric(iw) = EXP.irf.baric;
EXP_full.irf.data(:,iw) =  EXP.irf.data;
EXP_full.irf.peak(iw) = EXP.irf.peak;
EXP_full.irf.variance(iw) = EXP.irf.variance;
EXP_full.data.ref(:,(1:size(EXP.data.ref,2))+(iw-1)*size(EXP.data.ref,2)) = EXP.data.ref;
EXP_full.data.spc(:,(1:size(EXP.data.spc,2))+(iw-1)*size(EXP.data.spc,2)) = EXP.data.spc;
EXP_full.path = EXP.path; EXP_full.path.file_name = 'Tomo';
EXP_full.spc = EXP.spc;EXP_full.time = EXP.time;
EXP_full.bkg = EXP.bkg;EXP_full.grid = EXP.grid;
EXP_full.optp = EXP.optp; EXP_full.lambda = EXP.lambda;
end
clear EXP
EXP = EXP_full;
str_file = [fullfile(pwd,'EXP structs',folder_data),filesep,...
    'EXP_','Tomo',output_suffix];
save(str_file,'EXP');
