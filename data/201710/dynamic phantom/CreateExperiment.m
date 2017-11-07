%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a structure EXP containing all data of the experiment
% This structure is loaded in the reconstruction script
% Andrea Farina (CNR-Polimi)
% October 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clearvars;
addpath('../../../src/experimental/');
EXP_file = 'DYN_001_5rep';%'CYL_29_fit';%_noPN';
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
EXP.path.day = '201710';
EXP.path.data_folder = pwd;%['../../results/201710/dynamic phantom'];

EXP.path.data_file_spc = 'all_data';
% Time resolved windows
EXP.spc.gain = 4;
EXP.spc.n_chan = 4096;
EXP.spc.factor = 50e3/EXP.spc.n_chan/EXP.spc.gain;
%EXP.time.roi = [500,1500];
EXP.time.roi = [500,2500];  % just for info
EXP.bkg.ch_start = 300;
EXP.bkg.ch_end = 350;
% % time windows
%  a = 1:2001;
%  b = a + 0;
%  twin1 = [a', b'];
%  twin = twin1;
% EXP.time.twin = twin;
% irf file
EXP.path.irf_file = 'IRF.mat';
% irf shift
EXP.irf.t0 = 30; %40;%-15;%-30;   % ps
% sum repetitions
rep_int = 1:5;
% =========================================================================
%%                         Pre-process SPC and CCD
% =========================================================================
[spc,t,EXP.irf] = PreProcessing_SPC(EXP,EXP.path.data_file_spc);
ref = sum(spc(:,:,:,1,rep_int),5);   % homo
spc = sum(spc(:,:,:,2,rep_int),5);   % hetero
EXP.irf.data = sum(EXP.irf.data(:,rep_int),2);
% =========================================================================
%%                         Determine size of data
% =========================================================================
spc_size = size(spc);
ref_size = size(ref);


% =========================================================================
%%                            Save forward data
% =========================================================================
EXP.data.spc = int32(reshape(spc,size(spc,1),[]));
EXP.data.ref = int32(reshape(ref,size(ref,1),[]));
EXP.time.axis = t;
str_file = [pwd,filesep,...
    'EXP_',EXP_file,output_suffix];
save(str_file,'EXP');

