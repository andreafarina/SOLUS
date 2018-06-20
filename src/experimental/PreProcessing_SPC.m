 function [spc,time,irf] = PreProcessing_SPC(lprm,data_file_spc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process SPC raw data created with SaveCCD_SPC_Data.m
% 
% Allow to define a ROI and to subtract a background.
% Allow the gating of the SPC curves.
% It is possible to plot data as matrix.
% New dataset are created
% INPUT: parameter structure
% OUTPUT: spc data, time scale, irf structure containing the curve, peak,
% baricenter.
% July 2016: IRF loaded as raw mat file for bkg subtraction
% ROI and Bkg_ROI are absolute with respect to the whole timescale
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% =========================================================================
%%                          Options
% =========================================================================
SHOW_CURVES = 1;         % Show TR raw curves
POISS_DENOISE = 0;       % Anscombe transform and wavelet denoising
BKG_CONST = 1;           % subtract constant background for each measurements
APPLY_ROI = 0;           % decide whether or not apply ROI
APPLY_TWIN = 0;          % apply twin binning
SAVE = 0;                % save into a new data_file

%if nargin<1 DEFAULT PARAMETERS
% =========================================================================
%%                                Path
% =========================================================================
data_folder = '/Users/andreafarina/Documents/Politecnico/Work/Data/Compressive Imaging/';
% =========================================================================
%%                      Measurement info
% =========================================================================
% SPC settings
Gain = 3;
n_chan = 4096;
% Measurement folder and name
day = '201510';
%data_file_spc = 'CYL_19_SPC_raw';
%irf_file = 'IRF2_test.sdt';
% Constant background region
ch_bkg_start = 100;
ch_bkg_end = 150;
% IRF shift
t0 = 0;%-15;-30;       % ps
% =========================================================================
%%                          Region of interest
% =========================================================================
roi = [1,4096];   % channels
if nargin < 1
    disp('ATTENTION! No experiement specified! Using default');
end
% =========================================================================
%%                  DEFINE TEMPORAL WINDOWS
% =========================================================================
factor = 50e3/n_chan/Gain;               % ps/ch
% Manual definition of gates referred to the raw data channel numbers
%twin = [1,1000;1001,2000;2001,3000;3001,4000];
% Constant gates
% a = 300:100:2499;
% b = a + 99;
% twin1 = [a', b'];
% c = 1:(roi(2) - roi(1) +1);
c = 1:n_chan;
twin = [c',c'];
%end
    

if nargin>0

% =========================================================================
%%                        Read local paramenters
% =========================================================================
    % SPC settings
    if isfield(lprm,'spc')
        if isfield(lprm.spc,'gain')
            Gain = lprm.spc.gain;
        end
        if isfield(lprm.spc,'n_chan')
            n_chan = lprm.spc.n_chan;
        end
        if ((isfield(lprm.spc,'gain'))&&(isfield(lprm.spc,'n_chan')))
            factor = 50e3/lprm.spc.n_chan/lprm.spc.gain;
        end
    end
    % Path
    if isfield(lprm,'path')
        if isfield(lprm.path,'data_folder')
            data_folder = lprm.path.data_folder;
        end
        if isfield(lprm.path,'day')
            day = lprm.path.day;
        end
%         if isfield(lprm.path,'data_file_spc')
%             data_file_spc = lprm.path.data_file_spc;
%         end
        if isfield(lprm.path,'irf_file')
            irf_file = lprm.path.irf_file;
        end
    else
        disp('ATTENTION! No data file specified! Using default');
    end
    % Background
    if isfield(lprm,'bkg')
        if isfield(lprm.bkg,'use_measure')
            BKG_MEAS = lprm.bkg.use_measure;
        end
        if isfield(lprm.bkg,'const')
            BKG_CONST = lprm.bkg.const;
        end
        if isfield(lprm.bkg,'ch_start')
            ch_bkg_start = lprm.bkg.ch_start;
        end
        if isfield(lprm.bkg,'ch_end')
            ch_bkg_end = lprm.bkg.ch_end;
        end
    end
    if isfield(lprm,'time')
        if isfield(lprm.time,'roi')
            roi = lprm.time.roi;
            %% Adjust time window to ROI
            c = 1:(roi(2) - roi(1) +1);
            twin = [c',c'];
        end
        if isfield(lprm.time,'twin')
            twin = lprm.time.twin;
        end
    end
    % irf
    if isfield(lprm,'irf')
        if isfield(lprm.irf,'t0')
            t0 = lprm.irf.t0;
        end
    end
    
    %
end
% =========================================================================
%%                               BEGIN
% =========================================================================
%% Create ROI index
roi_vec = roi(1):roi(2);
% Check consistency between ROI and twin
if numel(roi_vec)<max(twin(:))
    disp('ROI and time windows not consistent!')
    spc = -1;
    time = -1;
    irf = -1;
    return;
end
%% Construct path
%data_path = [data_folder,day,filesep];
data_path = data_folder;
%% load data
load([data_path,filesep,data_file_spc]);
spc = all_data;
clear all_data
%% Determine the dataset format
s = size(spc);
spc = reshape(spc,s(1),[]);
if POISS_DENOISE==1
    spcf = zeros(size(spc));
    for i=1:size(spc,2)
        spcf(:,i) = AnscombeWaveletDenoising_AF(spc(:,i));
    end
    spc = spcf;
    clear spcf
end

    
%% load IRF 
if nargout > 2
    tmp  = load([data_path,filesep,irf_file]);
    irf.data = tmp.IRF;
    clear tmp
else
    irf.data = zeros(s(1),1);
    irf.data(1) = 1;
end

if POISS_DENOISE==1
    for i = 1:size(irf.data,2)
        irf.data(:,i) = AnscombeWaveletDenoising(irf.data(:,i));
    end
end
%% Background constant value
if BKG_CONST == 1
    %bkg=mean(spc(ch_bkg_start:ch_bkg_end,:,:,:));
    bkg = repmat(mean(spc(ch_bkg_start:ch_bkg_end,:)),n_chan,1);
    spc = spc-bkg;
end
%% exclude negative values
%spc(spc<0) = 0;
% =========================================================================
%%                      Temporal axis definition
% =========================================================================
time = (1:size(spc,1))*factor - factor/2;
%t = (t*factor-factor/2);
% =========================================================================
%%                              Apply ROI
% =========================================================================
if APPLY_ROI == 1
    %t = t(roi_vec);
    time = (1:length(roi_vec))*factor-factor/2;
    spc = spc(roi_vec,:);
    s(1) = size(spc,1);
end
% =========================================================================
%%                            Plot TR data
% =========================================================================
if SHOW_CURVES == 1
    spc = reshape(spc,s);
    spc_sum = sum(spc,5);
    spc_sum_resh = reshape(spc_sum,n_chan,[]);
    %hold on
    for i = 1:size(spc_sum_resh,2)
        %semilogy(time,[spc(:,i),irf.data(:,1)]),ylim([1 65000]),grid,xlabel('time (ps)')
        semilogy([spc_sum(:,i),sum(irf.data,2)]),ylim([1 65000]),grid,xlabel('chan')
        
        ylabel('Counts'),
        pause(0.1);
    end
end
% =========================================================================
%%                         Temporal windows
% =========================================================================
if APPLY_TWIN == 1
    nwin = size(twin,1);
    size_y = size(spc);
    size_y(1) = nwin;
    y = zeros(size_y);
    
    for w = 1:nwin
        y(w,:) = sum(spc(twin(w,1):twin(w,2),:),1);
    end
    time = (mean(twin,2)*factor-factor/2);% + (roi(1)-1) * factor;
    % DEBUG
    %semilogy(t,spc(:,1,10,9))
    %semilogy(t,spc(:,3,10,9)./max(spc(:,3,10,9)),time,y(:,3,10,9)./max(y(:,3,10,9)))
    %semilogy(t,spc(:,3,10,9),time,y(:,3,10,9)),legend('spc','gated')
    spc = y;
    clear y
end
% =========================================================================
%%                         Prepare dataset to be saved
% =========================================================================
spc = reshape(spc,[n_chan,s(2:end)]);
% =========================================================================
%%                          Save new data
% =========================================================================
if SAVE == 1
    save([data_path, data_file_spc(1:end-4)],'spc','time');
end

% =========================================================================
%%                          Process IRF
% =========================================================================
if nargout > 2
%    irf.data = double(f_read_sdt_01([data_path,irf_file]));
    figure;semilogy(irf.data),title('IRF')
    % Background subtraction (same channels of the measruements)
    if BKG_CONST == 1
        %bkg=mean(spc(ch_bkg_start:ch_bkg_end,:,:,:));
        bkg = mean(irf.data(ch_bkg_start:ch_bkg_end,:));
        irf.data = irf.data-repmat(bkg,n_chan,1);
    end
    % apply roi
    if APPLY_ROI == 1
        irf.data = irf.data(roi_vec,:);
    end
    % shift irf
    delta = round(t0./factor);
    irf.data = circshift(irf.data,delta);
    % apply ROI
    %irf.data = irf.data(roi_vec,:);
    %% calculate indicator summing up all the IRF repetitions
    irf_sum = sum(irf.data,2);
    % Calculate peak position
    [irf.peak.value,irf.peak.pos] = max(irf_sum);
    irf.peak.time = time(irf.peak.pos);
    % Calculate area
    irf.area = sum(irf_sum);
    % Calculate baricenter
    x = 1:size(irf_sum,1);
    irf_sum2 = reshape(irf_sum,[size(irf_sum,1) 1]);
    irf.baric.pos = round(x*irf_sum2./irf.area);
    irf.baric.time = time(irf.baric.pos);
    % Shift temoral axis
    time = time - irf.baric.time;
    % Calculate second centered moment
    irf.variance = ((x.^2)*irf_sum./irf.area - irf.baric.pos.^2).*factor;
    %irf.sigma2 = sqrt(((x-irf.baric.pos).^2*irf.data./irf.area).*factor);
    if APPLY_TWIN == 1
        for w = 1:nwin
            y(w,:) = sum(irf.data(twin(w,1):twin(w,2),:),1);
            time(w) = mean(t(twin(w,1):twin(w,2)));
        end
        irf.data = y;
    end
    semilogy(sum(irf.data,2)),title('IRF')
end



% =========================================================================
%%                              END
% =========================================================================