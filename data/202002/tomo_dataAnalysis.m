%% Create a file.mat with the data used for the tomographic reconstruction
%  version @ 11-June-2018
%% Content: 
%           heterogeneous => data(n_time,n_sd,n_lambda);
%           homogeneous   => ref(n_time,n_sd,n_lambda);
%           irf           => irf(n_time,n_lambda);
%           stand. dev.   => sd(n_time,n_sd,n_lambda);
%           Source Position => SOURCE_POS;
%           Det. Position   => DETECTOR_POS;
%           Logic Mask      => dmask(n_detector, n_source);
%=========================================================================
%% Clean workspace
InitScript

%=========================================================================
%% TCSPC board parameters
gain     = 6;
Nchannel = 4096;
factor   = 50000/gain/Nchannel;
phot_goal = 400000;            %Total photons per curve (goal)
repetitions = 5;               %# of repetitions for each lambda
roi = [500 3500];
%% Background removal
background = 1;     %background = 1 => remove the backgorund
                    %background = 0 => no removal
%=========================================================================
%% Repetitions sum

sum_rep = 1;        %sum_rep = 1 => the repetitions are summed
                    %sum_rep = 0 => takes just the first repetition
%=========================================================================
%% Read .DAT files
% Measurement file
filepath = 'Dati phantom silicone/Data_hete_10_0.2_6cm';
addpath(genpath(['./' filepath]));
%% LARD+HAM
% xs = [19.5 6.5 -6.5 -19.5];
% ys = [10 -10];
% zs = 0;
% N_homo = [473:1:500];  % #homog. measurements
% N_hete = [501:1:528];  % #heterog. measurements
% irf_file_name = 'SOLs0100';
% lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
% T.opt.homo.abs = [0.0038079316	0.0013426612	0.0013876201 0.0100602733	0.0075543996	0.0039425139	0.007613757	0.004933052];
% T.opt.homo.sca = [1.228420919	1.189904619	0.913511025	0.8030371548 0.7696766762	0.7560782619	0.7022870595	0.6751796929];
% T.opt.hete.abs = [0.037321025	0.021343425	0.0104115095	0.01574798	0.023087425	0.05156646	0.02941183	0.019642555];
% T.opt.hete.sca = [0.3711594	0.31473145	0.2147255	0.20802585	0.22279435	0.2970374	0.2203812	0.18975735];
%% VEAL+LARD
% xs = [19.5 6.5 -6.5 -19.5];
% ys = [10 -10];
% zs = 0;
% N_homo = [672:1:699];  % #homog. measurements
% N_hete = [700:1:727];  % #heterog. measurements
% irf_file_name = 'SOLs0111';
% lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
% T.opt.homo.abs = [0.02791161  0.0177973093    0.009091243 0.01586907  0.0260947261    0.053704805 0.0297806604    0.01985166];
% T.opt.homo.sca = [0.9385278214	0.8507553357    0.5872533214    0.5547145536    0.5709639036    0.6376319821    0.5172092786    0.4790727179];
% T.opt.hete.abs = [0.0036857775	0.001488738	0.001309731	0.010038185	0.0073729725	0.0050007045	0.0080228255	0.0052807725];
% T.opt.hete.sca = [1.761344	1.639484	1.4247535	1.349348	1.2949265	1.200157	1.215836	1.2253315];
%% LARD+TEND
% xs = [19.5 6.5 -6.5 -19.5];
% ys = [10 -10];
% zs = 0;
% N_homo = [614:1:641];  % #homog. measurements
% N_hete = [642:1:669];  % #heterog. measurements
% irf_file_name = 'SOLs0107';
% lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
% T.opt.homo.abs = [0.0038005485	0.0034049469	0.0015822979	0.0104953881	0.0075209223	0.0052855228	0.0090089653	0.0054133905];
% T.opt.homo.sca = [1.7920386071	1.6857158	1.4810270357	1.4032054643	1.3324607857	1.2932485	1.2526574786	1.24276525];
% T.opt.hete.abs = [0.005657467	0.00265227	0.003164864	0.008215303	0.01391677	0.02775662	0.01905252	0.01384971];
% T.opt.hete.sca = [2.777692	2.854745	2.679655	2.560578	2.494992	2.154414	2.323298	2.578156];
%% SIL_10_01_1cm
% xs = [20.0000   6.6667    -6.6667   -20.0000];
% ys = [14 -14];
% zs = 0;
% lambda = [620 670 740 800 910 1020 1050 1090];
% N_homo = [786:1:813];  % #homog. measurements
% N_hete = [842:1:869];  % #heterog. measurements
% irf_file_name = 'SOLs0116';
% T.opt.homo.abs = [0.00436839500000000 0.00399510520000000 0.00530931960000000 0.00449303100000000 0.0124849820000000 0.00865493720000000 0.00426235080000000 0.00704804460000000];
% T.opt.homo.sca = [1.12950240000000 0.996046560000000 0.920408340000000 0.858343020000000 0.656752800000000 0.537798400000000 0.524484660000000 0.472743960000000];
% T.opt.hete.abs = [0.0084826386	0.0079976316	0.0096794268	0.0086512542	0.016530138	0.012258804	0.0078621298	0.010433566];
% T.opt.hete.sca = [1.1072414	0.98465536	0.9394932	0.86275388	0.66062366	0.53074266	0.51240594	0.46274634];
%% SIL_10_01_6cm
% xs = [20.0000   6.6667    -6.6667   -20.0000];
% ys = [14 -14];
% zs = 0;
% lambda = [620 670 740 800 910 1020 1050 1090];
% N_homo = [786:1:813];  % #homog. measurements
% N_hete = [898:1:925];  % #heterog. measurements
% irf_file_name = 'SOLs0116';
% T.opt.homo.abs = [0.00436839500000000 0.00399510520000000 0.00530931960000000 0.00449303100000000 0.0124849820000000 0.00865493720000000 0.00426235080000000 0.00704804460000000];
% T.opt.homo.sca = [1.12950240000000 0.996046560000000 0.920408340000000 0.858343020000000 0.656752800000000 0.537798400000000 0.524484660000000 0.472743960000000];
% T.opt.hete.abs = [0.0084826386	0.0079976316	0.0096794268	0.0086512542	0.016530138	0.012258804	0.0078621298	0.010433566];
% T.opt.hete.sca = [1.1072414	0.98465536	0.9394932	0.86275388	0.66062366	0.53074266	0.51240594	0.46274634];
%% SIL_10_02_1cm
% xs = [20.0000   6.6667    -6.6667   -20.0000];
% ys = [14 -14];
% zs = 0;
% lambda = [620 670 740 800 910 1020 1050 1090];
% N_homo = [786:1:813];  % #homog. measurements
% N_hete = [814:1:841];  % #heterog. measurements
% irf_file_name = 'SOLs0116';
% T.opt.homo.abs = [0.00436839500000000 0.00399510520000000 0.00530931960000000 0.00449303100000000 0.0124849820000000 0.00865493720000000 0.00426235080000000 0.00704804460000000];
% T.opt.homo.sca = [1.12950240000000 0.996046560000000 0.920408340000000 0.858343020000000 0.656752800000000 0.537798400000000 0.524484660000000 0.472743960000000];
% T.opt.hete.abs = [0.01683116	0.01719449	0.019627084	0.017524408	0.024481872	0.019767674	0.015424918	0.017779734];
% T.opt.hete.sca = [0.98843248	0.9300621	0.91997866	0.7990178	0.5994491	0.47483964	0.4577962	0.4145961];
%% SIL_10_02_6cm
xs = [20.0000   6.6667    -6.6667   -20.0000];
ys = [14 -14];
zs = 0;
lambda = [620 670 740 800 910 1020 1050 1090];
N_homo = [786:1:813];  % #homog. measurements
N_hete = [870:1:897];  % #heterog. measurements
irf_file_name = 'SOLs0116';
T.opt.homo.abs = [0.00436839500000000 0.00399510520000000 0.00530931960000000 0.00449303100000000 0.0124849820000000 0.00865493720000000 0.00426235080000000 0.00704804460000000];
T.opt.homo.sca = [1.12950240000000 0.996046560000000 0.920408340000000 0.858343020000000 0.656752800000000 0.537798400000000 0.524484660000000 0.472743960000000];
T.opt.hete.abs = [0.01683116	0.01719449	0.019627084	0.017524408	0.024481872	0.019767674	0.015424918	0.017779734];
T.opt.hete.sca = [0.98843248	0.9300621	0.91997866	0.7990178	0.5994491	0.47483964	0.4577962	0.4145961];

meas_homo = [];
meas_hete = [];
B = [];
C = [];
TomoFolder = ['Tomo structs' filesep filepath];

% Homog. measurements
for i = 1:length(N_homo)
    B    = DatRead2(strcat('SOLm0',num2str(N_homo(i))),4096,5,8);
    if sum_rep == 0
        meas_homo(:,i,:) = B(:,1,:);  %takes just the first repetion.
    else
        meas_homo(:,i,:) = squeeze(sum(B,2)); %sums the 5 repetitions
    end
end 

% Heterog. measurements
for i = 1:length(N_hete)
    C    = DatRead2(strcat('SOLm0',num2str(N_hete(i))),4096,5,8);
    if sum_rep == 0
        meas_hete(:,i,:) = C(:,1,:);  %takes just the first repetion.
    else
        meas_hete(:,i,:) = squeeze(sum(C,2)); %sums the 5 repetitions
    end
end

meas_homo = permute(meas_homo,[3 2 1]);  %rearrange of the dimensions
meas_hete = permute(meas_hete,[3 2 1]);

%Irf file
irf = (irf_file_name,4096,5,8);  
    if sum_rep == 0 
        irf = squeeze(permute(irf(:,1,:),[3 2 1])); %rearrange of the dimensions
    else
        irf = squeeze(permute(sum(irf,2),[3 2 1]));
    end
% for i = 1: size(irf,2)
%     semilogy(irf(:,i));
%     hold on, grid on
% end
%=========================================================================
%% Subtraction of the background

if background == 0
    
else
    back_start = 450;  %first channel for the background
    back_stop  = 500;  %last channel for the background
    
    back_ref  = mean(meas_homo(back_start:back_stop,:,:));
    back_data = mean(meas_hete(back_start:back_stop,:,:));
    back_irf  = mean(irf(back_start:back_stop,:,:));
    
    meas_homo  = meas_homo  - back_ref;
    meas_homo(meas_homo < 0) = 0;
    meas_hete = meas_hete - back_data;
    meas_hete(meas_hete < 0) = 0;
    irf  = irf  - back_irf;
    irf(irf < 0) = 0;
end

% for i = 1: size(data,2)
%     semilogy(ref(:,i,1));
%     hold on, grid on
% end
%=========================================================================

%% Calculation of the standard deviation
sd = sqrt(meas_hete);
%=========================================================================

%% Create the S-D pos.
% SOURCE_POS = [xs1 ys1 zs1; xs2 ys2 zs2; ... ; xs8 ys8 zs8];
% SOURCE_POS = [xd1 yd1 zd1; xd2 yd2 zd2; ... ; xd8 yd8 zd8];

SOURCE_POS = [1.95 1 0; 0.65 1 0; -0.65 1 0; -1.95 1 0;...
              1.95 -1 0; 0.65 -1 0; -0.65 -1 0; -1.95 -1 0].*10;
DET_POS    = SOURCE_POS;
%=========================================================================

%% Create the logic mask
% dmask = [s1-d1 s2-d1 ... sn-dn; s1-d2 s2-d2 ... sn-d2 ...]
dmask = ones(size(SOURCE_POS,1),size(DET_POS,1));
dmask = dmask - triu(dmask); % we didn't perform rho = 0 and 1->2 , 2->1
CalcolaContrasto
%imagesc(dmask)
%=========================================================================
T.dmask = dmask;
T.sourcepos = SOURCE_POS;
T.detpos = DET_POS;
T.xs = xs;
T.ys = ys;
T.zs = zs;
T.gain = gain;
T.Nchannel = Nchannel;
T.factor = factor;
T.lambda = lambda;
T.roi = roi;
T.bkg.ch_start = back_start;
T.bkg.ch_stop = back_stop;
mkdir(TomoFolder);
cd(TomoFolder);
T.path.data_folder = pwd;

for iw=1:numel(T.lambda)
    T.irf.data = irf(:,iw);
    T.irf.t0 = 0;
    T.meas_hete = zeros(T.Nchannel,numel(T.dmask));
    T.meas_homo = zeros(T.Nchannel,numel(T.dmask));
    idm = 1;
    for im = 1:numel(T.dmask(:))
       if T.dmask(im)==1
            T.meas_hete(:,im)=meas_hete(:,idm,iw);
            T.meas_homo(:,im)=meas_homo(:,idm,iw);
            idm = idm +1;
       end
    end
    T.meas_hete = reshape(T.meas_hete,T.Nchannel,8,8);
    T.meas_homo = reshape(T.meas_homo,T.Nchannel,8,8);
    T.path.file_name = ['Tomo_wave_' num2str(T.lambda(iw))];
    save(['Tomo_wave_' num2str(T.lambda(iw))],'T');
end
cd(fileparts(mfilename('fullpath')))