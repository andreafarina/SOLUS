function EXP = createExp(DAT_STR, type_phantom)
%
% EXP = CreateMatlabExp(DAT_STR, type_phantom)
% Takes as input the path to a *.dat file and return a MAT variable
% suitable for operations in the SOLUS software
% DATA_STR is the pathname of a *.dat file
% type_phantom is a str defining from which phantom the experimental data should be created: 
% -'veallard'
% -'lardtend'
% -'silicon'
%  EXP is a Matlab variable with experimental data from the phantom
    old_pwd = pwd;
    TOMOFOLDER = from_dat2mat(DAT_STR, type_phantom);
    folder_data = DAT_STR;
    cd([TOMOFOLDER])
    %addpath((['./Tomo structs/' folder_data]));
    lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
    %lambda = [620 670 740 800 910 1020 1050 1090];

    for iw = 1:numel(lambda)
    EXP_file = strcat('Tomo_wave_',num2str(lambda(iw)));
    load(['Tomo_structs/',EXP_file,'.mat'])
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
    EXP.path.day = '202002';
    
    EXP.path.data_folder = [old_pwd,filesep, 'Tomo_structs'];%['../../results/201710/dynamic phantom'];

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
    cd(old_pwd)
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
    mkdir(fullfile(pwd,'EXP_structs',folder_data))
    str_file = [fullfile(pwd,'EXP_structs',folder_data),filesep,...
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
    str_file = [fullfile(pwd,'EXP_structs',folder_data),filesep,...
        'EXP_Tomo'];
    save(str_file,'EXP');

    return
end
% str_file = [fullfile(pwd,'EXP structs',folder_data),filesep,...
%     'EXP_','Tomo',output_suffix];
% save(str_file,'EXP');


function TomoFolder = from_dat2mat(STR, type_str)
%
% convert *.dat files into *.mat structures and returns the path to the
% folder where the structures are saved

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
    filepath = STR;
    %addpath(genpath(['./' filepath]));

    switch lower(type_str)
        %% LARD+HAM
        case 'lardham'
            xs = [19.5 6.5 -6.5 -19.5];
            ys = [10 -10];
            zs = 0;
            N_homo = [473:1:500];  % #homog. measurements
            N_hete = [501:1:528];  % #heterog. measurements
            irf_file_name = 'SOLs0100';
            lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
            T.opt.homo.abs = [0.0038079316	0.0013426612	0.0013876201 0.0100602733	0.0075543996	0.0039425139	0.007613757	0.004933052];
            T.opt.homo.sca = [1.228420919	1.189904619	0.913511025	0.8030371548 0.7696766762	0.7560782619	0.7022870595	0.6751796929];
            T.opt.hete.abs = [0.037321025	0.021343425	0.0104115095	0.01574798	0.023087425	0.05156646	0.02941183	0.019642555];
            T.opt.hete.sca = [0.3711594	0.31473145	0.2147255	0.20802585	0.22279435	0.2970374	0.2203812	0.18975735];
        %% VEAL+LARD
        case 'veallard'
            xs = [19.5 6.5 -6.5 -19.5];
            ys = [10 -10];
            zs = 0;
            N_homo = [672:1:699];  % #homog. measurements
            N_hete = [700:1:727];  % #heterog. measurements
            irf_file_name = 'SOLs0111';
            lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
            T.opt.homo.abs = [0.02791161  0.0177973093    0.009091243 0.01586907  0.0260947261    0.053704805 0.0297806604    0.01985166];
            T.opt.homo.sca = [0.9385278214	0.8507553357    0.5872533214    0.5547145536    0.5709639036    0.6376319821    0.5172092786    0.4790727179];
            T.opt.hete.abs = [0.0036857775	0.001488738	0.001309731	0.010038185	0.0073729725	0.0050007045	0.0080228255	0.0052807725];
            T.opt.hete.sca = [1.761344	1.639484	1.4247535	1.349348	1.2949265	1.200157	1.215836	1.2253315];
    %% LARD+TEND
        case 'lardtend'
            xs = [19.5 6.5 -6.5 -19.5];
            ys = [10 -10];
            zs = 0;
            N_homo = [614:1:641];  % #homog. measurements
            N_hete = [642:1:669];  % #heterog. measurements
            irf_file_name = 'SOLs0107';
            lambda   = [635, 670, 830, 915, 940, 980, 1030, 1065]; %lambda used
            T.opt.homo.abs = [0.0038005485	0.0034049469	0.0015822979	0.0104953881	0.0075209223	0.0052855228	0.0090089653	0.0054133905];
            T.opt.homo.sca = [1.7920386071	1.6857158	1.4810270357	1.4032054643	1.3324607857	1.2932485	1.2526574786	1.24276525];
            T.opt.hete.abs = [0.005657467	0.00265227	0.003164864	0.008215303	0.01391677	0.02775662	0.01905252	0.01384971];
            T.opt.hete.sca = [2.777692	2.854745	2.679655	2.560578	2.494992	2.154414	2.323298	2.578156];
    %% SIL_10_01_1cm
        case 'sil_10_01_1cm'
            xs = [20.0000   6.6667    -6.6667   -20.0000];
            ys = [14 -14];
            zs = 0;
            lambda = [620 670 740 800 910 1020 1050 1090];
            N_homo = [786:1:813];  % #homog. measurements
            N_hete = [842:1:869];  % #heterog. measurements
            irf_file_name = 'SOLs0116';
            T.opt.homo.abs = [0.00436839500000000 0.00399510520000000 0.00530931960000000 0.00449303100000000 0.0124849820000000 0.00865493720000000 0.00426235080000000 0.00704804460000000];
            T.opt.homo.sca = [1.12950240000000 0.996046560000000 0.920408340000000 0.858343020000000 0.656752800000000 0.537798400000000 0.524484660000000 0.472743960000000];
            T.opt.hete.abs = [0.0084826386	0.0079976316	0.0096794268	0.0086512542	0.016530138	0.012258804	0.0078621298	0.010433566];
            T.opt.hete.sca = [1.1072414	0.98465536	0.9394932	0.86275388	0.66062366	0.53074266	0.51240594	0.46274634];
    %% SIL_10_01_6cm
         case 'sil_10_01_6cm'
            xs = [20.0000   6.6667    -6.6667   -20.0000];
            ys = [14 -14];
            zs = 0;
            lambda = [620 670 740 800 910 1020 1050 1090];
            N_homo = [786:1:813];  % #homog. measurements
            N_hete = [898:1:925];  % #heterog. measurements
            irf_file_name = 'SOLs0116';
            T.opt.homo.abs = [0.00436839500000000 0.00399510520000000 0.00530931960000000 0.00449303100000000 0.0124849820000000 0.00865493720000000 0.00426235080000000 0.00704804460000000];
            T.opt.homo.sca = [1.12950240000000 0.996046560000000 0.920408340000000 0.858343020000000 0.656752800000000 0.537798400000000 0.524484660000000 0.472743960000000];
            T.opt.hete.abs = [0.0084826386	0.0079976316	0.0096794268	0.0086512542	0.016530138	0.012258804	0.0078621298	0.010433566];
            T.opt.hete.sca = [1.1072414	0.98465536	0.9394932	0.86275388	0.66062366	0.53074266	0.51240594	0.46274634];
    %% SIL_10_02_1cm
         case 'sil_10_02_1cm'
             xs = [20.0000   6.6667    -6.6667   -20.0000];
             ys = [14 -14];
             zs = 0;
             lambda = [620 670 740 800 910 1020 1050 1090];
             N_homo = [786:1:813];  % #homog. measurements
             N_hete = [814:1:841];  % #heterog. measurements
             irf_file_name = 'SOLs0116';
             T.opt.homo.abs = [0.00436839500000000 0.00399510520000000 0.00530931960000000 0.00449303100000000 0.0124849820000000 0.00865493720000000 0.00426235080000000 0.00704804460000000];
             T.opt.homo.sca = [1.12950240000000 0.996046560000000 0.920408340000000 0.858343020000000 0.656752800000000 0.537798400000000 0.524484660000000 0.472743960000000];
             T.opt.hete.abs = [0.01683116	0.01719449	0.019627084	0.017524408	0.024481872	0.019767674	0.015424918	0.017779734];
             T.opt.hete.sca = [0.98843248	0.9300621	0.91997866	0.7990178	0.5994491	0.47483964	0.4577962	0.4145961];
    %% SIL_10_02_6cm
        case 'sil_10_02_6cm'
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
    %% case silicon 200202
        case 'silicon_200202'
            xs = [20.0000   6.6667    -6.6667   -20.0000];
            ys = [14 -14];
            zs = 0;
            lambda = [635, 670, 830, 915, 940, 980, 1030, 1065];
            %:1:244];  % #homog. measurements
            %:1:244];  % #heterog. measurements
            irf_file_name = 'SOLs0211';
            str_phantom_homo = '10_0.2';
            str_phantom_hete = '10_0.4';
            switch lower(str_phantom_homo)
                case '10_0.1'
                    T.opt.homo.abs = [0.09118862 0.08252885 0.07457543 0.114975675 0.07321214 0.07269354 0.09152519 0.076381];
                    T.opt.homo.sca = [14.267035 13.06632 9.491283 7.8722685 7.707295 7.060895 6.352237 5.989931];
                    N_homo = [246];% #homog. measurements
                    
                    switch lower(str_phantom_hete)
                        case '10_0.05'
                            T.opt.hete.abs = [0.05246881 0.04837734 0.03696393 0.080474885 0.03518626 0.03308078 0.05274401 0.038200645];
                            T.opt.hete.sca = [11.46904 10.55009 7.656952 6.150009 6.002163 5.590301 4.857868 4.654679];
                            N_hete = [247];% #heterog. measurements
                        case '10_0.1'
                            T.opt.hete.abs = [0.094249395 0.09108484 0.07822272 0.12349727 0.07722901 0.07274856 0.09558736 0.081897605];
                            T.opt.hete.sca = [11.27842 10.56916 7.576539 6.138904 5.994108 5.46071 4.877412 4.694373];               
                            N_hete = [248];% #heterog. measurements
                        case '10_0.2'
                            T.opt.hete.abs = [0.18834215 0.1808696 0.1594672 0.2140688 0.1618171 0.1609721 0.1816612 0.1681408];
                            T.opt.hete.sca = [11.196185 10.30879 7.183976 6.0505265 5.763012 5.358386 4.799074 4.5417635];
                            N_hete = [249];% #heterog. measurements
                        case '10_0.4'
                            T.opt.hete.abs = [0.3244506 0.3090453 0.2850104 0.3521497 0.2937037 0.2913416 0.3152952 0.30309455];
                            T.opt.hete.sca = [9.540982 8.508501 5.992858 5.245241 4.992857 4.637014 4.149682 3.973167];
                            N_hete = [250];% #heterog. measurements
                        case '10_0.6'
                            T.opt.hete.abs = [0.46561525 0.4456669 0.422469 0.5012814 0.4393349 0.4392696 0.4642975 0.4596188];
                            T.opt.hete.sca = [9.5278705 8.628135 6.08301 5.301621 5.010855 4.710412 4.163082 4.0690635];
                            N_hete = [251];% #heterog. measurements
                        case '15_0.1'
                            T.opt.hete.abs = [0.080407055 0.07623264 0.06779367 0.110158655 0.06512177 0.06030613 0.07829765 0.066305675];
                            T.opt.hete.sca = [16.99245 15.68388 11.63303 9.3891785 9.218989 8.425847 7.417703 7.211785];
                            N_hete = [252];% #heterog. measurements                   
                        case '05_0.1'
                            T.opt.hete.abs = [0.087325025 0.08533506 0.07396581 0.1459342 0.08532022 0.08495646 0.1181968 0.10573535];
                            T.opt.hete.sca = [4.766733 4.387157 3.001634 2.639426 2.415703 2.165157 1.945762 1.762059];
                            N_hete = [253];% #heterog. measurements
                    end
                    
                case '10_0.2'
                    T.opt.homo.abs = [0.176269 0.1642548 0.1523256 0.1921913 0.148854 0.1484415 0.1653411 0.1511315];
                    T.opt.homo.sca = [13.239025 12.28287 8.910008 7.3108825 7.121039 6.537811 5.882884 5.576854];
                    N_homo = [254];% #homog. measurements
                    
                    switch lower(str_phantom_hete)
                        case '10_0.2'
                            T.opt.hete.abs = [0.18834215 0.1808696 0.1594672 0.2140688 0.1618171 0.1609721 0.1816612 0.1681408];
                            T.opt.hete.sca = [11.196185 10.30879 7.183976 6.0505265 5.763012 5.358386 4.799074 4.5417635];
                            N_hete = [255];% #heterog. measurements
                        case '10_0.4'
                            T.opt.hete.abs = [0.3244506 0.3090453 0.2850104 0.3521497 0.2937037 0.2913416 0.3152952 0.30309455];
                            T.opt.hete.sca = [9.540982 8.508501 5.992858 5.245241 4.992857 4.637014 4.149682 3.973167];
                            N_hete = [256];% #heterog. measurements

                    end
            end
 
                    
    end
    
    
    meas_homo = [];
    meas_hete = [];
    B = [];
    C = [];
    TomoFolder = [ filepath filesep 'Tomo_structs'];
    old_pwd = pwd;
    cd(filepath)
    % Homog. measurements
    for i = 1:length(N_homo)
        B    = DatRead3(strcat('SOLm1',num2str(N_homo(i))),4096,5,8);
        if sum_rep == 0
            meas_homo(:,i,:) = B(:,1,:);  %takes just the first repetion.
        elseif sum_rep == 1
            meas_homo(:,i,:) = squeeze(B); %sums the 5 repetitions
        else
            meas_homo(:,i,:) = squeeze(sum(B,2)); %sums the 5 repetitions
        end
    end 

    % Heterog. measurements
    for i = 1:length(N_hete)
        C    = DatRead3(strcat('SOLm1',num2str(N_hete(i))),4096,5,8);
        if sum_rep == 0
            meas_hete(:,i,:) = C(:,1,:);  %takes just the first repetion.
        elseif sum_rep == 1
            meas_hete(:,i,:) = squeeze(C); %sums the 5 repetitions
        else
            meas_hete(:,i,:) = squeeze(sum(C,2)); %sums the 5 repetitions
        end
    end

    meas_homo = permute(meas_homo,[3 2 1]);  %rearrange of the dimensions
    meas_hete = permute(meas_hete,[3 2 1]);

    %Irf file
    irf = DatRead3(irf_file_name,4096,5,8);  
        if sum_rep == 0 
            irf = squeeze(permute(irf(:,1,:),[3 2 1])); %rearrange of the dimensions
        elseif sum_rep == 1
            irf = squeeze(permute(irf(:,:,:),[3 2 1]));
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
        back_irf  = mean(irf(back_start:back_stop,:));

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
    dmask = dmask - diag(diag(ones(size(SOURCE_POS,1),size(DET_POS,1)))); % we didn't perform rho = 0 and 1->2 , 2->1

    %% rearrange measurements
    meas_homo0 = zeros(Nchannel,sum(dmask(:)) ,size(irf,2));
    meas_hete0 = zeros(Nchannel,sum(dmask(:)),size(irf,2));
    for lam= 1:size(irf,2)
        meas_homo0(:,:,lam) = squeeze(meas_homo(:, 1, lam:size(irf,2):end ));
        meas_hete0(:,:,lam) = squeeze(meas_hete(:, 1, lam:size(irf,2):end));
    end
    count = 1;
%     srdet = cat(2,[0,1:7]',...
%                     [8,0,9:14]',...
%                     [15:16, 0, 17:21]',...
%                     [22:24,0,25:28]',...
%                     [29:32,0,33:35]',...
%                     [36:40,0,41:42]',...
%                     [43:48,0,49]',...
%                     [50:56,0]');
%     srdetmeas = cat(2,[21,8:14]',...
%                       [27:28,15:20]',...
%                       [33:35,22:26]',...
%                       [39:42,29:32]',...
%                       [45:49,36:38]',...
%                       [51:56,43:44]',...
%                   [1:7,50]');
    srdetmeas = cat(2,[0,1:7]',...
                      [14,0,8:13]',...
                      [20:21,0,15:19]',...
                      [26:28,0,22:25]',...
                      [32:35,0,29:31]',...
                      [38:42,0,36:37]',...
                      [44:49,0,43]',...
                  [50:56,0]');
    %srdet = srdetmeas([2,3,4,5,1,8,7,6],[2,3,4,5,1,8,7,6]);        
    srdet = srdetmeas([5, 4, 3, 2, 6, 7, 8, 1],[5, 4, 3, 2, 6, 7, 8, 1])';          
    %srdet = srdetmeas([5,1,2,3,4,8,7,6],[5,1,2,3,4,8,7,6]); 


%     srdetmeas = cat(2,[8:14]',...
%                   [21,15:20]',...
%                   [27:28,22:26]',...
%                   [33:35,29:32]',...
%                   [39:42,36:38]',...
%                   [45:49,43:44]',...
%               [51:56,50]',...
%               [1:7]')';
              
    meas_homo00 = zeros(size(meas_homo0));
    meas_hete00 = meas_homo00;
    for lam= 1:size(irf,2)
        dc = 0;
        for d = 1:64
            if srdetmeas(d)>0 %& srdetmeas(d)>0
            dc = dc +1;
            disp(srdet(d))
            meas_homo00(:,(dc),lam) = meas_homo0(:,srdet(d), lam);
            meas_hete00(:,(dc),lam) = meas_hete0(:,srdet(d), lam);
            end
        end
    end

    meas_homo = meas_homo00;
    meas_hete = meas_hete00;
    clearvars meas_homo0 meas_hete0
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
    cd(old_pwd)
    mkdir(TomoFolder);
    cd(TomoFolder);
    T.path.data_folder = TomoFolder;

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
 cd(old_pwd)   
end
%cd(fileparts(mfilename('fullpath')))
