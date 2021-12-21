function [dt,nstep] = genEXP_fromSOLUS(savename,priorname, SOLUSfile,flag_segm,nom_vals)
% [dt,nstep]=genEXP_fromSOLUS(savename,priorname, SOLUSfile,flag_segm,nom_vals)
% It generates an EXP file and a Prior with the name specifiend in savename
% and priorname respectively starting form the file identified by SOLUSfile
% flag_segm: flag for segmentation:
% - 0 : use value from structure in SOLUSfile
% - 1 : segment and extract volume
% - 2 : load segmentation from file in nom_vals{5}
% - 3 : segmentation with snake fitting and extract volume
% - []: No value provided
% nom_vals: array of cells in the form {homo,hete,{log1,log2},{spe1,spe2},nom_vals{5}}
% where 
% --homo is the name of the field in the stucture of SOLUSfile assigned to be 
% the reference measurement in reconstruction, default is 'controlateral'
% -- hete is the name of the field in the stucture of SOLUSfile assigned to be
% the heterogeneous measurement in reconstruction, default is 'main'
% -- log1 and log2, when available, specify respectively the .xlsx file and the sheet therein
% with the logbook on the experiment
% -- spe1 and spe2, when available, specify respectively the .xlsx file and the sheet therein
% with the spectra of the media in the experiment
% -- nom_vals{5} extra parameter for the segmentation options above



homo = 'controlateral';
hete = 'main';

if exist('nom_vals','var') 
   
    if ~isempty(nom_vals)
        homo = nom_vals{1};%'controlateral';
        hete = nom_vals{2};
    end
    if numel(nom_vals) >2 % read ground truth from log
        logbook = nom_vals{3};
        spectrabook = nom_vals{4};
    else
        logbook = [];
        spectrabook = [];
    end
end


%% start
% load files
load(SOLUSfile);
orig = Measure;
ref.Measure = Measure.(homo).DOT;%load(refloadname);
irf.Measure = Measure.IRF.DOT;
goodq.Measure = Measure.(hete).DOT;
Measure = Measure.(hete).DOT;



%irf.Measure.Data = cat(4, irf.Measure.Data, 0*irf.Measure.Data);
goodq.Measure.Mask = goodq.Measure.Mask.*(sum(Measure.Data,4)~=0);
Measure.Mask = goodq.Measure.Mask.*(sum(Measure.Data,4)~=0);


%disp('loaded')

vec_perm = [3,1,2];
VEC_perm = [4,3,1,2];

strdsl1 = '(:,d,s,l)';
strdsl2 = '(:,d,s,l)';

permOptodesD = [5,6,7,8,4,3,2,1];
permOptodesS = permOptodesD;

%% data
EXP.data.ref = permute(ref.Measure.Data,VEC_perm);
EXP.data.ref = EXP.data.ref(:, permOptodesD, permOptodesS,:); 
EXP.data.spc = permute(Measure.Data,VEC_perm);
EXP.data.spc = EXP.data.spc(:, permOptodesD, permOptodesS,:);
irfmat = permute(irf.Measure.Data,VEC_perm);
irfmat = irfmat(:, permOptodesD, permOptodesS,:);
irfmat(isnan(irfmat) ) = 0;


%% oversample and downsample
timeD = [];
bins = Measure.BinWidth_ps(permOptodesD);
for i = 1:numel(bins)
    timeD(:,i) = 0: bins(i)  : bins(i)*(size(irfmat,1)-1);
end

time = 0:mean(bins(:)):mean(bins(:))*(size(irfmat,1)-1);
dt = mean(bins(:));

str= {'irfmat', 'EXP.data.spc','EXP.data.ref'};
for icell = 1:numel(str) 
    tmpV = eval([str{icell},';']);%#ok
    for d =1:8
        timeOVER = 0:1:timeD(end,i);
        for s = 1:8
            for l=1:8
               V = squeeze(eval(['tmpV',strdsl1,';'])); 
               yb = interp1(timeD(:,d), double(V), timeOVER);
               y = interp1(timeOVER, yb, time);
               y(isnan(y)) = 0;
               eval([str{icell},strdsl2,'= y;']);
               if icell ==1
                  eval(['tmpirfovertime',strdsl2,'= timeOVER;']);
                  eval(['tmpirfmat',strdsl2,'= yb;']);
               end
            end
        end
    end
end


%% path
[path,name] = fileparts(SOLUSfile);
EXP.path.day =datestr(now,'YYYYmmDD__hhMM');
EXP.path.data_folder = path;
EXP.path.file_name = name;
EXP.path.irf_file = SOLUSfile;

%% spc
EXP.spc.gain = [];
EXP.spc.n_chan = size(irfmat,1);
EXP.factor = [];
%% time

EXP.time.roi = [];
EXP.time.axis = time;
%% bkg
EXP.bkg.ch_start = [];
EXP.bkg.ch_end = [];
%% grid
EXP.grid.xs = [19.5000    6.5   -6.5  -19.5000];%[20.0000    6.6667   -6.6667  -20.0000];
EXP.grid.ys = [17.5 -17.5];%[14 -14]
EXP.grid.zs = [ 0 ];
EXP.grid.xd = [19.5000    6.5   -6.5  -19.5000];
EXP.grid.yd = [10.8 -10.8];
EXP.grid.zd = [ 0 ];
[xxs,yys,zzs] = ndgrid(EXP.grid.xs,EXP.grid.ys,EXP.grid.zs);
[xxd,yyd,zzd] = ndgrid(EXP.grid.xd,EXP.grid.yd,EXP.grid.zd);
EXP.grid.Source.Pos = [xxs(:),yys(:),zzs(:)];
EXP.grid.Detector.Pos = [xxd(:),yyd(:),zzd(:)];
Nd = numel(xxs);
Ns = numel(xxs);
EXP.grid.dmask = permute(Measure.Mask, vec_perm).*permute(irf.Measure.Mask, vec_perm);%.* logical(ones(Nd,Ns) - diag(ones(Nd,Ns)));
EXP.grid.dmask = permute(Measure.Mask, vec_perm).*permute(irf.Measure.Mask, vec_perm);
EXP.grid.dmask = EXP.grid.dmask(permOptodesD,permOptodesS,:);
%EXP.grid.dmask = permute(EXP.grid.dmask,[1,3,2,4]);
%% lambda
EXP.lambda = Measure.Wavelengths_nm;

%% optp
if isempty(logbook) || isempty(spectrabook)
    EXP.optp.homo.abs = repmat(0.01,[8,1]); %extractfromTab(tab,0,BULKA(i_control),BULKS(i_control));
    EXP.optp.homo.sca = repmat(1,[8,1]);
    EXP.optp.hete.abs = repmat(0.01,[8,1]);
    EXP.optp.hete.sca = repmat(1,[8,1]);%extractfromTab(tab,1,INCLA(i_control),INCLS(i_control));
else
    [EXP.optp.homo.abs, EXP.optp.homo.sca,EXP.optp.hete.abs, EXP.optp.hete.sca] = getTruthfromLog(SOLUSfile,logbook,spectrabook);
end
%% shift irf
%
Tshift3 = irf.Measure.Tshift_ps;
dist =  [0.67;1.4625;2.6849;2.83;3.1143;3.843;3.9571;4.8186]; %cm

%% irf
distmat = round(sqrt((EXP.grid.Source.Pos(:,1) - EXP.grid.Detector.Pos(:,1)').^(2) +...
    (EXP.grid.Source.Pos(:,2) - EXP.grid.Detector.Pos(:,2)').^(2)+...
    (EXP.grid.Source.Pos(:,3) - EXP.grid.Detector.Pos(:,3)').^(2)),1); 
distmat = repmat(distmat,[1,1,8]);
dists = sort(unique(round(distmat, 1)), 'ascend');

tmpirflin = reshape(tmpirfmat, [size(tmpirfmat,1),numel(tmpirfmat)/size(tmpirfmat,1)]);
%tmpirflin_ = tmpirflin;
irflintime = reshape(tmpirfovertime, [size(tmpirfmat,1),numel(tmpirfmat)/size(tmpirfmat,1)]);
a = 0;
for i =1:numel(dists)
    IDX = logical((distmat==dists(i))); 
    y = circshift(tmpirflin(:,IDX(:)'), Tshift3(i),1);
    y(end+Tshift3(i):end,:) = 0;
    tmpirflin(:,IDX(:)') =  y;
end
irflin = [];
for i = 1:size(tmpirflin,2) 
   irflin(:,i) = interp1(irflintime(:,i),tmpirflin(:,i), time); 
end

%irfmat = reshape(irflin, size(irfmat));
EXP.irf.data = irflin;
EXP.irf.area = squeeze(sum(irflin,1));
EXP.irf.variance = var(irflin,1);

for i = 1:size(irflin,2)
    tmp = irflin(:,i);  
    
    [~,I] = max(time .* tmp);
    EXP.irf.baric(i).pos = I;
    EXP.irf.baric(i).time = time(I);
    
    [V,I] = max(tmp(:));
    EXP.irf.peak(i).value = V; 
    EXP.irf.peak(i).pos =  I;
    EXP.irf.peak(i).time = time(I);   
end
%% adjust dimension
EXP.data.spc = double(reshape(EXP.data.spc, [size(EXP.data.spc,1), numel(EXP.data.spc)/size(EXP.data.spc,1)]));
EXP.data.ref = double(reshape(EXP.data.ref, [size(EXP.data.ref,1), numel(EXP.data.ref)/size(EXP.data.ref,1)]));

EXP.geom.homo = homo;
EXP.geom.hete = hete;

%% saving EXP
save(savename,'EXP');
%% handle prior flag_segm ==0 
if ~isempty(flag_segm)
    switch flag_segm
        case 0 % take from file
            Mask3D = orig.Segmentation.Mask3D;
            delta = str2double(orig.controlateral.US.Parameters.width_mm_per_px);
        case 1 % manual segmentation
            disp('-- segmentation...')
            Mask3D =  dicom_to_DT(orig.main.US.Image);
            delta = str2double(orig.main.US.Parameters.height_mm_per_px);
        case 2 %load
            load(nom_vals{5})
        case 3 % snake fitting
            disp('-- segmentation...')
            Mask3D =  dicom_to_DT(orig.main.US.Image,1);
            delta = str2double(orig.main.US.Parameters.heightmm_per_px);        
    end
    save(priorname,'Mask3D','delta', '-v7.3')
else
    if ~(exist([priorname],'file')>0) && ~(exist([priorname,'.mat'],'file')>0)
        warning('Variable ''priorname'' does not exist in folder, but no prior was chosen')
    else 
        disp('\t Using prior in folder')
    end
end
dt = EXP.time.axis(2)-EXP.time.axis(1);
nstep = size(EXP.data.spc,1);
return

end