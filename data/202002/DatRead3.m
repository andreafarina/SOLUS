function [ Data, varargout ] = DatRead3(FileName,varargin)
%DatRead3('FileName')
%Can be as input a selection of the following parameters
%DatRead3(...,'loop4',numloop4,'loop5',numloop5,'datatype','uint32','compilesub',true/false,'forcereading',true/false)
%[Data,Header,EasyReadableHead,SubHeaders,EasyReadableSubHead,UnSqueezedHeader]=DatRead3(...)

NumArgOut = nargout-1;
NumArgin = nargin-1;
ForceReading = false;
isCompileSubHeader = false;
HeadLen=764;
datatype = 'ushort';
nBoard = 1; nDet=1;
loop5=1;loop4=1;
ismandatoryarg = zeros(8,1);
m_loop5 = 1;m_loop4 = 2;m_loop3 = 3;m_loop2 = 4;m_loop1 = 5;
m_nsource = 6;m_ndet = 7;m_nboard = 8;m_nbin = 9;
isdatatypeargin=0;

FilePath = [FileName '.DAT'];
if isempty(fileparts(FileName))
    FilePath = fullfile(pwd,[FileName,'.DAT']);
end
fid=fopen(FilePath,'rb');
if fid<0, errordlg('File not found'); Data = []; return; end

for iN = 1:NumArgin
    if strcmpi(varargin{iN},'loop5')
        loop5 = varargin{iN+1};
        ismandatoryarg(m_loop5)=1;
    end
    if strcmpi(varargin{iN},'loop4')
        loop4 = varargin{iN+1};
        ismandatoryarg(m_loop4)=1;
    end
    if strcmpi(varargin{iN},'loop3')
        CompiledHeader.LoopNum(3) = varargin{iN+1};
        ismandatoryarg(m_loop3)=1;
    end
    if strcmpi(varargin{iN},'loop2')
        CompiledHeader.LoopNum(2) = varargin{iN+1};
        ismandatoryarg(m_loop2)=1;
    end
    if strcmpi(varargin{iN},'loop1')
        CompiledHeader.LoopNum(1) = varargin{iN+1};
        ismandatoryarg(m_loop1)=1;
    end
    if strcmpi(varargin{iN},'datatype')
        datatype = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'compilesub')
        isCompileSubHeader = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'forcereading')
        if(ischar(varargin{iN+1})||isstring(varargin{iN+1}))
            varargin{iN+1}=string2boolean(varargin{iN+1});
        end
        ForceReading = logical(varargin{iN+1});
    end
    if strcmpi(varargin{iN},'nSource')
        nDet = varargin{iN+1};
        ismandatoryarg(m_nsource)=1;
    end
    if strcmpi(varargin{iN},'nDet')
        nDet = varargin{iN+1};
        ismandatoryarg(m_ndet)=1;
    end
    if strcmpi(varargin{iN},'nBoard')
        nBoard = varargin{iN+1};
        ismandatoryarg(m_nboard)=1;
    end
    if strcmpi(varargin{iN},'nBin')
        nBin = varargin{iN+1};
        ismandatoryarg(m_nbin)=1;
    end
end

Head=fread(fid,HeadLen,'uint8');
if (numel(categories(categorical(Head)))==1||sum(Head)==0)&&sum(ismandatoryarg)~=numel(ismandatoryarg)
    errordlg('Please insert all loop values and nDet, nSource, nBoard, nBin'); Data = [];
    return;
end

CompiledHeader = FillHeader(Head);
SubLen=CompiledHeader.SizeSubHeader;
if SubLen == 0
    SkipSub = true;
    if(~all(ismandatory([m_nboard m_ndet m_nsource])))
        errordlg('Please insert nSource, nDet, nBoard'); Data = [];
    end
    return;
else
    SkipSub = false;
end
nBin = CompiledHeader.McaChannNum;
datasize = CompiledHeader.SizeData;
CompiledHeader.LoopNum(4) = loop4; CompiledHeader.LoopNum(5) = loop5;
CompiledHeader.LoopFirst(4) = 0; CompiledHeader.LoopFirst(5) = 0;
CompiledHeader.LoopLast(4) = 0; CompiledHeader.LoopLast(5) = 0;
CompiledHeader.LoopDelta(4) = 1; CompiledHeader.LoopDelta(5) = 1;
CompiledHeader.LoopHome(4) = 0; CompiledHeader.LoopHome(5) = 0;
CompiledHeader.LoopNum(1:3)=flip(CompiledHeader.LoopNum(1:3));
CompiledHeader.LoopFirst(1:3)=flip(CompiledHeader.LoopFirst(1:3));
CompiledHeader.LoopLast(1:3)=flip(CompiledHeader.LoopLast(1:3));
CompiledHeader.LoopDelta(1:3)=flip(CompiledHeader.LoopDelta(1:3));
CompiledHeader.LoopHome(1:3)=flip(CompiledHeader.LoopHome(1:3));

for iN = 1:NumArgin
    if strcmpi(varargin{iN},'loop5')
        loop5 = varargin{iN+1};
        ismandatoryarg(m_loop5)=1;
    end
    if strcmpi(varargin{iN},'loop4')
        loop4 = varargin{iN+1};
        ismandatoryarg(m_loop4)=1;
    end
    if strcmpi(varargin{iN},'loop3')
        CompiledHeader.LoopNum(3) = varargin{iN+1};
        ismandatoryarg(m_loop3)=1;
    end
    if strcmpi(varargin{iN},'loop2')
        CompiledHeader.LoopNum(2) = varargin{iN+1};
        ismandatoryarg(m_loop2)=1;
    end
    if strcmpi(varargin{iN},'loop1')
        CompiledHeader.LoopNum(1) = varargin{iN+1};
        ismandatoryarg(m_loop1)=1;
    end
    if strcmpi(varargin{iN},'datatype')
        datatype = varargin{iN+1};
        switch datatype
            case 'ushort'
                datasize = 2;
            case 'uint32'
                datasize = 4;
            case 'double'
                datasize = 8;
        end
        isdatatypeargin=1;
    end
    if strcmpi(varargin{iN},'compilesub')
        isCompileSubHeader = varargin{iN+1};
    end
    if strcmpi(varargin{iN},'forcereading')
        if(ischar(varargin{iN+1})||isstring(varargin{iN+1}))
            varargin{iN+1}=string2boolean(varargin{iN+1});
        end
        ForceReading = logical(varargin{iN+1});
    end
    if strcmpi(varargin{iN},'nSource')
        nDet = varargin{iN+1};
        ismandatoryarg(m_nsource)=1;
    end
    if strcmpi(varargin{iN},'nDet')
        nDet = varargin{iN+1};
        ismandatoryarg(m_ndet)=1;
    end
    if strcmpi(varargin{iN},'nBoard')
        nBoard = varargin{iN+1};
        ismandatoryarg(m_nboard)=1;
    end
    if strcmpi(varargin{iN},'nBin')
        nBin = varargin{iN+1};
        ismandatoryarg(m_nbin)=1;
    end
end

if(~isdatatypeargin)
    datatry = {'ushort','uint32','double'};
else
    datatry={datatype};
end
if ~all(ismandatoryarg([m_nboard m_ndet m_nsource])) && SkipSub==0
    for itry = 1:numel(datatry)
        frewind(fid);
        fread(fid,HeadLen,'uint8');
        out = false;
        BuffParms = []; Parms = [];
        while(out==false)
            SubRaw=fread(fid,SubLen,'uint8');
            if isempty(SubRaw), break; end
            CompSub = FillSub(SubRaw); fread(fid,nBin,datatry{itry});
            ActParms = [CompSub.Source CompSub.Det CompSub.Board];
            if(isequal(BuffParms,ActParms))
                out = true;
            else
                if isempty(BuffParms)
                    BuffParms = ActParms;
                end
                if isempty(Parms)
                    Parms = ActParms;
                else
                    Parms(end+1,:) = ActParms;
                end
            end
        end
        nSource = numel(categories(categorical(Parms(:,1))));
        nDet = numel(categories(categorical(Parms(:,2))));
        nBoard = numel(categories(categorical(Parms(:,3))));
        NumLoop=CompiledHeader.LoopNum;
        info=dir(FilePath);
        datatype = datatry{itry};
        if info.bytes == (HeadLen + prod(NumLoop)*(nBoard*nDet*nSource)*(SubLen+nBin*2))
            datasize = 2;
            break;
        end
        if info.bytes == (HeadLen + prod(NumLoop)*(nBoard*nDet*nSource)*(SubLen+nBin*4))
            datasize = 4;
            break;
        end
        if info.bytes == (HeadLen + prod(NumLoop)*(nBoard*nDet*nSource)*(SubLen+nBin*8))
            datasize = 8;
            break;
        end
        if(ForceReading==true&&isdatatypeargin)
            break;
        end
        if (ForceReading==false)&&itry==numel(datatry)
            errordlg({'Can''t handle sizemismatch. Insert more argin' 'Or use (...''forcereading'',''true'') argin'}); Data = [];
            fclose(fid);
            return;
        end
        if (ForceReading==true)&&itry==numel(datatry)
            fh = figure('NumberTitle','off','Name','Choose type','Toolbar','none','menubar','none','HandleVisibility','off','Units','normalized','Position',[0.5 0.5 0.1 0.3]);
            movegui(fh,'center');
            uph = uipanel(fh,'Title','Choose type','units','normalized','position',[0 0 1 1]);
            bg = uibuttongroup(uph,'Visible','on','Position',[0 0 1 1]);
            uicontrol(bg,'style','radiobutton','String','ushort','units','normalized','position',[0 0 1 0.5]);
            uicontrol(bg,'style','radiobutton','String','uint32','units','normalized','position',[0 1/3 1 0.5]);
            uicontrol(bg,'style','radiobutton','String','double','units','normalized','position',[0 2/3 1 0.5]);
            uicontrol(uph,'style','pushbutton','String','Ok','units','normalized','position',[0.5 0 0.5 0.1],'Callback',@AssignDataType);
            waitfor(fh,'Tag');
            close(fh);
            
            fclose(fid);
            fid=fopen(FilePath,'rb');
            fread(fid,HeadLen,'uint8');
            out = false;
            BuffParms = []; Parms = [];
            while(out==false)
                SubRaw=fread(fid,SubLen,'uint8');
                if isempty(SubRaw), break; end
                CompSub = FillSub(SubRaw); fread(fid,nBin,datatype);
                ActParms = [CompSub.Source CompSub.Det CompSub.Board];
                if(isequal(BuffParms,ActParms))
                    out = true;
                else
                    if isempty(BuffParms)
                        BuffParms = ActParms;
                    end
                    if isempty(Parms)
                        Parms = ActParms;
                    else
                        Parms(end+1,:) = ActParms;
                    end
                end
            end
            nSource = numel(categories(categorical(Parms(:,1))));
            nDet = numel(categories(categorical(Parms(:,2))));
            nBoard = numel(categories(categorical(Parms(:,3))));
            break;
        end
    end
    
    
    
end
fclose(fid);
NumLoop=CompiledHeader.LoopNum;
CompiledHeader.NumBoard = nBoard;
CompiledHeader.NumDet = nDet;
CompiledHeader.NumSource = nSource;

fid=fopen(FilePath,'rb');
Head=fread(fid,HeadLen,'uint8');

A=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),nBoard,nDet,nSource,nBin);
Sub=zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1),nBoard,nDet,nSource,SubLen);
%CompiledSub = zeros(NumLoop(5),NumLoop(4),NumLoop(3),NumLoop(2),NumLoop(1));
isbreakcond = false;
try
    for il5= 1:NumLoop(5)
        for il4= 1:NumLoop(4)
            for il3= 1:NumLoop(3)
                for il2=1:NumLoop(2)
                    for il1=1:NumLoop(1)
                        for iB = 1:nBoard
                            for iD = 1: nDet
                                for iS = 1:nSource
                                    if SkipSub==0
                                        BuffSub = fread(fid,SubLen,'uint8');
                                        if ~isempty(BuffSub)&&numel(BuffSub)==SubLen
                                            Sub(il5,il4,il3,il2,il1,iB,iD,iS,:)=BuffSub;
                                        else
                                            warning('backtrace','off')
                                            warning('Reading interrupted at:');
                                            warning(strcat('Loop5: ',num2str(il5),'/',num2str(NumLoop(5))));
                                            warning(strcat('Loop4: ',num2str(il4),'/',num2str(NumLoop(4))));
                                            warning(strcat('Loop3: ',num2str(il3),'/',num2str(NumLoop(3))));
                                            warning(strcat('Loop2: ',num2str(il2),'/',num2str(NumLoop(2))));
                                            warning(strcat('Loop1: ',num2str(il1),'/',num2str(NumLoop(1))));
                                            warning(strcat('NumBoard: ',num2str(iB),'/',num2str(nBoard)));
                                            warning(strcat('NumDet: ',num2str(iD),'/',num2str(nDet)));
                                            warning(strcat('NumSource: ',num2str(iS),'/',num2str(nSource)));
                                            warning('Output data will have the dimension specified in TRS settings (Header)');
                                            warning('backtrace','on')
                                            isbreakcond = true;
                                            break;
                                        end
                                        if (isCompileSubHeader)
                                            CompiledSub(il5,il4,il3,il2,il1,iB,iD,iS) = FillSub(squeeze(Sub(il5,il4,il3,il2,il1,iB,iD,iS,:)));
                                        end
                                    else
                                        BuffSub = 0;
                                    end
                                    TrashData = fread(fid,nBin,datatype);
                                    if ~isempty(TrashData)&&numel(TrashData)==nBin
                                        Sub(il5,il4,il3,il2,il1,iB,iD,iS,:)=BuffSub;
                                    else
                                        warning('backtrace','off')
                                        warning('Reading interrupted at:');
                                        warning(strcat('Loop5: ',num2str(il5),'/',num2str(NumLoop(5))));
                                        warning(strcat('Loop4: ',num2str(il4),'/',num2str(NumLoop(4))));
                                        warning(strcat('Loop3: ',num2str(il3),'/',num2str(NumLoop(3))));
                                        warning(strcat('Loop2: ',num2str(il2),'/',num2str(NumLoop(2))));
                                        warning(strcat('Loop1: ',num2str(il1),'/',num2str(NumLoop(1))));
                                        warning(strcat('NumBoard: ',num2str(iB),'/',num2str(nBoard)));
                                        warning(strcat('NumDet: ',num2str(iD),'/',num2str(nDet)));
                                        warning(strcat('NumSource: ',num2str(iS),'/',num2str(nSource)));
                                        warning('Output data will have the dimension specified in TRS settings (Header)');
                                        warning('backtrace','on')
                                        isbreakcond = true;
                                        break;
                                    end
                                    A(il5,il4,il3,il2,il1,iB,iD,iS,:)=TrashData;
                                end
                                if isbreakcond, break; end
                            end
                            if isbreakcond, break; end
                        end
                        if isbreakcond, break; end
                    end
                    if isbreakcond, break; end
                end
                if isbreakcond, break; end
            end
            if isbreakcond, break; end
        end
        if isbreakcond, break; end
    end
    
    fclose(fid);
catch ME
    fclose(fid);
    %ME = MException('FileReading:GeneralError',{'Error reading the file:',ME.message, ' Encountered error at Loop1: ', num2str(il1), ' Loop2: ', num2str(il2), ' Loop3: ', num2str(il3), 'Loop4: ', num2str(il4),' Loop5: ', num2str(il5)});
    throw(ME);
    return
end

Data = squeeze(A);
switch NumArgOut
    case 1
        varargout{1}= Head;
    case 2
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
    case 3
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
    case 4
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
        if isCompileSubHeader == false
            varargout{4} = [];
        else
            varargout{4} = squeeze(CompiledSub);
        end
    case 5
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
        if isCompileSubHeader == false
            varargout{4} = [];
        else
            varargout{4} = squeeze(CompiledSub);
        end
        varargout{5} = Sub;
    case 6
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
        if isCompileSubHeader == false
            varargout{4} = [];
        else
            varargout{4} = squeeze(CompiledSub);
        end
        varargout{5} = Sub;
        varargout{6} = datasize;
    case 7
        varargout{1} = Head;
        varargout{2} = CompiledHeader;
        varargout{3} = squeeze(Sub);
        if isCompileSubHeader == false
            varargout{4} = [];
        else
            varargout{4} = squeeze(CompiledSub);
        end
        varargout{5} = Sub;
        varargout{6} = datasize;
        varargout{7} = datatype;
end

end

function FS = FillSub(Sub)
FieldNames = {'Geom','Source','Fiber','Det','Board','Coord','Pad','Xf','Yf','Zf','Rf','Xs','Ys','Zs','Rs','Rho','TimeNom','TimeEff'...
    'n','Loop','Acq','Page','RoiNum','RoiFirst','RoiLast','RoiLambda','RoiPower'};
C = 1; CT = 1;
L = 4; LT = 2;
D = 8; DT = 3;
S = 2; ST = 4;
Unit = 1; %byte
FieldSize = [C,C,C,C,C,C,C,D,D,D,D,D,D,D,D,D,D,D,D,3*L,L,L,C,4*S,4*S,4*D,4*D]./Unit;
FieldType = [CT,CT,CT,CT,CT,CT,CT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,DT,LT,LT,LT,CT,ST,ST,DT,DT];
iS = 1;
iF = 1;

while iS<=length(Sub)
    iStart = iS;
    iStop = iS + FieldSize(iF);
    iStop = iStop-1;
    iS = iStop+1;
    RawData = Sub(iStart:iStop);
    
    switch FieldType(iF)
        case CT
            Data = RawData;
            if strcmp(FieldNames{iF},'Geom')
                if Data == 0
                    Data = 'REFL';
                else
                    Data = 'TRASM';
                end
            end
            if strcmp(FieldNames{iF},'Coord')
                if Data == 0
                    Data = 'CART';
                else
                    Data = 'POLAR';
                end
            end
        case LT
            Data = cast(typecast(uint8(RawData),'int32'),'double');
        case DT
            Data = typecast(uint8(RawData),'double');
        case ST
            Data = cast(typecast(uint8(RawData),'int16'),'double');
    end
    
    FS.(FieldNames{iF}) = Data;
    iF = iF +1;
end

end


function FH = FillHeader(Head)
FieldNames = {'Ver','SubHeader','SubHeaderVer','SizeHeader','SizeSubHeader','SizeData','Kind','Appl'...
    ,'Oma','Date','Time','LoopHome','LoopFirst','LoopLast','LoopDelta','LoopNum','McaChannNum','PageNum'...
    'FrameNum','RamNum','McaTime','McaFactor','MeasNorm','LabelName','LabelContent','Constn','ConstRho'...
    'ConstThick','MammHeader','MammIdxFirst','MammIdxLast','MammIdxTop','MammRateMid','MammRateHigh'};
C = 1; CT = 1;
L = 4; LT = 2;
D = 8; DT = 3;
S = 2; ST = 4;
Unit = 1; %byte
FieldSize = [2*S,L,L,L,L,L,L,L,L,11*C,9*C,3*L,3*L,3*L,3*L,3*L,L,L,L,L,D,D,L,192*C,352*C,D,D,D,L,2*L,2*L,2*L,2*L,2*L]./Unit;
FieldType = [ST,LT,LT,LT,LT,LT,LT,LT,LT,CT,CT,LT,LT,LT,LT,LT,LT,LT,LT,LT,DT,DT,LT,CT,CT,DT,DT,DT,LT,LT,LT,LT,LT,LT];
iS = 1;
iF = 1;

while iS<=length(Head)
    iStart = iS;
    iStop = iS + FieldSize(iF);
    iStop = iStop-1;
    iS = iStop+1;
    RawData = Head(iStart:iStop);
    %isCHAR = false;
    switch FieldType(iF)
        case CT
            if strcmp(FieldNames{iF},'LabelName')
                Data = reshape(RawData,[12,16])';
                Data = char(Data);
                Data = string(Data);
                %isCHAR = true;
            end
            if strcmp(FieldNames{iF},'LabelContent')
                Data = reshape(RawData,[22,16])';
                Data = char(Data);
                Data = string(Data);
                %isCHAR = true;
            end
            if strcmp(FieldNames{iF},'Date')
                Data = char(RawData');
                %isCHAR = true;
            end
            if strcmp(FieldNames{iF},'Time')
                Data = char(RawData');
                %isCHAR = true;
            end
            %isCHAR = false;
        case LT
            Data = cast(typecast(uint8(RawData),'int32'),'double');
            if strcmp(FieldNames{iF},'Kind')
                switch Data
                    case 0
                        Data = 'Measure';
                    case 1
                        Data = 'System';
                    case 2
                        Data = 'Simul';
                end
            end
            if strcmp(FieldNames{iF},'Appl')
                switch Data
                    case 0
                        Data = 'Diff';
                    case 1
                        Data = 'Mamm';
                    case 2
                        Data = 'Oxym';
                    case 3
                        Data = 'Fluo';
                    case 4
                        Data = 'Spec';
                end
            end
            if strcmp(FieldNames{iF},'SubHeader')
                switch Data
                    case 0
                        Data = false;
                    case 1
                        Data = true;
                end
            end
            if strcmp(FieldNames{iF},'Oma')
                switch Data
                    case 0
                        Data = false;
                    case 1
                        Data = true;
                end
            end
            if strcmp(FieldNames{iF},'MammHeader')
                switch Data
                    case 0
                        Data = false;
                    case 1
                        Data = true;
                end
            end
        case DT
            Data = typecast(uint8(RawData),'double');
        case ST
            if strcmp(FieldNames{iF},'Ver')
                Data = cast(typecast(uint8(RawData(1:2)),'int16'),'double');
                Data = [Data cast(typecast(uint8(RawData(3:4)),'int16'),'double')];
            else
                Data = cast(typecast(uint8(RawData),'int16'),'double');
            end
    end
    
    FH.(FieldNames{iF}) = Data;
    iF = iF +1;
end

end
function [output]=string2boolean(string)
if strcmp(string,'false')
    output = false;
else
    output = true;
end
end
function AssignDataType(src,~)
ph = src.Parent;
rbh = findobj(ph,'style','radiobutton');
assignin('caller','datatype',rbh(logical([rbh.Value])).String);
fh = ancestor(src,'figure');
fh.Tag = '1';
end