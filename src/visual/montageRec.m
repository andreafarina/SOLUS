function montageRec(X,varargin)
%==========================================================================
%%                         OPTIONS   
%==========================================================================
displayRange = [0,1];
for ii=1:size(varargin,2)/2
   if strcmpi(varargin{1,2*ii-1},'DisplayRange'),  displayRange = varargin{1,2*ii}; end
end
m = displayRange(1);
M = displayRange(2);
%==========================================================================
%%                           MAIN
%==========================================================================
ind0 = find(X==0);

% % reaffect the min
X(X< m + 1/64*(M-m)) = m + 1/64*(M-m);

% reaffect 0
%X(X==0) = 0;

%
% X = flipdim(X,3);
% X = permute(X,[1,2,4,3]);
% montage(X,varargin{:})


% reaffect the min
%X(X< m) = m;

% reaffect 0
X(ind0) = 0;%m-1/64*(M-m);

for ii=1:size(varargin,2)/2
   if strcmpi(varargin{1,2*ii-1},'DisplayRange'),  
       varargin{1,2*ii} = [m-1/64*(M-m) M];
   end
end

%
%X = flipdim(X,3);
X = permute(X,[1,2,4,3]);
montage(X,varargin{:})

jet2 = jet;%parula;%hot;%jet,parula;
jet2(1,:) = 1;

colormap(gca,jet2)