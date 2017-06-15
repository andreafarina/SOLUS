function y = WindowTPSF(TPSF,twin)
%
% take a set of TPSFs and window them
%
% dimTPSF = numel(size(TPSF));
% nt = size(TPSF,1);
% nwin = size(twin,1);
% switch dimTPSF
%     case 2
%         ndat = size(TPSF,2);
%         y = zeros(nwin,ndat);
%         for w = 1:nwin
%             y(w,:) = sum(TPSF(twin(w,1):twin(w,2),:));
%         end
%    
%     case 4
%         ndat = [size(TPSF,2) size(TPSF,3) size(TPSF,4)];
%         y = zeros([nwin ndat]);
%         for w = 1:nwin
%             y(w,:,:,:) = sum(TPSF(twin(w,1):twin(w,2),:,:,:));
%         end
% end

size_y = size(TPSF);
nwin = size(twin,1);
size_y(1) = nwin;
y = zeros(size_y);
for w = 1:nwin
             y(w,:) = sum(TPSF(twin(w,1):twin(w,2),:),1);
end

    
        
       
    
        
        