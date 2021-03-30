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
if size(twin,3) > 1
    for i = 1:size_y(2)
        for w = 1:nwin
            y(w,i) = sum(TPSF(squeeze(twin(w,1,i)):squeeze(twin(w,2,i)),i),1,'omitnan');
        end
    end

else
    
    for w = 1:nwin
        y(w,:) = sum(TPSF(squeeze(twin(w,1)):squeeze(twin(w,2)),:),1,'omitnan');
    end

end
    
        
       
    
        
        