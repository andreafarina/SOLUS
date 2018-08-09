% meas_homo = meas_homo - mean(meas_homo(bkg.start:bkg.stop,:,:),1);
% meas_homo(meas_homo<0) = 0;
% meas_hete = meas_hete - mean(meas_hete(bkg.start:bkg.stop,:,:),1);
% meas_hete(meas_hete<0) = 0;
full_raw_contrast_cw = zeros(64,numel(lambda));
full_raw_contrast_td = zeros(Nchannel,64,numel(lambda));
raw_contrast_cw = (squeeze(sum(meas_homo,1))-squeeze(sum(meas_hete,1)))./squeeze(sum(meas_homo,1));
raw_contrast_cw(isnan(raw_contrast_cw)) = 0;
raw_contrast_td = (meas_homo - meas_hete)./meas_homo;
raw_contrast_td(isnan(raw_contrast_td)) = 0;
for iw = 1:numel(lambda)
    idm = 1;
    for id = 1:numel(dmask(:))
       if dmask(id)==1
           full_raw_contrast_cw(id,iw) = raw_contrast_cw(idm,iw);
           full_raw_contrast_td(:,id,iw) = raw_contrast_td(:,idm,iw); 
            idm = idm +1;
       end
    end
end
full_raw_contrast_cw=reshape(full_raw_contrast_cw,8,8,numel(lambda));
full_raw_contrast_td=reshape(full_raw_contrast_td,Nchannel,8,8,numel(lambda));
clim_cw = [min(full_raw_contrast_cw(:)) max(full_raw_contrast_cw(:))];
clim_td = [min(full_raw_contrast_td(:)) max(full_raw_contrast_td(:))];
h=figure(3);
h.Name = 'CW';
for iw = 1:numel(lambda)
    subplot(2,4,iw)
    imagesc(full_raw_contrast_cw(:,:,iw),clim_cw)
    axis image
    colormap pink, shading interp; colorbar
    title(num2str(lambda(iw)))
end
%%
% figure(4)
% for iw = 1:numel(lambda)
%     subplot(2,4,iw)
%     [val,pos]=max(irf(:,iw));
%     full_raw_contrast_td = movmean(full_raw_contrast_td,[20 20],1);
%     pos = 1;
%     plot(1:numel(pos:Nchannel),squeeze(full_raw_contrast_td(pos:end,:,1,iw)))
%     ylim([-1 1]);
%     title(num2str(lambda(iw)))
%     legend(num2str([1:8]'),'location','best');
% end
%%
return
figure
sel_roi = SelectROI(meas_homo(:,:,1),irf(:,1));
num_tw = round(numel(sel_roi(1):sel_roi(2))/round(500/factor));
twin = CreateTimeWindows(4096,sel_roi,'even',num_tw);
ShowTimeWindows(meas_homo(:,:,1),twin,factor);
for iw = 1:numel(lambda)
    for id = 1:28
        for itw = 1:num_tw
            raw_contrast_gate(itw,id,iw) = squeeze((sum(meas_homo(twin(itw,1):twin(itw,2),id,iw),1)-sum(meas_hete(twin(itw,1):twin(itw,2),id,iw),1))./sum(meas_homo(twin(itw,1):twin(itw,2),id,iw),1));
        end
    end
end
full_raw_contrast_gate = zeros(num_tw,64,numel(lambda));
for iw = 1:numel(lambda)
    idm = 1;
    for id = 1:numel(dmask(:))
       if dmask(id)==1
           full_raw_contrast_gate(:,id,iw) = raw_contrast_gate(:,idm,iw); 
            idm = idm +1;
       end
    end
end
full_raw_contrast_gate = reshape(full_raw_contrast_gate,num_tw,8,8,numel(lambda));
clim_gate = [min(full_raw_contrast_gate(:)) max(full_raw_contrast_gate(:))];
for itw = 1:num_tw
    h=FFS;
    h.Name = num2str(twin(itw,:));
    
for iw = 1:numel(lambda)
    subplot(2,4,iw)
    imagesc(squeeze(full_raw_contrast_gate(itw,:,:,iw)),clim_gate)
    axis image
    colormap pink, shading interp; colorbar
    title(num2str(lambda(iw)))
end
end
%tilefigs;
%%
return
full_raw_contrast_gate = zeros(64,numel(lambda));
for iw = 1:numel(lambda)
    for id = 1:28
        [val,pos]=max(meas_homo(:,id,iw));
        normalized_data = meas_homo(:,id,iw)./max(meas_homo(:,id,iw));
        idxs = find((and(logical(normalized_data<70/100),logical(normalized_data>20/100))));
        roi = idxs(idxs > pos);
        semilogy(1:Nchannel,meas_homo(:,id,iw),1:Nchannel,meas_hete(:,id,iw))
        vline([roi(1) roi(end)],{'-k' '-k'});
        pause(0.05);
        Wave(iw).Pos(id).ROI = roi;
        raw_contrast_gate(id,iw) = squeeze((sum(meas_homo(roi,id,iw),1)-sum(meas_hete(roi,id,iw),1))./sum(meas_homo(roi,id,iw),1));
    end
end
for iw = 1:numel(lambda)
    idm = 1;
    for id = 1:numel(dmask(:))
       if dmask(id)==1
           full_raw_contrast_gate(id,iw) = raw_contrast_gate(idm,iw); 
            idm = idm +1;
       end
    end
end
full_raw_contrast_gate = reshape(full_raw_contrast_gate,8,8,numel(lambda));
figure(5)
clim_gate = [min(full_raw_contrast_gate(:)) max(full_raw_contrast_gate(:))];
for iw = 1:numel(lambda)
    subplot(2,4,iw)
    imagesc(full_raw_contrast_gate(:,:,iw),clim_gate)
    axis image
    colormap pink, shading interp; colorbar
    title(num2str(lambda(iw)))
end
%%
return
semilogy(1:Nchannel,squeeze(meas_hete(:,:,1)));
hold on
semilogy(1:Nchannel,irf(:,1));
[x,~] = ginput(2);
roi = round(x);
roi_vect(iw).pos(id).vect = roi(1):roi(2);
roi = roi(1):roi(2);
hold off
full_raw_contrast_gate = zeros(64,numel(lambda));
raw_contrast_gate = squeeze((sum(meas_homo(roi,:,:),1)-sum(meas_hete(roi,:,:),1))./sum(meas_homo(roi,:,:),1));
raw_contrast_gate(isnan(raw_contrast_gate)) = 0;
for iw = 1:numel(lambda)
    idm = 1;
    for id = 1:numel(dmask(:))
       if dmask(id)==1
           full_raw_contrast_gate(id,iw) = raw_contrast_gate(idm,iw); 
            idm = idm +1;
       end
    end
end
full_raw_contrast_gate = reshape(full_raw_contrast_gate,8,8,numel(lambda));
hold off
figure(5)
clim_gate = [min(full_raw_contrast_gate(:)) max(full_raw_contrast_gate(:))];
for iw = 1:numel(lambda)
    subplot(2,4,iw)
    imagesc(full_raw_contrast_gate(:,:,iw),clim_gate)
    axis image
    colormap pink, shading interp; colorbar
    title(num2str(lambda(iw)))
end