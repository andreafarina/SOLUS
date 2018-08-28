function [factor,y,std] = CutCounts(ch,dt,data,MaxCount,radiometry,...
            modality, NumDelay)
% calculate the multiplicative factor for converting 
% TPSFs curves(without Poisson noise) expressed in photon-counts into 
% different kind of data representative of gated and non-gated measurements
% size(data) = [numel(t),num_TPSF]
% modality = 0 
%   classic mode: the highest power TPSF is identified, a constant 
%   factor is calculated and all curves are rescaled. In this case
%   radiometry, NumDelay are not effective.
% modality = 1
%   gated mode: A Time-gated measurement is simulated.
%   Measurements on each dealay are rescaled accordingly 
%   to MaxCount. In particular, if:
%        radiometry==1, the count-rate for each delay is cut to 
%             MaxCount if photons are available, otherwise no;
%        radiometry==0, the count-rate for each delay is cut to 
%             MaxCounts in any case. 
% 
%   CutCounts(time,data,MaxCount):--> modalty = 0
% 
% A. Farina - CNR-IFN  - Dip di Fisica - Politecnico di Milano 27/09/17
% 
% =========================================================================

%dt = t(2)-t(1); % assimung uniform sampling of the time axis
t = ch * dt;

verbosity = 0;
%% remodulate MaxCounts accordingly to the number of delays
if modality == 1
    MaxCount = MaxCount/NumDelay;
end
%% identify the highest power TPSF
[m,~] = max(sum(data,1));
%yref = data(:,j);
if verbosity == 1
    figure(1),subplot(2,2,1),semilogy(t,data),
end
if nargin<5
    factor = 1./m .* MaxCount;
    y = bsxfun(@times,data,factor);
    y = poissrnd(round(y));
    std = sqrt(y);
   return
end

switch modality 
    case 0 % classic mode
        factor = 1./m .* MaxCount;
        y = bsxfun(@times,data,factor);
        y = poissrnd(round(y));
        std = sqrt(y);
    case 1 % gated
        twin = CreateTimeWindows(numel(t),[ch(1),ch(end)],'even',NumDelay);
        if NumDelay == 1
            twin = [ch(1),ch(end)];
        end
               
        factor = zeros(size(data));
        for i = 1:NumDelay
            idx = (twin(i,1):twin(i,2));
            area = sum(data(twin(i,1):end,:),1);
            %if radiometry == 1
            %    area(i) = max(sum(yref(twin(i,1):end)),MaxCount);
            %else
            %    area(i) = sum(yref(twin(i,1):end));
            %end
            if sum(area==0) > 0
                disp('CutCounts: area=0 error');
                break
            end
            if radiometry == 1
                factor(idx,:) = repmat(min(area,MaxCount)./area,...
                    [numel(idx) 1]);
            else
                factor(idx,:) = repmat((MaxCount)./area,...
                    [numel(idx) 1]);
            end
        end 
        % manage zero elements
        %idz = find(factor>0, 1, 'last' );
        y = poissrnd(round(data.*factor));
        %y = bsxfun(@times,y1,factor(idz)./factor');
        std = sqrt(y);%bsxfun(@times, sqrt(y1), (factor(idz)./factor'));
        
        if verbosity == 1
            figure%,subplot(2,2,2),semilogy(t,y1),title('sliced')
            subplot(2,2,3),semilogy(t,y,'b',t,std,'r'),
            title('reconstructed'),legend('recon','sigma');
            subplot(2,2,4),plot(t,factor),title('factor')
        end
end


