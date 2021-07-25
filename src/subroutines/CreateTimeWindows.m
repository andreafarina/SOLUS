function twin = CreateTimeWindows(nstep,roi,kind,param)

% Synopsis: twin = CreateTimeWindows(nstep,roi,kind,param)
% Create a n x 2 matrix containing for each row the boundary of the time
% window
% nstep: number of steps set in the TOAST forward problem
% roi:  [step_in,step_last] region of interest
% kind:
%   'even'      --> evenly spaced time windows in roi
%        param  --> number of windows
%   'integral'  --> integral of all the curve in the roi
%   'delay'     --> the first time is taken as 'even' and the last time is
%                   the end of the vector.
% Andrea Farina - CNR/IFN - Politecnico di Milano   6/10/17

%% consistency check
%if ((roi(1) > nstep) || (roi(2) > nstep) || (roi(1) * roi(2) <0) || (roi(1)>roi(2)))
if ((roi(1) * roi(2) <0) || (roi(1)>roi(2)))

    error('error! check time_windows');
end
if param == 1
    kind = 'integral';
end   
switch lower(kind)
    case ('even')
        if nargin < 3
            error('Set the number of windows!')
        end
        while rem(diff(roi),param)~=0
            if rem(diff(roi),param) < floor(0.01*param) && diff(roi) > param
                if mod(rem(diff(roi),param),2) == 1 
                    roi(2) = roi(2) - 1;
                else
                    roi(1) = roi(1) + 1;
                end
            elseif rem(diff(roi),param) >= floor(0.01*param) || diff(roi) <= param
                if mod(rem(diff(roi),param),2) == 1 
                    roi(2) = roi(2) + 1;
                else
                    roi(1) = roi(1) - 1;
                end
            else
                roi(2) = roi(2) + 1;
            end           
        end
        if roi(1) <= 0 %translate roi(1) if its negative or roi(2) bigger than nstep            
            roi(2) = roi(2)-roi(1)+1;
            roi(1) = 1;
        end
        if roi(2) > nstep
            roi(1) = roi(1) - (roi(2)-nstep);
            roi(2) = nstep;
        end
        if roi(1)<=0 && roi(2)>nstep
            error('Decrease number of time windows or relax roi conditions!')
        end
        
        fprintf(['<strong>Actual roi: ', num2str(roi),'</strong>\n'])
        nt_win = round(diff(roi)/param);
        a = roi(1):nt_win:roi(2);
        b = a + (nt_win-1);
        %a = a(a<=roi(2));
        a = min(a,roi(2));
        %b = b(b<=roi(2));
        b = min(b,roi(2));
        %minl = min(length(a),length(b));
        %a = a(1:minl);
        %b = b(1:minl);
        nmax = min(param,length(a));
        a = a(1:nmax);
        b = b(1:nmax);
    case ('delay')
        if nargin < 3
            error('Set the number of windows!')
        end
        nt_win = round(diff(roi)/param);
        a = roi(1):nt_win:roi(2);
        b = a + (nt_win-1);
        %a = a(a<=roi(2));
        a = min(a,roi(2));
        %b = b(b<=roi(2));
        b = min(b,roi(2));
        %minl = min(length(a),length(b));
        %a = a(1:minl);
        %b = b(1:minl);
        b = nstep*ones(size(b));
        
    case 'integral'
        a = roi(1);
        b = roi(2);
end
twin = uint16([a',b']);
