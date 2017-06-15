function twin = CreateTimeWindows(nstep,roi,kind,param)

% Synopsis: twin = CreateTimeWindows(nstep,roi,kind,param)
% Create a n x 2 matrix containing for each row the boundary of the time
% window
% nstep: number of steps set in the TOAST forward problem
% roi:  [step_in,step_last] region of interest
% kind:
% 'even' --> evenly spaced time windows in roi
%        param --> number of windows
% 'integral'   --> integral of all the curve in the roi


%% consistency check
%if ((roi(1) > nstep) || (roi(2) > nstep) || (roi(1) * roi(2) <0) || (roi(1)>roi(2)))
if ((roi(1) * roi(2) <0) || (roi(1)>roi(2)))

    error('error! check time_windows');
end
   
switch lower(kind)
    case 'even'
        if nargin < 3
            error('Set the number of windows!')
        end
        nt_win = round(diff(roi)/param);
        a = roi(1):nt_win:roi(2);
        b = a + (nt_win-1);
        a = a(a<=roi(2));
        b = b(b<=roi(2));
        minl = min(length(a),length(b));
        a = a(1:minl);
        b = b(1:minl);
    case 'integral'
        a = roi(1);
        b = roi(2);
end
twin = uint16([a',b']);
