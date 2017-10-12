warning('off','backtrace');
if (~EXP_DATA) == 1
    if exist('errXP','var'), F(errXP).Value(iL) = Q.COM.error(1); end
    if exist('errYP','var'), F(errYP).Value(iL) = Q.COM.error(2);end
    if exist('errZP.Value','var'), F(errZP).Value(iL) = Q.COM.error(3); end
    if exist('relMUAP','var'), F(relMUAP).Value(iL) = Q.max.error; end
    if exist('relVOL','var'), F(relVOL).Value(iL) = Q.volume.rel_error; end
else
    warning('No quantification without reference');
    for iF=1:numel(F), F(iF).Value(iL) = 0; end %#ok<SAGROW>
end
warning('on','backtrace');