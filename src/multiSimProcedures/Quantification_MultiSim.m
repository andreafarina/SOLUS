if exist('isComingFromInterface','var')
    warning('off','backtrace');
    if ~isempty(Q)
    if (~EXP_DATA) == 1
        if exist('errXP','var'), F(errXP).Value(iL) = Q.COM.error(1); else, warning('errXP not calculated'); end
        if exist('errYP','var'), F(errYP).Value(iL) = Q.COM.error(2);else, warning('errYP not calculated'); end
        if exist('errZP','var'), F(errZP).Value(iL) = Q.COM.error(3);else, warning('errZP not calculated'); end
        if exist('relMUAP','var'), F(relMUAP).Value(iL) = Q.max.rel_error;else, warning('relMUAP not calculated'); end
        if exist('relVOL','var'), F(relVOL).Value(iL) = Q.volume.rel_error;else, warning('relVOL not calculated'); end
        if exist('relVOLGauss','var'), F(relVOLGauss).Value(iL) = Q.volumeG.rel_error;else, warning('relVOLGauss not calculated'); end
    else
        warning('No quantification without reference');
        for iF=1:numel(F), F(iF).Value(iL) = 0; end %#ok<SAGROW>
    end
    else
        warning('No quantification due to error in fittin procedure');
    end
    warning('on','backtrace');
end