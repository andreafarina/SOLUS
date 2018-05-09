if exist('isComingFromInterface','var')
    warning('off','backtrace');
    FN = fieldnames(Q);
    for in = 1:numel(FN)
        if ~isempty(Q.(FN{in}))
            if (~EXP_DATA) == 1
                if exist([FN{in} 'errXP'],'var'), F(eval([FN{in} 'errXP'])).Value(iL) = Q.(FN{in}).COM.error(1); else, warning('errXP not calculated'); end
                if exist([FN{in} 'errYP'],'var'), F(eval([FN{in} 'errYP'])).Value(iL) = Q.(FN{in}).COM.error(2);else, warning('errYP not calculated'); end
                if exist([FN{in} 'errZP'],'var'), F(eval([FN{in} 'errZP'])).Value(iL) = Q.(FN{in}).COM.error(3);else, warning('errZP not calculated'); end
                if exist([FN{in} 'relP'],'var'), F(eval([FN{in} 'relP'])).Value(iL) = Q.(FN{in}).max.rel_error;else, warning('relMUAP not calculated'); end
                if exist([FN{in} 'relVOL'],'var'), F(eval([FN{in} 'relVOL'])).Value(iL) = Q.(FN{in}).volume.rel_error;else, warning('relVOL not calculated'); end
                if exist([FN{in} 'relVOLGauss'],'var'), F(eval([FN{in} 'relVOLGauss'])).Value(iL) = Q.(FN{in}).volumeG.rel_error;else, warning('relVOLGauss not calculated'); end
            else
                warning('No quantification without reference');
                for iF=1:numel(F), F(iF).Value(iL) = 0; end 
            end
        else
            warning('No quantification due to error in fittin procedure');
        end
    end
    warning('on','backtrace');
end