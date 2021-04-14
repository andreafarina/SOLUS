function ext_c = LoadSpectra(spectra_file,lambda,conc,density)
    nC = numel(conc);
    warning off
    %Data = readtable(spectra_file);
    Data=dlmread([spectra_file '.txt'],'\t',1,0);
    warning on
    %Dlambda = Data{:,1};
    Dlambda = Data(:,1);
    lambdaIdx=ismember(Dlambda,lambda);
    %ext_c = Data{lambdaIdx,2:nC+1};
    ext_c = Data(lambdaIdx,2:nC+1);
    ext_c = ext_c./10; %da cm-1 a mm-1
end