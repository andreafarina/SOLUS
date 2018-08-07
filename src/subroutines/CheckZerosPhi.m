function CheckZerosPhi(phi,radiometry,phi_name)
spacer = '-----'; 
nmeas = size(phi,2)/radiometry.nL;
disp(['Checking for zeros in ', phi_name])
for inl = 1:radiometry.nL
    disp(spacer);
    disp(['Wavelength: ' num2str(radiometry.lambda(inl)) ' nm']);
    meas_set = (1:nmeas)+(inl-1)*nmeas;
    dummy_phi = phi(:,meas_set);
    if sum(dummy_phi(:)<0) > 0
        warning([num2str(sum(dummy_phi(:)<0)),' elements of ' phi_name ' < 0!']);
    end
end
end