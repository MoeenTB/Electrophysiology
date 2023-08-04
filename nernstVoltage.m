function E_i = nernstVoltage(C_i, C_e, z_i)
    % RT/F = 25.8 mV
    E_i = 25.8/abs(z_i)*log(C_e/C_i);
end

