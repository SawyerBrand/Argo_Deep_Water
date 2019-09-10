function CO2_calculation(fn,titl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was made by Sawyer Brand, SIO/SOCCOM            %
% during research into the Argentine Basin mesoscale activities %
% and CO2/Heat Flux in that region                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        f = load(fn);
        d = load('ERA_Interim_Winds_062019.mat');
        c = load('CapeGrimObs_062019.mat');
        
        pCO2_ocn = f.data.pCO2; 
        
        % Function used by Alison Gray:
        % F = k*Ko*(pCO2_ocn - pCO2_atm)
        
        F1 = (0.251.*((d.wind.testWIND).^2).*((668/660).^(-0.5))).*(pCO2_ocn(1,5)-c.pCO2(9,:));
        

        for i = 1:1:size(f.idat.pCO2)
            F_CO2(i,:) = (0.251.*((d.wind.testWIND).^2).*((668/660).^(-0.5))).*(pCO2_ocn(i,1)-c.pCO2(9,:));
        end
        
        save([titl,'_CO2Flux.mat'])

end