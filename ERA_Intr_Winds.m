function ERA_Intr_Winds(filename)

        wind = ncreadall(filename);
        
        for i = 1:1:size(wind.time)
            wind.time(i,:) = (wind.time(i,:)/24)+693962;
            timy = double(wind.time(i,:));
            wind.dattime(i,:) = datestr(timy);
        end
        
        wind.wind = sqrt((wind.u10).^2+(wind.v10).^2);
        wind.longitude = wind.longitude-180;
        wind.testWIND = (wind.wind(193,167,1));
        
        F12700.May1 = (wind.wind(193,167,1));
        F12700.May10 = (wind.wind(193,167,10));
        F12700.May20 = (wind.wind(192,167,20));
        F12700.May31 = (wind.wind(192,167,31));
        
        outfile = ['ERA_Interim_Winds_062019.mat'];
        
        save(outfile,'wind','F12700')



end