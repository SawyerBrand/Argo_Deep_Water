function interpolation_comparison(file)

    F = load(file); 
    
    subplot(1,2,1)
        
   
    hold on 
    title('First Profile (12700)')
    set(gca,'ydir','reverse')
    set(gca,'FontSize',20)
    
    xlabel('Real - Interpolated [C]')
    ylabel('Pressure [dbars]')
    for i = 1:2:2
        plot((F.idat.T(i,1:500)-F.data.T(i,:))+(i*0.5),F.data.P(i,:))
    end

    hold off
    
    
    %All profiles compariosn
    subplot(1,2,2)
    
    hold on 
    title('All Profiles (12700)')
    set(gca,'ydir','reverse')
    set(gca,'FontSize',20)
    
    xlabel('Real - Interpolated [C]')
    ylabel('Pressure [dbars]')
    
    for i = 1:2:size(F.data.T)
        plot((F.idat.T(i,1:500)-F.data.T(i,:))+(i*0.5),F.data.P(i,:))
    end

    hold off
    
end