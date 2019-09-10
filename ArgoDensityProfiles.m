function ArgoDensityProfiles(fn,titl)

Data = load(fn);

figure(1)
hold on

title([titl,' Raw Temperature vs Potential Density'])
xlabel('In-Situ Temp [C]')
ylabel('Potential Density')

set(gca,'Ydir','reverse')
set(gca,'FontSize',16)

for i = 1:size(Data.data.T,1)
    plot((Data.data.T(i,63:69)),Data.data.pden(i,63:69))
end 

saveas(gcf,[titl,'_RawT_Pden.jpg'])

hold off


figure(2)
hold on

title([titl,' Raw Abs Salinity vs Potential Density'])
xlabel('Absolute Salinity [g/kg]')
ylabel('Potential Density')

set(gca,'Ydir','reverse')
set(gca,'FontSize',16)

for i = 1:size(Data.data.SP,1)
    plot((Data.data.SP(i,63:69)),Data.data.pden(i,63:69))
end 

saveas(gcf,[titl,'_RawSA_Pden.jpg'])

hold off


figure(3)
hold on

title([titl,' Raw Dissolved Oxygen vs Potential Density'])
xlabel('Dissolved Oxygen [micro-mol/kg]')
ylabel('Potential Density')

set(gca,'Ydir','reverse')
set(gca,'FontSize',16)

for i = 1:size(Data.data.DO,1)
    plot((Data.data.DO(i,43:69)),Data.data.pden(i,43:69))
end 

saveas(gcf,[titl,'_RawDO_Pden.jpg'])

hold off