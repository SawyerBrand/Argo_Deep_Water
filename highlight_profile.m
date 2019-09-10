function highlight_profile(fn,profile)

% This script highlights a specific profile within a whole argo .mat file
% you input, highlighting the profile number you specify 
%for 12700, the weird profile correlates to 33

h = profile
d = load(fn)

figure(1)
hold on
title('Temperature Profiles')

set(gca,'fontsize',18)
set(gca,'ydir','reverse')
for i = (profile-8):2:(profile+8)
    x = plot(d.data.T(i,:),d.data.P(i,:),'k.','markersize',4)
    plot(d.data.T(i,:),d.data.P(i,:),'k-','linewidth',0.5)
end

y = plot(d.data.T(profile,:),d.data.P(profile,:),'r.','markersize',6)
plot(d.data.T(profile,:),d.data.P(profile,:),'r-','linewidth',0.5)

z = plot(d.data.T(profile+2,:),d.data.P(profile+2,:),'c.','markersize',6)
plot(d.data.T(profile+2,:),d.data.P(profile+2,:),'c-','linewidth',0.5)

plot(d.data.T(profile-2,:),d.data.P(profile-2,:),'c.','markersize',6)
plot(d.data.T(profile-2,:),d.data.P(profile-2,:),'c-','linewidth',0.5)

xlabel('Temperature [C]')
ylabel('Pressure [dbar]')
legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

acc_plots

hold off

figure(2)
hold on

title('Absolute Salinity Profiles')

set(gca,'fontsize',18)
set(gca,'ydir','reverse')
for i = (profile-8):2:(profile+8)
    x= plot(d.data.sa(i,:),d.data.P(i,:),'k.','markersize',4)
    plot(d.data.sa(i,:),d.data.P(i,:),'k-','linewidth',0.5)
end

y = plot(d.data.sa(profile,:),d.data.P(profile,:),'r.','markersize',6)
plot(d.data.sa(profile,:),d.data.P(profile,:),'r-','linewidth',0.5)

z = plot(d.data.sa(profile+2,:),d.data.P(profile+2,:),'c.','markersize',6)
plot(d.data.sa(profile+2,:),d.data.P(profile+2,:),'c-','linewidth',0.5)

plot(d.data.sa(profile-2,:),d.data.P(profile-2,:),'c.','markersize',6)
plot(d.data.sa(profile-2,:),d.data.P(profile-2,:),'c-','linewidth',0.5)

ylabel('Pressure [dbar]')
xlabel('Absolute Salinity [g/kg]')
legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

acc_plots

hold off

figure(3)
hold on

title('Dissolved Oxygen Profiles')

set(gca,'fontsize',18)
set(gca,'ydir','reverse')

profile = profile +1; 

for i = (profile-8):2:(profile+8)
    x= plot(d.data.DO(i,:),d.data.P(i,:),'k.','markersize',4)
    plot(d.data.DO(i,:),d.data.P(i,:),'k-','linewidth',0.5)
end

y = plot(d.data.DO(profile,:),d.data.P(profile,:),'r.','markersize',6)
plot(d.data.DO(profile,:),d.data.P(profile,:),'r-','linewidth',0.5)

z = plot(d.data.DO(profile+2,:),d.data.P(profile+2,:),'c.','markersize',6)
plot(d.data.DO(profile+2,:),d.data.P(profile+2,:),'c-','linewidth',0.5)

plot(d.data.DO(profile-2,:),d.data.P(profile-2,:),'c.','markersize',6)
plot(d.data.DO(profile-2,:),d.data.P(profile-2,:),'c-','linewidth',0.5)

ylabel('Pressure [dbar]')
xlabel('Dissolved Oxygen [mol/m^2]')
legend([x z y],'Normal Profiles','3 Before and 3 After Profiles','Weird Profile','location','northwest')

acc_plots

hold off

end
