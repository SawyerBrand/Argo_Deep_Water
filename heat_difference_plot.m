function heat_difference_plot(fn,titl)

figure(1)
%shows the difference between the eha 
f = load(fn);
clf

hold on

title([titl,': Differences between Heat Profiles'],'FontSize',20)
xlabel('Delta Q','FontSize',18)
ylabel('Pressure [dbar]','FontSize',18)

set(gca,'ydir','reverse')

for i = 2:1:size(f.idat.heat,1)
    plot((f.idat.heat(i,:)-f.idat.heat(i-1,:)),f.idat.P);
end

%For Float 12700:
%bbb = plot(f.idat.heat(17,:)-f.idat.heat(16,:),f.idat.P,'b*','MarkerSize',2);
%rrr = plot(f.idat.heat(18,:)-f.idat.heat(17,:),f.idat.P,'r*','MarkerSize',2);
%ggg = plot(f.idat.heat(6,:)-f.idat.heat(5,:),f.idat.P,'g*','MarkerSize',2);
%yyy = plot(f.idat.heat(5,:)-f.idat.heat(4,:),f.idat.P,'y*','MarkerSize',2);

%legend('Profile 17 minus 16','LineColor','r*', 'Profile 18 minus 17','Profile 6 minus 5','Profile 7 minus 6','FontSize',18)

%For Float 12881:
%rrr = plot(f.idat.heat(21,:)-f.idat.heat(20,:),f.idat.P,'b*','MarkerSize',2);

%legend('Profile 21 minus 20','FontSize',15)

hold off


figure(2)

%shows the difference between the eha 
f = load(fn);
clf

hold on

title([titl,': Differences between Heat Profiles'],'FontSize',20)
xlabel('Delta Q','FontSize',18)
ylabel('PDen []','FontSize',18)

set(gca,'ydir','reverse')

for i = 2:1:size(f.ipdat.heat,1)
    plot((f.ipdat.heat(i,:)-f.ipdat.heat(i-1,:)),f.ipdat.pden);
end

%For Float 12700:
if titl == 'Float 12700'
    bbb = plot(f.ipdat.heat(17,:)-f.ipdat.heat(16,:),f.ipdat.pden,'b-','MarkerSize',2.5);
    rrr = plot(f.ipdat.heat(18,:)-f.ipdat.heat(17,:),f.ipdat.pden,'r-','MarkerSize',2.5);
    ggg = plot(f.ipdat.heat(6,:)-f.ipdat.heat(5,:),f.ipdat.pden,'g-','MarkerSize',2.5);
    %yyy = plot(f.ipdat.heat(5,:)-f.ipdat.heat(4,:),f.ipdat.P,'y*','MarkerSize',2);
    legend('Profile 17 minus 16','Profile 18 minus 17','Profile 6 minus 5','FontSize',18)
end

%For Float 12881:
if titl == 'Float 12881'
   rrr = plot(f.ipdat.heat(21,:)-f.ipdat.heat(20,:),f.ipdat.pden,'b-','MarkerSize',2.5);
   legend('Profile 21 minus 20','FontSize',15)
end

hold off

end