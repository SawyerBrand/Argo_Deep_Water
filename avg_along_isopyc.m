function avg_along_isopyc(file)



f = load(file);

for i = 1:16
    T_avg(i,:) = mean(f.data.T(:,i));
    SA_avg(i,:) = mean(f.data.sa(:,i));
    DO_avg(i,:) = mean(f.data.DO(:,i));
    PDEN(i,:) = f.data.pden(:,i);
end

for i = 17
    T_avg(i,:) = NaN;
    SA_avg(i,:) = NaN;
    DO_avg(i,:) = NaN;
    PDEN(i,:) = f.data.pden(:,i);
end

for i = 18:size(f.data.DO,2) 
    T_avg(i,:) = mean(f.data.T(:,i));
    SA_avg(i,:) = mean(f.data.sa(:,i));
    DO_avg(i,:) = mean(f.data.DO(:,i));
    PDEN(i,:) = f.data.pden(:,i);
end

for i = 1:16
    T_std(i,:) = std(f.data.T(:,i));
    SA_std(i,:) = std(f.data.sa(:,i));
    DO_std(i,:) = std(f.data.DO(:,i));
    PDEN(i,:) = f.data.pden(:,i);
end

for i = 6
    T_std(i,:) = NaN;
    SA_std(i,:) = NaN;
    DO_std(i,:) = NaN;
    PDEN(i,:) = f.data.pden(:,i);
end

for i = 7:size(f.data.DO,2) 
    T_std(i,:) = std(f.data.T(:,i));
    SA_std(i,:) = std(f.data.sa(:,i));
    DO_std(i,:) = std(f.data.DO(:,i));
    PDEN(i,:) = f.data.pden(:,i);
end

eyo = (DO_avg'-f.data.DO(6,:));

hold on
set(gca,'ydir','reverse')
plot((DO_avg'-f.data.DO(17,:)),f.data.P(17,:))
plot((DO_avg'-f.data.DO(17,:)+T_std'),f.data.P(17,:),'r-')
plot((DO_avg'-f.data.DO(17,:)-T_std'),f.data.P(17,:),'r-')
hold off

keyboard
end