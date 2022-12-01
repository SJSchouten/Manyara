clear all
close all

tic

data.names = {'Mbulu','Babati','Karatu','Arusha','Monduli'};

rmeanwind = 5;

ss = xlsread('sunspots.xlsx','sel');

h = waitbar(0,'Running data')

figure,
for n = 1:length(data.names)
    
waitbar(n/length(data.names),h)
    
data.station{n} = xlsread('Monthlyrain.xlsx',string(data.names{n})+'_P');
data.time = data.station{n}(:,1);
data.station{n} = data.station{n}(:,2:13);

    for i = 1:length(data.station{1})-1
        data.dryseason{n}(i,:) = sum(data.station{n}(i,(6:10)));
        data.yearsum{n}(i,:) = sum(data.station{n}(i,(1:12)));
        data.hydryearsum{n}(i,:) = sum(data.station{n}(i,8:12)) + sum(data.station{n}(i+1,(1:7))); %august to august hydrological year
        data.shortrains{n}(i,:) = sum(data.station{n}(i,(11:12))) + data.station{n}(i+1,(1));
        data.longrains{n}(i,:) = sum(data.station{n}(i+1,(2:5)));
        data.wetseason{n}(i,:) = data.station{n}(i,12) + sum(data.station{n}(i+1,(1:5)));
        int.stat{n} = data.station{n}';
        int.range{n} = max(data.station{n})-min(data.station{n});
        int.monthlyanomaly{n}(i,:) = data.station{n}(i,:) - nanmean(data.station{n});
        int.monthlyrelanomaly{n}(i,:) = int.monthlyanomaly{n}(i,:)./int.range{n};
        int.conv{n} = int.monthlyanomaly{n}'; 
    end
    
    int.runmeanH{n} = movmean(data.hydryearsum{n},rmeanwind);
    int.runmeanD{n} = movmean(data.dryseason{n},rmeanwind);
    int.runmeanS{n} = movmean(data.shortrains{n},rmeanwind);
    int.runmeanL{n} = movmean(data.longrains{n},rmeanwind);
    int.anomalyH{n} = data.hydryearsum{n} - nanmean(data.hydryearsum{n});
    int.anomalyD{n} = data.dryseason{n} - nanmean(data.dryseason{n});
    int.anomalyS{n} = data.shortrains{n} - nanmean(data.shortrains{n});
    int.anomalyL{n} = data.longrains{n} - nanmean(data.longrains{n});
    int.anomaly{n} = data.station{n} - nanmean(data.station{n});
    int.time = data.time(1:length(data.station{1})-1);
    
    for b = 1:(12*(length(data.station{1})-1))
        full.station{n}(b) = int.stat{n}(b);
        full.anomaly{n}(b) = int.conv{n}(b);
        full.timem = (1:1:length(full.station{n}));
        full.time = 1922+(full.timem*(100/1200));
        full.plottime = (data.time(1):(100/1200):(max(data.time)-(100/1200)));
    end
    
% if n == 1
%     windowfit = 1
% else
%     windowfit = 0
% end

%% Harmonic fitting section 
% 
input_time = full.timem;
input_record = full.station{n};

ref = find(isnan(input_record));
for k = 1:(length(ref)-1)
if ref(k+1) - ref(k) > 1
    rep(k) = ref(k+1);
    bound = [min(rep(rep > 0)),max(rep)];
end
end

input_record = input_record(1:bound(2));
input_time = input_time(1:bound(2));
full.plottimes{n} = full.plottime(1:bound(2));
full.plottime = full.plottime(1:bound(2));

% Initiating climatic harmonics

global wn

pers{n} = sort([1:3:600]);

for q = 1:length(pers{n})
wn(q)=2*pi/(pers{n}(q));
end

% We will fit the function: T=a0+Cn*sin(wn*t)+Dn*cos(wn*t) by the whole dataset

[simpfitting.ff_simple{n},simpfitting.Mean_comp{n},simpfitting.comp1_amp{n},simpfitting.comp2_amp{n},simpfitting.comp3_amp{n},simpfitting.comp4_amp{n},simpfitting.RMSE{n},simpfitting.diff{n},output.comp{n}]=Hfitsimple(input_time,input_record);
simpfitting.ff_simple{n}(isnan(input_record)) = NaN;

% simpfitting.ff_simple{n}(find(simpfitting.ff_simple{n} < 0)) = 0;

% We will fit the function: T=a0+Cn*sin(wn*t)+Dn*cos(wn*t) using fitting windows

comps{n}(2,:) = output.comp{n}
comps{n}(1,(2:length(pers{n})+1)) = pers{n}

% if windowfit == 1
% 
% window = 5;
% range = length(input_time);
% 
% [fitting.ff{n},fitting.nonper_comp{n},fitting.comp1_amp{n},fitting.comp2_amp{n},fitting.comp3_amp{n},fitting.comp4_amp{n},fitting.comp5_amp{n},fitting.comp6_amp{n},fitting.diff{n}]=Hfit(window,input_time,input_record);
% 
% fitting.ff{n}(isnan(input_record),1) = NaN;
% fitting.nonper_comp{n}(isnan(input_record)) = NaN;
% fitting.comp1_amp{n}(isnan(input_record)) = NaN;
% fitting.comp2_amp{n}(isnan(input_record)) = NaN;
% fitting.comp3_amp{n}(isnan(input_record)) = NaN;
% fitting.comp4_amp{n}(isnan(input_record)) = NaN;
% fitting.comp5_amp{n}(isnan(input_record)) = NaN;
% fitting.comp6_amp{n}(isnan(input_record)) = NaN;
% 
% end
% 
% Whole dataset fit
% 
subplot(((length(data.names)+1)),1,n),plot(full.plottimes{n}',simpfitting.ff_simple{n},':r'),hold on, plot(full.plottimes{n}',input_record,':k');
hold on, plot(full.plottimes{n}',movmean(input_record,12),'k'), hold on, plot(full.plottimes{n}',movmean(simpfitting.ff_simple{n},12),'r');legend('Fit RMSE = '+string(simpfitting.RMSE{n}),'Raw data','Yearly runmean data','Yearly runmean fit')
grid;title(data.names{n});xlim([1922,2018]);ylim([-200 200]);

subplot(((length(data.names)+1)),1,6), plot((ss(:,1)/1000),ss(:,2),'r:'), hold on, plot((ss(:,1)/1000),movmean(ss(:,2),365),'k');
grid;title('Sunspot cycles');ylabel('nr. of sunspots');xlabel('Year');legend('Daily observations','Yearly mean')
xlim([1922,2018]);
% 
% if windowfit == 1
% 
% figure,    
% plot(input_time(1:(range)),fitting.ff{n}(:,1),'r'),hold on, plot(input_time((range-window):range),fitting.ff{n}(50,:),'r')
% hold on,plot(input_time,input_record,':k'),grid;ylim([0,500])
% 
% figure,
% % Amplitude changes over time
% subplot(3,1,1),yyaxis left,plot(input_time(1:(range)),fitting.nonper_comp{n}','k');
% hold on,yyaxis right, plot(input_time(1:(range)),fitting.comp2_amp{n}','b');
% grid; legend('nonper','half yearly');title(string(data.names{n}))
% 
% subplot(3,1,2),yyaxis left,plot(input_time(1:(range)),fitting.comp3_amp{n}','k');
% hold on,yyaxis right, plot(input_time(1:(range)),fitting.comp4_amp{n}','b');
% grid; legend('yearly','five-yearly');
% 
% subplot(3,1,3),yyaxis left,plot(input_time(1:(range)),fitting.comp5_amp{n}','k');
% hold on,yyaxis right, plot(input_time(1:(range)),fitting.comp6_amp{n}','b');
% grid; legend('10-yearly','50-yearly');
% 
% figure,plot(input_time(1:(range)),fitting.nonper_comp{n}');title(string(data.names{n}))
% hold on,plot(input_time(1:(range)),fitting.comp2_amp{n}');
% plot(input_time(1:(range)),fitting.comp3_amp{n}');
% plot(input_time(1:(range)),fitting.comp4_amp{n}');
% plot(input_time(1:(range)),fitting.comp5_amp{n}');
% plot(input_time(1:(range)),fitting.comp6_amp{n}');
% legend('nonper','half yearly','yearly','five-yearly','10-yearly','50-yearly');grid
% 
% end

windowtime = (find(data.time == 2009):1:find(data.time == 2018));

    for t = 1:length(windowtime)-1
        wind.stat{n}(t,:) = data.station{n}(windowtime(t),:);
        wind.conv{n} = wind.stat{n}';
%         wind.fft{n} = simpfitting.ff_simple{n}(windowtime(t),:);
        wind.short{n}(t,:) = data.shortrains{n}(windowtime(t),:);
        wind.long{n}(t,:) = data.longrains{n}(windowtime(t),:);
        wind.dry{n}(t,:) = data.dryseason{n}(windowtime(t),:);
        wind.anomaly{n}(t,:) = int.anomaly{n}(windowtime(t),:);
        wind.anomalyconv{n} = wind.anomaly{n}';
        wind.llv = xlsread('Rdis.xls','Lake LVL');
        wind.llvtime = wind.llv(:,7)'; 
        wind.llv = wind.llv(:,8)';
    end

    for c = 1:(12*(length(windowtime)-1))
        wind.full{n}(c) = wind.conv{n}(c);
        wind.fullan{n}(c) = wind.anomalyconv{n}(c);
        wind.time(c) = c;    
    end

end

close(h)

%% Plotting

figure,
for n = 1:length(data.names)
    loglog(comps{n}(1,:),comps{n}(2,:)); grid; hold on,
    legend(data.names)
end

data.names2 = {'Mbulu','Mbulu fit','Babati','Babati fit','Karatu','Karatu fit','Arusha','Arusha fit','Monduli','Monduli fit'};

figure,
for n = 1:length(data.names)
    subplot(2,1,1),plot(full.time,full.station{n});grid;xlabel('Month nr.');ylabel('Rainfall monthlysum (mm)');title('Rainfall of all stations over time');
    hold on, plot(full.plottimes{n},movmean(simpfitting.ff_simple{n},12),'k');
    hold on, legend(data.names,'fitfunctions')
    hold on,
    subplot(2,1,2),plot(full.time,full.anomaly{n});grid;xlabel('Month nr.');ylabel('Rainfall monthly anomaly (mm)');title('Rainfall anomaly of all stations over time');
    hold on, legend(data.names2)
end

figure,
for n = [2,3,5]
    subplot(2,1,1),plot(wind.llvtime(28:length(wind.llvtime)),wind.full{n});grid;xlabel('Month nr.');ylabel('Rainfall monthlysum (mm)');title('Rainfall of all stations over time');
    hold on, legend(data.names)
    hold on,
end
plot(wind.llvtime,wind.llv,'k-x');
legend(data.names,'Lake level')

data.names2 = {'Mbulu','Mbulu','Babati','Babati','Karatu','Karatu','Arusha','Arusha','Monduli','Monduli'};

figure,
for n = 1:length(data.names)
    subplot(3,1,1),plot(int.time,int.runmeanD{n});grid;xlabel('Year');ylabel('Rainfall dryseasons (mm)');title('Dry');
    hold on, plot(int.time,data.dryseason{n},':'); legend(data.names2)
    hold on,
    subplot(3,1,2),plot(int.time,int.runmeanS{n});grid;xlabel('Year');ylabel('Rainfall shortrains (mm)');title('Short');
    hold on, plot(int.time,data.shortrains{n},':');
    hold on,
    subplot(3,1,3),plot(int.time,int.runmeanL{n});grid;xlabel('Year');ylabel('Rainfall longrains (mm)');title('Long');
    hold on, plot(int.time,data.longrains{n},':');
    hold on,
end

figure,
for n = 1:length(data.names)
    subplot(3,1,1),plot(int.time,int.anomalyD{n});grid;xlabel('Year');ylabel('Rainfall anomaly dryseasons (mm)');title('Dry');
    hold on, legend(data.names);ylim([-500,600])
    hold on,
    subplot(3,1,2),plot(int.time,int.anomalyS{n});grid;xlabel('Year');ylabel('Rainfall anomaly shortrains (mm)');title('Short');
    hold on, legend(data.names);ylim([-500,600])
    hold on,
    subplot(3,1,3),plot(int.time,int.anomalyL{n});grid;xlabel('Year');ylabel('Rainfall anomaly longrains (mm)');title('Long');
    hold on, legend(data.names);ylim([-500,600])
    hold on,
end

figure,
for n = 1:length(data.names)
subplot(2,1,1),plot([1:12],data.station{n},':xk'), hold on,
xticks(0:1:12);xlim([1 12]);xlabel('Month nr.');ylabel('Rainfall(mm)');grid;title('Yearly rainfall distribution all stations');
subplot(2,1,2),plot([1:12],int.anomaly{n}',':xk'), hold on, ;title('Monthly anomaly all stations');grid;xlabel('Month nr.');ylabel('P anomaly (mm)')
xticks(0:1:12);xlim([1 12])
end

toc

