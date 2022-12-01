clear
close all

% reading rainfall data for wavelet analysis
inf = 'Monthlyrain.xlsx';

Arushadat = xlsread(inf,'Arusha_P');
Arushadat = Arushadat(:,1:13);
Mbuludat = xlsread(inf,'Mbulu_P');
Mbuludat = Mbuludat(:,1:13);
Mbuludat = xlsread(inf,'Mbulu_P');
Mbuludat = Mbuludat(:,1:13);
Mbuludat = xlsread(inf,'Mbulu_P');
Mbuludat = Mbuludat(:,1:13);
Mbuludat = xlsread(inf,'Mbulu_P');
Mbuludat = Mbuludat(:,1:13);

startMbulu = Mbuludat(1,1);
startArusha = Arushadat(1,1);
k = startMbulu - startArusha;
Arushadata = Arushadat;
Mbuludata((1:k),(1:13)) = NaN;
Mbuludata((1:k),1) = Arushadat((1:k),1);
Mbuludata(((k+1):length(Arushadat)),:) = Mbuludat;

for i = 1:length(Arushadata)-1
    dryseasonsumA(i,:) = sum(Arushadata(i,(7:11)));
    dryseasonsumM(i,:) = sum(Mbuludata(i,(7:11)));
    yearsumA(i,:) = sum(Arushadata(i,(2:13)));
    yearsumM(i,:) = sum(Mbuludata(i,(2:13)));
    HyearsumA(i,:) = sum(Arushadata(i,9:13))+sum(Arushadata(i+1,(2:8))); %august to august hydrological year
    HyearsumM(i,:) = sum(Mbuludata(i,9:13))+sum(Mbuludata(i+1,(2:8)));
    dryseasonsumMix(i,:) = (dryseasonsumA(i,:)+dryseasonsumM(i,:))/2;
    wetseasonsumA(i,:) = Arushadata(i,13)+sum(Arushadata(i+1,(2:6)));
    shortwetseasonsumA(i,:) = Arushadata(i,(12:13))+sum(Arushadata(i+1,(2)));
    longwetseasonsumA(i,:) = sum(Arushadata(i+1,(3:6)));
    wetseasonsumM(i,:) = Mbuludata(i,13)+sum(Mbuludata(i+1,(2:6)));
    wetseasonsumMix(i,:) = (wetseasonsumA(i,:)+wetseasonsumM(i,:))/2;
    if i >= 6
    runningsumM(i,:) = sum(HyearsumM((i-5):i));
    runningsumA(i,:) = sum(HyearsumA((i-5):i));
    runningmeanM(i,:) = mean(HyearsumM((i-5):i));
    runningmeanA(i,:) = mean(HyearsumA((i-5):i));
    runningmeandryM(i,:) = mean(dryseasonsumM((i-5):i));
    runningmeandryA(i,:) = mean(dryseasonsumA((i-2):i));
    runningmeanwetM(i,:) = mean(wetseasonsumM((i-2):i));
    runningmeanwetA(i,:) = mean(wetseasonsumA((i-2):i));
    runningsumMix(i,:) = (runningsumA(i)+runningsumM(i))/2;
    runningmeanMix(i,:) = (runningmeanwetA(i)+runningmeanwetM(i))/2;
    end
end

ind = 1975 - Arushadata(1,1)+1
wind = Arushadata((ind:ind+13),(2:13))'
windM = Mbuludata((ind:ind+13),(2:13))'

anomalyA = Arushadata(:,[2:13])'-nanmean(Arushadata(:,[2:13]))';
anomalyM = Mbuludata(:,[2:13])'-nanmean(Mbuludata(:,[2:13]))';
ManomalyA = Arushadata(:,[2:13])-nanmean(Arushadata(:,[2:13]));
ManomalyM = Mbuludata(:,[2:13])-nanmean(Mbuludata(:,[2:13]));

for h = 1:(12*13)
fullr(h) = wind(h)
fullrM(h) = windM(h)
end
windk = Arushadata(:,(2:13))';
windkM = Mbuludata(:,(2:13))';
for h = 1:(12*92)
fullanomalyA(h) = anomalyA(h);    
fullanomalyM(h) = anomalyM(h);
MfullanomalyA(h) = ManomalyA(h);    
MfullanomalyM(h) = ManomalyM(h);
fullrangeA(h) = windk(h);
fullrangeM(h) = windkM(h);
end
fullrangemnr = [1:1:length(fullrangeA)];

% %% Harmonic fitting section 
% input_time = fullrangemnr(1:870);
% input_record = fullanomalyA(1:870);
% % Initiating record data
% 
% % input_time = Arushadata((1:96),1)'
% % input_record = longwetseasonsumA(1:96)'
% 
% % Initiating climatic harmonics
% 
% global wn
% 
% wn(1)=2*pi/1;
% wn(2)=2*pi/6;
% wn(3)=2*pi/12;
% wn(4)=2*pi/60;
% wn(5)=2*pi/120;
% wn(6)=2*pi/600;
% 
% wn(7)=2*pi/2;
% wn(8)=2*pi/3;
% wn(9)=2*pi/4;
% wn(10)=2*pi/5;
% wn(11)=2*pi/7;
% wn(12)=2*pi/8;
% wn(13)=2*pi/9;
% wn(14)=2*pi/10;
% wn(15)=2*pi/11;
% 
% % wn(16)=2*pi/90;
% % wn(17)=2*pi/900;
% % wn(18)=2*pi/300;
% % wn(19)=2*pi/1200;
% % wn(20)=2*pi/6000;
% 
% % We will fit the function: T=a0+Cn*sin(wn*t)+Dn*cos(wn*t) by the whole dataset
% 
% [ff_simple,Mean_comp,Decadal,Centenial,Millenial,Semi_Decadal,RMSE,difference_simple]=Hfitsimple(input_time,input_record);
% ff_simple(isnan(input_record)) = NaN;
% 
% % We will fit the function: T=a0+Cn*sin(wn*t)+Dn*cos(wn*t) using fitting windows
% 
% window = 5;
% range = length(input_time);
% 
% [ff,Nonperiodic_component,component1_amplitude,component2_amplitude,component3_amplitude,component4_amplitude,component5_amplitude,component6_amplitude,difference]=Hfit(window,input_time,input_record);
% 
% ff(isnan(input_record),1) = NaN;
% Nonperiodic_component(isnan(input_record)) = NaN;
% component1_amplitude(isnan(input_record)) = NaN;
% component2_amplitude(isnan(input_record)) = NaN;
% component3_amplitude(isnan(input_record)) = NaN;
% component4_amplitude(isnan(input_record)) = NaN;
% component5_amplitude(isnan(input_record)) = NaN;
% component6_amplitude(isnan(input_record)) = NaN;
% 
% figure,
% % Window fit plot
% subplot(2,1,1),plot(input_time(1:(range-window)),ff(:,1),'r'),hold on, plot(input_time((range-window):range),ff(50,:),'r')
% hold on,plot(input_time,input_record,':k'),grid;ylim([-500,500])
% 
% % Whole dataset fit
% subplot(2,1,2), plot(input_time',ff_simple,'r'),hold on, plot(input_time,input_record,':k'),grid;legend('RMSE = '+string(RMSE))
% hold on, plot(input_time,runningmeandryA(1:96),'k')
% 
% figure,
% % Amplitude changes over time
% subplot(3,1,1),yyaxis left,plot(input_time(1:(range-window)),Nonperiodic_component','k');
% hold on,yyaxis right, plot(input_time(1:(range-window)),component2_amplitude','b');
% grid; legend('nonper','half yearly');
% 
% subplot(3,1,2),yyaxis left,plot(input_time(1:(range-window)),component3_amplitude','k');
% hold on,yyaxis right, plot(input_time(1:(range-window)),component4_amplitude','b');
% grid; legend('yearly','five-yearly');
% 
% subplot(3,1,3),yyaxis left,plot(input_time(1:(range-window)),component5_amplitude','k');
% hold on,yyaxis right, plot(input_time(1:(range-window)),component6_amplitude','b');
% grid; legend('10-yearly','50-yearly');
% 
% figure,plot(input_time(1:(range-window)),Nonperiodic_component');
% hold on,plot(input_time(1:(range-window)),component2_amplitude');
% plot(input_time(1:(range-window)),component3_amplitude');
% plot(input_time(1:(range-window)),component4_amplitude');
% plot(input_time(1:(range-window)),component5_amplitude');
% plot(input_time(1:(range-window)),component6_amplitude');
% legend('nonper','half yearly','yearly','five-yearly','10-yearly','50-yearly');grid
%% plotting

% Wet and dry season figures
figure,subplot(3,1,1),plot(Arushadata((1:96),1),dryseasonsumA,':k');
hold on, plot(Arushadata((1:96),1),dryseasonsumM,':r');
hold on, plot(Arushadata((1:96),1),runningmeandryA,'k');xlim([1920,1995])
hold on, plot(Arushadata((1:96),1),runningmeandryM,'r');grid;xlabel('year');ylabel('mm');set(gca, 'Xdir', 'reverse');title('Dryseason');legend('Arusha','Mbulu','Run-Mean Arusha','Run-Mean Mbulu')

subplot(3,1,2), plot(Arushadata((1:96),1),wetseasonsumA,':k');
hold on, plot(Arushadata((1:96),1),wetseasonsumM,':r');
hold on, plot(Arushadata((1:96),1),runningmeanwetA,'k');xlim([1920,1995])
hold on, plot(Arushadata((1:96),1),runningmeanwetM,'r');grid;xlabel('year');ylabel('mm');set(gca, 'Xdir', 'reverse');title('Wetseason')

subplot(3,1,3), plot(Arushadata((1:96),1),HyearsumA,':k');
hold on, plot(Arushadata((1:96),1),HyearsumM,':r');
hold on, plot(Arushadata((1:96),1),runningmeanA,'k');xlim([1920,1995])
hold on, plot(Arushadata((1:96),1),runningmeanM,'r');grid;xlabel('year');ylabel('mm');set(gca, 'Xdir', 'reverse');title('Hydrological year')

% Anomaly figures
figure, plot(fullrangemnr,fullanomalyA);
hold on, plot(fullrangemnr,fullanomalyM);ylabel('P anomaly (mm)');xlabel('Year');grid;title('Timeseries of P Anomaly')
xticks(0:48:length(fullanomalyM))

% fullange rainfall
figure, plot(fullrangemnr,fullrangeA);
hold on, plot(fullrangemnr,fullrangeM);;ylabel('P (mm)');xlabel('Year');grid;title('Timeseries of P Rainfall')
xticks(0:48:length(fullanomalyM))

% per month rainfall
figure,plot(Arushadata(:,1),windk)
legend('Month nr. '+string([1:12]));grid;title('Month pracipover time');xlabel('Year');ylabel('rainfall (mm)');

% yearly distribution
figure,subplot(2,1,1),plot([2:13],windk,':xk')
hold on, plot([2:13],windkM,':xr')
xlabel('Month nr.');ylabel('Rainfall(mm)');grid;title('Yearly rainfall distribution');legend('black = arusha', 'red = mbulu')
subplot(2,1,2), plot([2:13],ManomalyA,':xk');
hold on, plot([2:13],ManomalyM,':xr');title('Monthly anomaly');grid;xlabel('Month nr.');ylabel('P anomaly (mm)')
xticks(0:1:12)

% combining with discharge for time frame
table = xlsread('Rdis.xls','Series');
X = (table(:,3)); A = (table(:,4)); B = (table(:,5)); C = (table(:,6)); 
D = (table(:,7)); E = (table(:,8)); F = (table(:,9));

figure,
subplot(6,1,1),plot(X,A,'rx-');grid; legend('MTO WA MBU');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(6,1,2),plot(X,B,'bx-');grid; legend('SIMBA');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(6,1,3),plot(X,C,'gx-'); ylabel('discharge (m^3/s)');grid; legend('KIRURUMO');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(6,1,4),plot(X,D,'mx-');grid; legend('MAKUYUNI');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(6,1,5),plot(X,E,'kx-');grid; legend('');xlabel('Month number(nr.)');xlim([0,156]);xticks(0:12:156)
subplot(6,1,6),plot(X,fullr,'k'); hold on, plot(X,fullrM,'r') ;grid; legend('');xlabel('Month number(nr.)');xlim([0,156]);xticks(0:12:156)
createdischarge(X,A,B,D,E,C)

