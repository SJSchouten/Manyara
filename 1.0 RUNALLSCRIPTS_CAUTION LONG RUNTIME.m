function Africafunct(inputmeasurments)
%% Data import

[T17,T18,fig] = readmeasurments(inputmeasurments);

kiliman = 1;
gisp = 1;
solar = 1;
enso = 1;
lakes = 1;

for kiliman = 1
%Kilimanjaro ice core (solar cycles from antarctica)
kili = xlsread('ice.xlsx');
kili_depth = kili(:,1);
kili_age = kili(:,8);              
kili_SIF1 = kili(:,9);
kili_SIF2 = kili(:,10);              
kili_NIF2 = kili(:,11); 
kili_NIF3 = kili(:,12);              
kili_FURT = kili(:,13); 
kili_mix = kili(:,14);
kili_AVG = kili(:,14);
kili_long_age = kili(:,16);
kili_long_NIF2 = kili(:,17);
kili_long_NIF3 = kili(:,18);
kili_long_mix = kili(:,19);
kili_long_dust = kili(:,20);
husc_long_dust = kili(:,21);
% particles KNIF2 short -
par_age = kili(:,24);
par_Cl = kili(:,25);
par_F = kili(:,26);
par_Na = kili(:,27);
par_NO3 = kili(:,28);
par_NH4 = kili(:,29);
par_Ca = kili(:,30);
par_K = kili(:,31);
par_Mg = kili(:,32);
par_SO4 = kili(:,33);
% particles KNIF3 long -
par_long_age = kili(:,37);
par_long_Cl = kili(:,38);
par_long_F = kili(:,39);
par_long_Na = kili(:,40);
par_long_NO3 = kili(:,41);
par_long_NH4 = kili(:,42);
par_long_Ca = kili(:,43);
par_long_K = kili(:,44);
par_long_Mg = kili(:,45);
par_long_SO4 = kili(:,46);
end
for gisp = 1
gisp = xlsread('GISPII_oxygen_isotopes.xls');
gisp_age = gisp(:,1)*1000;
gisp_dO18 = gisp(:,2);
end
for solar = 1
%DELAYGUE (solar cycles from antarctica)
delaygue = xlsread('delaygue2010be10.xls','data');
delaygue_age = -1*(delaygue(:,1)-2000);
delaygue_Be10 = delaygue(:,2);              % Be10 record from ice core
delaygue_Phimix = delaygue(:,3);            % cosmic ray intensity mixed 
delaygue_Philoc = delaygue(:,4);            % cosmic ray intensity local
delaygue_TSI = delaygue(:,5);               % Total solar irradiance

%STEINHILBER (proxy composite solar data)
steinhilber = xlsread('steinhilber2012','data','A23:G451');
steinhilber_age = steinhilber(:,1)+50;
steinhilber_PC = steinhilber(:,2);          % common production rate of radionuclei
steinhilber_PCerr = steinhilber(:,3);       % error band
steinhilber_Phi = steinhilber(:,4);         % cosmic ray intensity
steinhilber_Phierr = steinhilber(:,5);      % error band
steinhilber_TSI = steinhilber(:,6);         % total solar irradiance (solar activity)
steinhilber_TSIerr = steinhilber(:,7);      % error band

%USOSKIN (sunspot reconstruction)
usoskin = importdata('usoskin2014solar.txt','\t');
usoskin_age = (((usoskin.data(:,1))-2000).*-1);
usoskin_SN = usoskin.data(:,2);             % number of sunspots 
usoskin_SN95 = usoskin.data(:,3);           % error band -
usoskin_SN_95 = usoskin.data(:,4);          % error band +
usoskin_Phi = usoskin.data(:,5);            % cosmic ray intensity
usoskin_Phi95 = usoskin.data(:,6);          % error band -
usoskin_Phi_95 = usoskin.data(:,7);         % errorband +
end
for enso = 1
% ENSO data
yandata = xlsread('yan2011soipr.xls','SOIpr');
yan_age = -1*(yandata(:,1)-2000);   %age
yan_SOI = yandata(:,2);             %ENSO index

lidata = xlsread('enso-li2011.xls','ENSO index');
li_age = -1*(lidata(:,1)-2000);     %age
li_SOI = lidata(:,2);               %ENSO index
li_RM = movmean(li_SOI,10);
li_Var = lidata(:,3);               %ENSO index variance
li_yan = zeros(1103,1);
for k=1:1055
    li_yan(k) = li_SOI(k)+yan_SOI(851+k);
end
end
for lakes = 1 
nvic = xlsread('Composite_african_lakes.xlsx','Stager_2005_N-Victoria');
nvic_age = nvic(:,3);
nvic_swd = nvic(:,4);
nvic_cond = nvic(:,5);
vic = xlsread('Composite_african_lakes.xlsx','Berke_2012_Victoria');
vic_age = vic(:,3);
vic_bit = vic(:,6);
vic_c13 = vic(:,10);
tang = xlsread('Composite_african_lakes.xlsx','Tierney_2010_Tanganyika');
tang_temp_age = tang(:,2);
tang_temp = tang(:,4);
tang_temp_errup = tang(:,5);
tang_temp_errlow = tang(:,6);
tang_char_age = tang(:,8);
tang_char = tang(:,9);
tang_bsi_age = tang(:,12);
tang_bsi = tang(:,13);
tang2 = xlsread('Composite_african_lakes.xlsx','Stager_2009_Tanganyika');
tang_cond_age = tang2(:,2);
tang_cond = tang2(:,3);
tang_cond_errhigh = tang_cond+(tang2(:,5));
tang_cond_errlow = tang_cond-(tang2(:,5));
tang_pol_age = tang2(:,8);
tang_aqpol = tang2(:,9);
naiv = xlsread('Composite_african_lakes.xlsx','Verschuren_2000_Naivasha');
naiv_depth = naiv(:,1);
naiv_llv_age = naiv(:,3);
naiv_llv = naiv(:,5);
naiv_cond_age = naiv(:,9);
naiv_cond1 = naiv(:,11); %Diatom inferred
naiv_cond2 = naiv(:,13); %Chironomid inferred
chal = xlsread('Composite_african_lakes.xlsx','Tierney_2011_Challa');
chal_d_age = chal(:,3);
chal_d = chal(:,5);
chal_bit_age = chal(:,10);
chal_bit = chal(:,12);
sacr = xlsread('Composite_african_lakes.xlsx','Konecky_2015_Sacred');
sacr_d_age = sacr(:,3);
sacr_d = sacr(:,4);
sacr_d_errhigh = sacr_d+(sacr(:,5));
sacr_d_errlow = sacr_d-(sacr(:,5));
sacr_c13_age = sacr(:,10);
sacr_c13 = sacr(:,11);
sacr_c13_err = sacr(:,12);
edwa = xlsread('Composite_african_lakes.xlsx','Russel_2005_2003_Edward');
edwa_age = edwa(:,2);
edwa_mg = edwa(:,3);
edwa_age2_mg = edwa(:,6);
edwa_mg2 = edwa(:,7);
edwa_detrend = edwa(:,8);
ugan = xlsread('Composite_african_lakes.xlsx','Mills_2013_Ugandacrater');
ugan_depth = ugan(:,1);
ugan_age = ugan(:,3); % Lake Nyamogusingiri
ugan_llv = ugan(:,4);
ugan_cond = ugan(:,7);
ugan_depth2 = ugan(:,11);
ugan_age2 = ugan(:,13); % Lake Kyasanduka
ugan_llv2 = ugan(:,14);
ugan_cond2 = ugan(:,17);
bogo = xlsread('Composite_african_lakes.xlsx','Verschuren_2018_Bogoria');
bogo_age = bogo(:,1); % Lake Nyamogusingiri
bogo_llvtop = bogo(:,2);
bogo_llvbot = bogo(:,3);
bogo_ll = bogo(:,4);
MCEOF = xlsread('Composite_african_lakes.xlsx','Tierney_2013_MCEOF1');
MCEOF_age = MCEOF(:,2); % Lake Nyamogusingiri
MCEOF_mean = MCEOF(:,3); 
MCEOF_high = MCEOF(:,4); % 95 confidence
MCEOF_low = MCEOF(:,7); % 95 confidence
maso = xlsread('Composite_african_lakes.xlsx','Garcin_2006_Masoko');
maso_age = maso(:,3); 
maso_mag = maso(:,4); 
end

figure,yyaxis left,subplot(2,1,1),plot(T17(:,3)+27,T17(:,8))
hold on, yyaxis right,plot(naiv_depth,naiv_llv)
% hold on, yyaxis right,plot(MCEOF_age,MCEOF_mean)
% hold on, yyaxis right,plot(ugan_depth,ugan_llv)
hold on, yyaxis right,plot(ugan_depth2,ugan_llv2)
xlim([1,600]);

subplot(2,1,2),plot(T17(:,3)+17,T17(:,8))
hold on,plot(naiv_llv_age,naiv_llv)
xlim([1,600])

%% Resampling the kilimanjaro ice core (staircase data)

for i = 1:length(kili_age)
kili_age2(i+(i-1)+1) = kili_age(i);
kili_age2((2*i)+1) = kili_age(i);
kili_mix2(i+(i-1)) = kili_mix(i);
kili_mix2(2*i) = kili_mix(i);
kili_age2(1) = 0;
kili_mix2(1059) = NaN;
end
for m = 1:length(kili_long_age)
kili_long_age2(m+(m-1)+1) = kili_long_age(m);
kili_long_age2((2*m)+1) = kili_long_age(m);
kili_long_mix2(m+(m-1)) = kili_long_mix(m);
kili_long_mix2(2*m) = kili_long_mix(m);
kili_long_age2(1) = 0;
kili_long_mix2(1059) = NaN;
end
for c = 1:length(bogo_age)
    bogo_age_res(c+(c-1)+1) = bogo_age(c);
    bogo_age_res((2*c)+1) = bogo_age(c);
    bogo_llv_top(c+(c-1)) = bogo_llvtop(c);
    bogo_llv_top(2*c) = bogo_llvtop(c);
    bogo_llv(c+(c-1)) = bogo_ll(c);
    bogo_llv(2*c) = bogo_ll(c);
    bogo_llv_bot(c+(c-1)) = bogo_llvbot(c);
    bogo_llv_bot(2*c) = bogo_llvbot(c);
    bogo_age_res(1) = 0;
    bogo_llv(107) = NaN;
    bogo_llv_top(107) = NaN;
    bogo_llv_bot(107) = NaN;
end
for a = 1:length(chal_d_age)
    chal_d_age_res(a+(a-1)+1) = chal_d_age(a);
    chal_d_age_res((2*a)+1) = chal_d_age(a);
    chal_d_res(a+(a-1)) = chal_d(a);
    chal_d_res(2*a) = chal_d(a);
    chal_d_age_res(1) = 0;
    chal_d_res(403) = NaN;
    chal_bit_age_res(a+(a-1)+1) = chal_bit_age(a);
    chal_bit_age_res((2*a)+1) = chal_bit_age(a);
    chal_bit_res(a+(a-1)) = chal_bit(a);
    chal_bit_res(2*a) = chal_bit(a);
    chal_bit_age_res(1) = 0;
    chal_bit_res(403) = NaN;
end

%% Resampling & periodicity plots

% Resampling
[kili_resampled,kili_age_resampled] = resample(kili_long_mix(1:234),kili_long_age(1:234),1,3,1);
[edwa_resampled,edwa_age_resampled] = resample(edwa_mg,edwa_age,1,3,1);
[chal_resampled,chal_age_resampled] = resample(chal_d,chal_d_age,1,3,1);

% Periodicity
kili_mix_p = kili_resampled - nanmean(kili_resampled);
[pxx,f] = periodogram(kili_resampled,[],[],1);
f1=1./f;

edwa_resampled_p = edwa_resampled - mean(edwa_resampled);
[pxx1,f2] = periodogram(edwa_resampled_p,[],[],1);
f3=1./f2;

chal_resampled_p = chal_resampled - mean(chal_resampled);
[pxx2,f4] = periodogram(chal_resampled_p,[],[],1);
f5=1./f4;

figure,
subplot(2,1,1),plot(kili_age_resampled,kili_resampled,'xg'),hold on,plot(kili_long_age,kili_long_mix,'.-k')
hold on, yyaxis right, plot(edwa_age_resampled,edwa_resampled,'xb'), hold on, plot(edwa_age,edwa_mg,'.-k')
hold on, yyaxis left, plot(chal_age_resampled,chal_resampled+90,'xr'), hold on, plot(chal_d_age,chal_d+90,'.-k')
legend('\delta^{18}O_{kili}','Mg_{Edward}','\deltaD_{Chala}');ylim([-30,0])
grid
subplot(2,1,2),plot(f1,10*log10(pxx),'-xk')
ax = gca;
ax.XLim = [5 2000];
hold on, plot(f3,10*log10(pxx1),'-xb')
hold on, plot(f5,10*log10(pxx2),'-xr')
title('Frequency of cyclicity analysis on both \delta^{18}O records')
xlabel('years/cycles')
ylabel('Amplification magnitude')
legend('\delta^{18}O_{kili}','Mg_{Edward}','\deltaD_{Chala}')
grid;ylim([0,50])

%% Harmonic fitting section 

% Initiating record data
input_time = kili_age(1:150)';
input_record = kili_mix(1:150)';

% Initiating climatic harmonics

global wn

wn(1)=2*pi/50;
wn(2)=2*pi/100;
wn(3)=2*pi/1000;
wn(4)=2*pi/500;
wn(5)=2*pi/200;
wn(6)=2*pi/400;
wn(7)=2*pi/300;
wn(8)=2*pi/800;

% We will fit the function: T=a0+Cn*sin(wn*t)+Dn*cos(wn*t) by the whole dataset

[ff_simple,Mean_comp,Decadal,Centenial,Millenial,Semi_Decadal,difference_simple]=Hfitsimple(input_time,input_record);

% We will fit the function: T=a0+Cn*sin(wn*t)+Dn*cos(wn*t) using fitting windows

window = 10
range = length(input_time)

[ff,Nonperiodic_component,Decadal_component_amplitude,Centenial_component_amplitude,Millenial_component_amplitude,Cust_1_component_amplitude,difference]=Hfit(window,input_time,input_record);

figure,
% Window fit plot
subplot(3,1,1), plot(input_time(1:(range-window)),ff(:,1),'r'),hold on, plot(input_time((range-window):range),ff(50,:),'r')
hold on, plot(input_time,input_record,':k'),grid;
% Amplitude changes over time
subplot(3,1,2),plot(input_time(1:(range-window)),Decadal_component_amplitude','r');
hold on, plot(input_time(1:(range-window)),Centenial_component_amplitude','k');
hold on, plot(input_time(1:(range-window)),Millenial_component_amplitude','b');
hold on, plot(input_time(1:(range-window)),Nonperiodic_component','g');
grid; legend('Dec','Cent','Mil','ext');
% Whole dataset fit
subplot(3,1,3), plot(input_time',ff_simple,'r'),hold on, plot(input_time,input_record,':k'),grid;

%% --------------------------- NICE PLOTS -------------------------------
% % Glacial
% createfigure6(gisp_age,gisp_dO18,kili_age2,kili_mix2,kili_long_age2(61:1059),kili_long_mix2(61:1059),edwa_age,edwa_mg,chal_d_age_res,chal_d_res,chal_bit_age_res,chal_bit_res,vic_age,vic_bit);
% 
% % Millenial I
% sacrd = [sacr_d,sacr_d_errhigh,sacr_d_errlow];
% bogollv = [bogo_ll,bogo_llvbot,bogo_llvtop];
% createfigure5(kili_age2,kili_mix2,kili_age,kili_SIF1,kili_age,kili_SIF2,kili_age,kili_NIF2,kili_age,kili_NIF3,kili_age,kili_FURT, sacr_d_age,sacrd,chal_d_age,chal_d,ugan_age,ugan_llv,ugan_age2,ugan_llv2, naiv_llv_age,naiv_llv,bogo_age,bogollv,MCEOF_age,MCEOF_mean,MCEOF_age,MCEOF_high,MCEOF_age,MCEOF_low);
% 
% % Millenial II
% createfigure3(tang_temp_age,tang_temp,tang_temp_age,tang_temp_errlow,tang_temp_age,tang_temp_errup,naiv_cond_age,naiv_cond2,naiv_cond1,ugan_age2,ugan_cond2,ugan_age,ugan_cond,tang_cond_age,tang_cond,nvic_age,nvic_cond,edwa_age,edwa_mg,par_age,par_Ca,par_age,par_Cl,par_age,par_F,par_age,par_Na,par_age,par_Mg,par_age,par_K);
% 
% % Millenial III
% steintsi = [(steinhilber_TSI+1365.57),((steinhilber_TSI-(0.5*steinhilber_TSIerr))+1365.57),((steinhilber_TSI+(0.5*steinhilber_TSIerr))+1365.57)];
% steinphi = [steinhilber_Phi,(steinhilber_Phi-(0.5*steinhilber_Phierr)),(steinhilber_Phi+(0.5*steinhilber_Phierr))];
% uso = [usoskin_Phi,usoskin_Phi_95,usoskin_Phi95];
% createfigure4(tang_pol_age,tang_aqpol,tang_char_age,tang_char,yan_age,yan_SOI,li_age,li_yan,steinhilber_age,steintsi,delaygue_age,delaygue_Phimix,steinphi,usoskin_age,uso,gisp_age,gisp_dO18);

%% --------------------------- RAW plots -------------------------------
% ---------------------------------------------------------------------
% Glacial FIG
% ---------------------------------------------------------------------

figure, 
subplot(6,1,1), plot(gisp_age,gisp_dO18,'r');ylim([-43,-33])
xlim([1,25000]);title('Glacial proxies compared');ylabel('\deltaO_{18} GISPII');
xticks(0:2000:25000);

subplot(6,1,2),plot(kili_age2,kili_mix2,'k-');ylim([-15,-4]);ylabel('\deltaO_{18} Kilimanjaro');
hold on, plot(kili_long_age2(61:1059),kili_long_mix2(61:1059),'k'); 
xlim([1,25000]);xticks(0:2000:25000);
ax = gca; ax.XAxis.Visible = 'off';

subplot(6,1,3), plot(edwa_age,edwa_mg,'r'), axis ij, ylabel('Lake Edward Mg in calcium (%)');
xlim([1,25000]); xlabel('Time (B2k)');xticks(0:2000:25000);
ax = gca; ax.XAxis.Visible = 'off';

subplot(6,1,4), plot(chal_d_age_res,chal_d_res,'k'); 
axis ij; ylabel('Lake Challa \deltaD_{wax}');
xlim([1,25000]);
ax = gca; ax.XAxis.Visible = 'off';xticks(0:2000:25000);

subplot(6,1,5), plot(chal_bit_age_res,chal_bit_res,'b');ylabel('lake ChallaBIT-index');
xlim([1,25000]);
ax = gca; ax.XAxis.Visible = 'off';xticks(0:2000:25000);

subplot(6,1,6), plot(vic_age,vic_bit,'r');
xlim([1,25000]); ylabel('lake Victoria BIT-index');
xticks(0:2000:25000);

% ------------------------------------------------------------------------
% Millenial FIG part 1
% ------------------------------------------------------------------------

figure, 
% isotopes
subplot(8,1,1), plot(kili_age2,kili_mix2,'k'); 
hold on, plot(kili_age,kili_SIF1,'r:')
hold on, plot(kili_age,kili_SIF2,'g:')
hold on, plot(kili_age,kili_NIF2,'b:')
hold on, plot(kili_age,kili_NIF3,'m:')
hold on, plot(kili_age,kili_FURT,'c:')
set(gca, 'Ydir', 'reverse');ylabel('Kilimanjaro ice core \deltaO_{18}')
xlim([0,2000]);xticks(0:200:2000);
legend('Kilimanjaro \deltaO_{18} mixed','SIF1','SIF2','NIF2','NIF3','FURT','Lake level naivasha')
xlabel('Time (B2k)'); title('Comparisson between millenial proxies')

subplot(8,1,2), plot(sacr_d_age,sacr_d,'k'); 
hold on, plot(sacr_d_age, sacr_d_errhigh,'k:');
hold on, plot(sacr_d_age, sacr_d_errlow,'k:'); ylabel('Sacred lake Deuterium_{wax}');axis ij;
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,3),plot(chal_d_age,chal_d,'b'); axis ij;
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('lake Challa Deuterium_{wax}');

% lake levels
subplot(8,1,4),plot(ugan_age,ugan_llv,'b');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('Nyamogusingiri Lake level (m)')

subplot(8,1,5),plot(ugan_age2,ugan_llv2,'r');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('Kyasanduka Lake level (m)')

subplot(8,1,6),plot(naiv_llv_age, naiv_llv,'k');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('Naivasha Lake level (m)')

subplot(8,1,7),plot(bogo_age,bogo_ll,'g-');
hold on, plot(bogo_age,bogo_llvbot,'g:');
hold on, plot(bogo_age,bogo_llvtop,'g:');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('Bogoria Lake level (m)')

subplot(8,1,8),plot(MCEOF_age,MCEOF_mean,'r');
hold on, plot(MCEOF_age,MCEOF_high,'r:');
hold on, plot(MCEOF_age,MCEOF_low,'r:');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('MCEOF Weighted lake levels')

% ---------------------------------------------------------------------
% Millenial FIG part II
% ---------------------------------------------------------------------

figure,
% SST
subplot(8,1,1),plot(tang_temp_age,tang_temp,'r')
hold on,plot(tang_temp_age,tang_temp_errlow,'r:')
hold on,plot(tang_temp_age,tang_temp_errup,'r:')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('Lake Tanganyika surface temperature (celcius)')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

% Seccondary proxies
subplot(8,1,2),plot(naiv_cond_age,naiv_cond2,'k-');
yyaxis right, hold on, plot(naiv_cond_age,naiv_cond1,'k--');
ylabel('Naivasha Conductivity (uS/cm)')
legend('Chironomid inferred','Diatom inferred')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,3),plot(ugan_age2,ugan_cond2,'r-');
ylabel('Nyamogusingiri Conductivity (uS/cm)')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,4),plot(ugan_age,ugan_cond,'b-');
ylabel('Kyasanduka Conductivity (uS/cm)')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,5),plot(tang_cond_age,tang_cond,'g-');
ylabel('Tanganyika Conductivity (uS/cm)')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,6),plot(nvic_age,nvic_cond,'c-');
ylabel('North-Cictoria Conductivity (uS/cm)');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,7),plot(edwa_age,edwa_mg,'k');
ylabel('Lake Edward Mg in calcium (%)')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,8),plot(par_age,par_Ca,'k'); 
hold on, plot(par_age,par_Cl); 
hold on, plot(par_age,par_F); 
hold on, plot(par_age,par_Na);
hold on, plot(par_age,par_Mg); 
hold on, plot(par_age,par_K);
legend('Ca','Cl','F','Na','Mg','K'); ylabel('Particles in Kilimanjaro ice core (ppb)');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

% -------------------------------------------------------------------
% Millenial FIG part III
% -------------------------------------------------------------------

figure,
subplot(8,1,1), plot(tang_pol_age,tang_aqpol,'k'); ylabel('Tanganyika aquatic Pollen count');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,2), plot(tang_char_age,tang_char,'b'); ylabel('Tanganyika charcoal count');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

% Causes
subplot(8,1,3),plot(yan_age,yan_SOI,'k');
hold on,plot(li_age,li_yan,'r:');
ylabel('SO_{index}');ylim([-3,3])
legend('ENSO from \deltaO_{18}','ENSO from droughts USA');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);

subplot(8,1,4),plot(steinhilber_age,steinhilber_TSI+1365.57,'k')
hold on,plot(steinhilber_age,(steinhilber_TSI-(0.5*steinhilber_TSIerr))+1365.57,'k:')
hold on,plot(steinhilber_age,(steinhilber_TSI+(0.5*steinhilber_TSIerr))+1365.57,'k:')
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('TSI (Total solar irradiance) (W/m^2)')

subplot(8,1,5),plot(delaygue_age,delaygue_Phimix,'g');
hold on,plot(steinhilber_age,steinhilber_Phi,'b');
hold on,plot(steinhilber_age,(steinhilber_Phi-(0.5*steinhilber_Phierr)),'b:');
hold on,plot(steinhilber_age,(steinhilber_Phi+(0.5*steinhilber_Phierr)),'b:');
hold on,plot(usoskin_age,usoskin_Phi,'r');
hold on,plot(usoskin_age,usoskin_Phi_95,'r:');
hold on,plot(usoskin_age,usoskin_Phi95,'r:');
xlim([0,2000]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:200:2000);
ylabel('\Phi (Cosmic ray intensity) (MV)');
legend('\Phi_{antarctica}','\Phi_{composite}','\Phi_{95-error band}','\Phi_{95-error band}','\Phi_{usoskin}','\Phi_{95-error band}','\Phi_{95-error band}');

subplot(8,1,6); plot(gisp_age,gisp_dO18,'k'); 
xlim([1,2000]);xlabel('Time (B2k)');xticks(0:200:2000);
ylabel('\deltaO_{18} GISPII');

% --------------------------
% END
% --------------------------
%% -------------------------- VERY RAW --------------------------
figure,
subplot(2,1,1),plot(bogo_age_res,bogo_llv,'k'); 
hold on, plot(bogo_age_res,bogo_llv_top,'k:'); hold on, plot(bogo_age_res,bogo_llv_bot,'k:'); ylabel('lake level (m)')
% hold on, plot(bogo_age,bogo_ll,'r');hold on, plot(bogo_age,bogo_llvbot,'r:');hold on, plot(bogo_age,bogo_llvtop,'r:'); 
yyaxis right, hold on, plot(MCEOF_age,MCEOF_mean,'r');
hold on, plot(MCEOF_age,MCEOF_high,'r:');
hold on, plot(MCEOF_age,MCEOF_low,'r:');
hold on, plot(yan_age,yan_SOI,'b-'); ylim([-3,4]); xlabel('Time (B2k)'); title('SOI and lake proxies')
grid;legend('Lake Bogoria','Bogoria range','','MCEOF','95 Confidence','','Southern oscilation index');ylabel('Master curve and SOI');

subplot(2,1,2),plot(edwa_age,edwa_mg,'k');grid;xlim([0,2000]);ylim([14,23]);ylabel('%Mg in Calcite')
yyaxis right, plot(maso_age,maso_mag);xlabel('Time (B2k)');ylabel('Magnetic susceptebility');title('Magnetic suscepts')
legend('Lake Maso','Lake Edward')
% ---------------------------------------------------------------------------
figure,
subplot(3,1,1), plot(ugan_age,ugan_llv,'b'), hold on, plot(ugan_age2,ugan_llv2,'r');
grid; xlim([0,2000]); ylabel('Lake level (m) Nyamogusingiri & Kyasanduka')
yyaxis right, hold on, plot(naiv_llv_age, naiv_llv,'k');
hold on, plot(bogo_age,bogo_ll,'g-');hold on, plot(bogo_age,bogo_llvbot,'g:');hold on, plot(bogo_age,bogo_llvtop,'g:');
legend('Lake Nyamogusingiri','Lake Kyasanduka','Lake Naivasha','Lake Bogoria','Lake Bogoria margins')
title('Lake levels compared');xlabel('Time (B2k)');ylabel('Lake level (m) Naivasha & Bogoria')

subplot(3,1,2), plot(naiv_cond_age,naiv_cond2,'k-');
% hold on, plot(ugan_age2,ugan_cond2,'r-');
% hold on,plot(ugan_age,ugan_cond,'b-');
hold on, plot(tang_cond_age,tang_cond,'g-');
hold on, plot(nvic_age,nvic_cond,'c-');ylabel('Conductivity (uS/cm)');ylim([0,600])
yyaxis right, hold on, plot(naiv_cond_age,naiv_cond1,'k--'); grid; xlim([0,2000])
legend('Lake Naivasha','Lake Tangyayika','North victoria','Lake Naivasha')
title('Conductivities compared');xlabel('Time (B2k)');ylabel('Conductivity (uS/cm)')
% 'Lake Nyamogusingiri','Lake Kyasanduka'
subplot(3,1,3),plot(edwa_age,edwa_mg,'k');grid;xlim([0,2000]);ylim([14,23])
xlabel('Time (B2k)');ylabel('%Mg in Calcite');title('Lake Edward records');
%----------------------------------------------------------------------------
figure,
subplot(3,1,1),plot(yan_age,yan_SOI,'k');
hold on,plot(li_age,li_yan,'r:')
xlim([1,2000]); title('ENSO cycles & ice core');xlabel('Time (B2k)');ylabel('SO_{index}');ylim([-3,3])
hold on, yyaxis right, plot(kili_long_age2(61:1059),kili_long_mix2(61:1059),'b-'); 
hold on, plot(kili_age2,kili_mix2,'b-'); ylabel('\deltaO_{18}'); ylim([-13,-7])
grid; legend('ENSO from \deltaO_{18}','ENSO from droughts USA','\deltaO_{18} Kilimanjaro'); title('Kilimanjaro ice core and ENSO signals')

subplot(3,1,2),plot(steinhilber_age,steinhilber_TSI+1365.57,'k')
hold on,plot(steinhilber_age,(steinhilber_TSI-(0.5*steinhilber_TSIerr))+1365.57,'k:')
hold on,plot(steinhilber_age,(steinhilber_TSI+(0.5*steinhilber_TSIerr))+1365.57,'k:')
ylabel('TSI (Total solar irradiance) (W/m^2)')
hold on, yyaxis right, plot(kili_long_age2(61:1059),kili_long_mix2(61:1059),'b'); 
hold on, plot(kili_age2,kili_mix2,'b-'); ylabel('\deltaO_{18}'); ylim([-13,-7])
xlim([1,2000]);xlabel('Time (B2k)'); title('Solar insolation & Kilimanjaro ice core')
grid; legend('TSI_{composite}','TSI_{95-error band}','','\deltaO_{18} Kilimanjaro')

subplot(3,1,3),plot(yan_age,yan_SOI,'k');
% hold on,plot(li_age,li_yan,'k:')
xlim([0,2000]); title('ENSO cycles & lakes');xlabel('Time (B2k)');ylabel('SO_{index}');ylim([-3,3])
yyaxis right, plot(ugan_age,ugan_llv,'b-'), hold on, plot(ugan_age2,ugan_llv2,'r-');
hold on, plot(naiv_llv_age, naiv_llv,'g-'); 
ylabel('Lake level (m)'); title('lake levels & ENSO signal')
grid; legend('ENSO from \deltaO_{18}','Lake Nyamogusingiri','Lake Kyasanduka','Lake Naivasha');
% -------------------------------------------------------------------------
% Composite fig
figure,subplot(3,1,1),plot(delaygue_age,delaygue_Phimix,'g')
hold on,plot(steinhilber_age,steinhilber_Phi,'b')
hold on,plot(steinhilber_age,(steinhilber_Phi-(0.5*steinhilber_Phierr)),'b:')
hold on,plot(steinhilber_age,(steinhilber_Phi+(0.5*steinhilber_Phierr)),'b:')
hold on,plot(usoskin_age,usoskin_Phi,'r')
hold on,plot(usoskin_age,usoskin_Phi_95,'r:')
hold on,plot(usoskin_age,usoskin_Phi95,'r:')
xlim([1,2000]);xlabel('Time (B2k)');ylabel('\Phi (Cosmic ray intensity) (MV)')
hold on, yyaxis right, plot(kili_long_age2(61:1059),kili_long_mix2(61:1059),'k'); 
hold on, plot(kili_age2,kili_mix2,'k-'); ylabel('\deltaO_{18}'); title('Solar proxies and Kilimanjaro ice')
grid; legend('\Phi_{antarctica}','\Phi_{composite}','\Phi_{95-error band}','\Phi_{95-error band}','\Phi_{usoskin}','\Phi_{95-error band}','\Phi_{95-error band}','\deltaO_{18} Kilimanjaro')

subplot(3,1,2),plot(par_age,par_Ca,'k'); 
hold on, plot(par_age,par_Cl); 
hold on, plot(par_age,par_F); 
hold on, plot(par_age,par_Na)
hold on, plot(par_age,par_Mg); 
hold on, plot(par_age,par_K);
legend('Ca','Cl','F','Na','Mg','K'); grid; ylabel('ppb of chemicals');xlim([1,2000]); title('Particles in kilimanjaro'); xlabel('Time(B2k)');

subplot(3,1,3); title('GISPII \deltaO_{18}'); ylabel('\deltaO_{18}')
plot(gisp_age,gisp_dO18,'k'); xlim([1,2000])
hold on, yyaxis right, plot(kili_age2,kili_mix2,'r');ylim([-13,-7]); title('Holocene ice cores compared')
ylabel('\deltaO_{18}');xlabel('Time (B2k)');legend('\deltaO_{18} GISPII','\deltaO_{18} Kilimanjaro.');grid

end
