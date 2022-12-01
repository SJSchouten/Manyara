% Script fourier transformers

close all
clear

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
ugan_age = ugan(:,3); % Lake Nyamogusingiri
ugan_llv = ugan(:,4);
ugan_cond = ugan(:,7);
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

% you can change dataset (kilimanjaro, lake eward) as example

% input_time = kili_age(1:150);
% input_record = kili_mix(1:150);
input_time = edwa_age;
input_record = edwa_mg;

y = fft(input_record);
y(1) = [];
figure, plot(y,'ro')
xlabel('real(y)')
ylabel('imag(y)')
title('Fourier Coefficients')

n = length(y);
power = abs(y(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
figure, plot(freq,power)
xlabel('Cycles/Year')
ylabel('Power')

period = 1./freq;
figure, plot(period,power);
xlim([0 1000]); %zoom in on max power
xlabel('Years/Cycle')
ylabel('Power')