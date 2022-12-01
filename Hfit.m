function [ff,Nonperiodic_component,Decadal_component_amplitude,Centenial_component_amplitude,Millenial_component_amplitude,Cust_1_component_amplitude,Cust_2_component_amplitude,Cust_3_component_amplitude,difference]=Hfit(window,input_time,input_record,wn)
%% Looping

global wn

for n=1:(length(input_time)-window)
%% defining time

in = 1+(n-1);
out = (window+1)+(n-1);
t = input_time(in:out);
O = input_record(in:out);

% We will fit the function: T=a0+Cn*sin(wn*t)+Dn*cos(wn*t)
% This are the first estimates of the parameters

a0in = 0.1;
Cnin = zeros(1,length(wn))+0.1;
Dnin = zeros(1,length(wn))+0.1;

% coefin contains all the parameters that have to be optimized.
coefin=[a0in Cnin Dnin];

% Use optimized amps and phases to make the best fit function
[coefout]=nlinfit(t,O,@harmfit,coefin);

% Let nlinfit optimize the amplitudes and phaes.
ff(n,:)=harmfit(coefout,t);

Nonperiodic_component(n) = coefout(1);
Decadal_component_amplitude(n)=sqrt((coefout(2)^2+(coefout(length(wn)+2)^2)));
Centenial_component_amplitude(n)=sqrt((coefout(3)^2+(coefout(length(wn)+3)^2)));
Millenial_component_amplitude(n)=sqrt((coefout(4)^2+(coefout(length(wn)+4)^2)));
Cust_1_component_amplitude(n)=sqrt((coefout(5)^2+(coefout(length(wn)+5)^2)));
Cust_2_component_amplitude(n)=sqrt((coefout(6)^2+(coefout(length(wn)+6)^2)));
Cust_3_component_amplitude(n)=sqrt((coefout(7)^2+(coefout(length(wn)+7)^2)));

%finding the RMSE

% RMSE(n,:) = sqrt(sum((ff-O').^2)./length(O));

end

difference = abs(input_record(1:length(input_time)-window) - ff(:,1));

end