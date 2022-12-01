function [ff_simple,Mean_comp,cust_1,cust_2,cust_3,cust_4,RMSE,difference_simple,comps]=Hfitsimple(input_time,input_record,wn)

global wn

% defining time

t = input_time;
O = input_record;

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
ff_simple = harmfit(coefout,t);

Mean_comp = coefout(1);
cust_1 = sqrt((coefout(2)^2+(coefout(length(wn)+2)^2)));
cust_2 = sqrt((coefout(3)^2+(coefout(length(wn)+3)^2)));
cust_3 = sqrt((coefout(4)^2+(coefout(length(wn)+4)^2)));
cust_4 = sqrt((coefout(5)^2+(coefout(length(wn)+5)^2)));

%finding the RMSE
RMSE = sqrt((nansum((ff_simple-O).^2))/length(O))/(nanmean(ff_simple));

difference_simple = abs(input_record(1:length(input_time)) - ff_simple(:,1));

b = length(wn)+1
comps = zeros([1 b])

for g = 1:b
    if g == 1
        comps(g) = coefout(g)
    else
        comps(g) = sqrt((coefout(g)^2+(coefout(length(wn)+g)^2)))
    end
end

end