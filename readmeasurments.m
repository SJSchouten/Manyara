function [T17,T18,fig] = readmeasurments(inputmeasurments);

[~,~, data] = xlsread(inputmeasurments,'Sheet2');
arraydat = xlsread(inputmeasurments,'Sheet2');
arraydata = arraydat((11:60),:);
cond = arraydata(:,(29:38));
sampleinf = arraydata(:,(1:5));
LOI = arraydata(:,(14:15));

T17(:,(1:5)) = sampleinf((2:26),:);
T17(:,(6:7)) = LOI((2:26),:);
T17(:,8) = cond((2:26),10);

T18(:,(1:5)) = sampleinf((27:42),:);
T18(:,(6:7)) = LOI((27:42),:);
T18(:,8) = cond((27:42),10);

head = data(15,:);
headers = head([1,2,3,4,5,14,15,38]);

fig = figure,
subplot(2,1,1), plot(T18(:,3),T18(:,8),':k');
hold on, plot(T17(:,3),T17(:,8),'rx-');
hold on, plot(T18(:,4),T18(:,8),'kx-');
grid; title('Conductivity T17/T18');legend('T18','T17','T18 revised depth');xlabel('Depth(cm)');ylabel('Conductivity (qS/cm)');

subplot(2,1,2),plot(T18(:,3),T18(:,7),':k');
hold on, plot(T17(:,3),T17(:,7),'rx-');
hold on, plot(T18(:,4),T18(:,7),'kx-');
grid;title('LOI T17/T18');legend('T18','T17','T18 revised depth');xlabel('Depth(cm)');ylabel('Weight percentage organics (%)');

end


