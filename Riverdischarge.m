clear
close all

table = xlsread('Rdis.xls','Series');
X = (table(:,3)); A = (table(:,4)); B = (table(:,5)); C = (table(:,6)); 
D = (table(:,7)); E = (table(:,8)); F = (table(:,9));

figure,
subplot(5,1,1),plot(X,A,'rx-');grid; legend('MTO WA MBU');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(5,1,2),plot(X,B,'bx-');grid; legend('SIMBA');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(5,1,3),plot(X,C,'gx-'); ylabel('discharge (m^3/s)');grid; legend('KIRURUMO');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(5,1,4),plot(X,D,'mx-');grid; legend('MAKUYUNI');xlim([0,156]);ax = gca; ax.XAxis.Visible = 'off';xticks(0:12:156)
subplot(5,1,5),plot(X,E,'kx-');grid; legend('');xlabel('Month number(nr.)');xlim([0,156]);xticks(0:12:156)

createdischarge(X,A,B,D,E,C)
