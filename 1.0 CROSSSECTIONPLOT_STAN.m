%% Core plotting script using for plotting cross sections

% Welcome to this core plotting tool
% Use the excel sheet to insert the coring data
% Specify core plotting properties in section below
% By default it plots the full length of the core as a black line
% -  The program reads texture(colors in plot), organics to plot in cross-section
%(+ = Bulk organics; * = Remains,spots,org. layers)
% -  Redox is plotted in the filled dots 
%(red = oxidized, blue is = oxidized and reduced, full reduction is not plotted)
% -  The text in the "interpretation list" in excel sheet is plotted as text next to core
% -  With heightset, one can set a constant height to cores 
%(by default this is turned off and the cross-section is plotted to Zact in excelsheet)

%% clearing the workspace

clear
close all

%% --- Modify here ---

% input settings
inputfile = 'Corings_example_crosssection_stan';        % Specify file name, has to be in same folder as script
cores = (1:18);         % Specify which cores to use, WORKS ONLY IF THE SHEETS ARE DEFINED (in excel) AS: 'Core 1', 'Core 2' etc.

% figure settings
sortdir = 'Ver';       % sorting of cores, Ver == N to S, Hor == E to W
heightset = 0;         % Set an absolute height, if the base level is higher then Sea level
horizonsplot = 'l';    % 0 = no horizons, p = points, l = lines, Plotting interpretation horizons in some way (insert as string)
nrhor = 5;             % nr of interpretation horizons
sandplotting = 0;      % 0 = off, 1 = on; (plotting the grainsize (dotsize in plot) and sorting (as text next to core)
top = 95850; bot = 95350;                                % specify core plot top and bottom
plotname = 'Corings tanzania lake manyara NP 2019';      % specify plot title

% laboratory settings
lab = 1;                    % 0 = off, 1 = on, Plotting labdata asside to cross section
labcores = [17,18];         % specify lab core numbers

%% --- runs automatically from here ---

numberofcores = length(cores);

% Checking whether the combi's are possible

if lab == 1
for n = 1:length(labcores)
if sum(cores == labcores(n)) == 0
   disp('Core 17 is lab core, lab = on and T17 is not there');
   return
end
end
end

%% Loading the data into a Cell stucture

infolength = nrhor+6;

for i = 1:numberofcores
Sheets(i) = ('Core '+ string(cores(i)));
[~,~, data.core{i}] = xlsread(inputfile,Sheets{i});
for k = 1:infolength
infohold(i,k) = data.core{i}(5,(2+k));
end
end

nrchar = size(data.core{1});
nrchar = nrchar(2);
data.title = 'Corings';

% Combining the Metadata from the cores into one array (i.e. coordinates,depth, accuracy)
info(:,1) = 1:numberofcores;
info(:,2:(infolength+1)) = cell2mat(infohold);
if heightset > 0
    info(:,5) = heightset;                 % Core height settings
end

if sortdir == 'Ver'
    infosort = sortrows(info,3,'descend'); % Sorting this array from N to S
elseif sortdir == 'Hor'
    infosort = sortrows(info,4,'descend'); % Sorting this array from W to E
end

%% Calculating cross sectional lengthpath and relative positioning

positioning = zeros(1,numberofcores);        

for n = 1:(numberofcores-1)
    diff(n) = sqrt(((single(infosort(n,3))-single(infosort(n+1,3)))^2)+((single(infosort(n,2))-single(infosort(n+1,2)))^2));
end % Differences between cores

for u = 1:(numberofcores)
    if u == 1
        positioning(u) = (sum(diff)/12.5);
    else
        positioning(u) = positioning(u-1) + diff(u-1);
    end
end % Net positioning

orderhold(:,1) = (1:1:numberofcores);
orderhold(:,2) = infosort(:,1);
ordertwo = sortrows(orderhold,2);
order = ordertwo(:,1);               % Core plotting order (North to South)
xplotmarge = sum(diff)/12.5;         % Defining the core plotting margins
dist = sum(diff) + (2*xplotmarge);   
positioning(1) = xplotmarge;

%% Reformatting data into cells --> array structure

for m = 1:nrchar
    names{:,m} = data.core{1}{9,m};
    type{:,m} = data.core{1}{11,m};
    for j = 1:numberofcores
        char{m,(j+2)} = data.core{j}(12:length(data.core{j}(:,m)),m); char{m,1} = names{m};char{m,2} = type{m};
        if m == 2
        maxdepthhold(j) = max(cell2mat(char{2,(j+2)}));
        end
    end
end

maxdepth = max(maxdepthhold);

%% Defining the top, bottom and middle depth arrays used for plotting

plotyaxis_top = [0:-10:-(maxdepth-10)];
plotyaxis_bot = [-10:-10:-maxdepth];
plotyaxis_mid = [-5:-10:-maxdepth];

% Reformatting Texture and Organics data into nummerical arrays
% P = 1; CP = 2; PC = 3; C = 4; SiC = 5; Si = 6; SiL = 7
% L = 8; SL = 9; LS = 10; vfS = 11; fS = 12; mS = 13; cS = 14

% Organics, M1 = 1; M2 = 2; Plr or vlk or both = 3; lyr = 4

% Redox, O = 1; OR = 2; R = 3; lyr = 4
% Groundwater, GHG = 1; GW = 2; GLG = 3;

textnum = nan((((maxdepth/10))),(numberofcores));  
orgnum = nan((((maxdepth/10))),(numberofcores));  
rednum = nan((((maxdepth/10))),(numberofcores));  
gwnum = nan((((maxdepth/10))),(numberofcores));  

for h = 1:numberofcores
    for p = 1:length(char{3,(h+2)})
        switch char{3,(h+2)}{p}
            case 'C'
                textnum(p,h) = 4.0;
            case 'P'
                textnum(p,h) = 1.0;
            case 'CP'
                textnum(p,h) = 2.0;
            case 'PC'
                textnum(p,h) = 3.0;     
            case 'SiC'
                textnum(p,h) = 5.0;
            case 'CL'
                textnum(p,h) = 6.0;
            case 'SiCL'
                textnum(p,h) = 6.0;
            case 'Si'
                textnum(p,h) = 7.0;
            case 'SiL'
                textnum(p,h) = 8.0;
            case 'L'
                textnum(p,h) = 9.0;
            case 'SL'
                textnum(p,h) = 10.0;
            case 'LS'
                textnum(p,h) = 11.0;
            case 'vfS'
                textnum(p,h) = 12.0;
            case 'fS'
                textnum(p,h) = 13.0;
            case 'mS'
                textnum(p,h) = 14.0;
            case 'cS'
                textnum(p,h) = 15.0;
        end

        switch char{5,(h+2)}{p}
            case 'vlk'
                orgnum(p,h) = 3;
            case 'plr'
                orgnum(p,h) = 3;
            case 'lyr'
                orgnum(p,h) = 4;
            case 'vlk, plr'
                orgnum(p,h) = 3;
            case 'plr, lyr'
                orgnum(p,h) = 4;
            case 'vlk, lyr'
                orgnum(p,h) = 4;
            case 'vlk, plr, lyr'
                orgnum(p,h) = 4;
            case 'M1'
                orgnum(p,h) = 1;     
            case 'M2'
                orgnum(p,h) = 2;
            case 'M1, vlk'
                orgnum(p,h) = 1;
            case 'M1, plr'
                orgnum(p,h) = 1;
            case 'M1, lyr'
                orgnum(p,h) = 1;
            case 'M2, vlk'
                orgnum(p,h) = 2;
            case 'M2, plr'
                orgnum(p,h) = 2;
            case 'M2, lyr'
                orgnum(p,h) = 2;
            case 'M1, vlk, plr'
                orgnum(p,h) = 1;
            case 'M1, plr, lyr'
                orgnum(p,h) = 1;
            case 'M1, vlk, lyr'
                orgnum(p,h) = 1;
            case 'M1, vlk, plr, lyr'
                orgnum(p,h) = 1;
            case 'M2, vlk, plr'
                orgnum(p,h) = 2;
            case 'M2, plr, lyr'
                orgnum(p,h) = 2;
            case 'M2, vlk, lyr'
                orgnum(p,h) = 2;
            case 'M2, vlk, plr, lyr'
                orgnum(p,h) = 2;
        end
        
        switch char{16,(h+2)}{p}
            case 'o'
                rednum(p,h) = 1;
            case 'or'
                rednum(p,h) = 2;
            case 'r'
                rednum(p,h) = 3;
            case 'O'
                rednum(p,h) = 1;
            case 'O/R'
                rednum(p,h) = 2;
            case 'R'
                rednum(p,h) = 3;
            case 'R/O'
                rednum(p,h) = 2;
            case 'oxidized'
                rednum(p,h) = 1;     
            case 'Oxidized'
                rednum(p,h) = 1;
            case 'reduced'
                rednum(p,h) = 3;
            case 'Reduced'
                rednum(p,h) = 3;
        end
        
        switch char{17,(h+2)}{p}
            case 'GW'
                gwnum(p,h) = 2;
            case 'Gw'
                gwnum(p,h) = 2;
            case 'gw'
                gwnum(p,h) = 2;
            case 'GHG'
                gwnum(p,h) = 1;
            case 'ghg'
                gwnum(p,h) = 1;
            case 'GLG'
                gwnum(p,h) = 3;
            case 'glg'
                gwnum(p,h) = 3;
            case 'Ghg'
                gwnum(p,h) = 1;     
            case 'Glg'
                gwnum(p,h) = 3; 
        end
    end
    
    % Defining core arrays for plotting, merging texture, depth, organics and positioning data
    coreplot{1,h} = orderhold(h,2);
    coreplot{2,h} = zeros((maxdepth/10),4);
    coreplot{2,h}((1:(maxdepth/10)),2) = plotyaxis_mid(1,:);
    coreplot{2,h}((1:(maxdepth/10)),3) = plotyaxis_bot(1,:);
    coreplot{2,h}((1:(maxdepth/10)),4) = textnum(:,h);
    coreplot{2,h}((1:(maxdepth/10)),5) = orgnum(:,h);
    coreplot{2,h}((1:(maxdepth/10)),6) = rednum(:,h);
    coreplot{2,h}((1:(maxdepth/10)),7) = gwnum(:,h);
    coreplot{2,h}((1:(maxdepth/10)),1) = positioning(order(h));
    for t = 1:(maxdepth/10)
        if (-plotyaxis_bot(1,t)) <= (info(h,7))
            coreplot{2,h}(t,2) = plotyaxis_mid(1,t) + (info(h,5)*100);
            coreplot{2,h}(t,3) = plotyaxis_bot(1,t) + (info(h,5)*100);
        else
            coreplot{2,h}(t,2) = NaN;
            coreplot{2,h}(t,3) = NaN;
        end
    end
end

%% Horizon interpretation

horizons(:,1) = positioning';
horizons(:,2) = infosort(:,5);
for p = 1:nrhor
    horizons(:,(p+2)) = infosort(:,5)-(infosort(:,(p+7))*10^-2);
end

labels = 'Core '+string(ordertwo(:,2));

%% Reading function data from the lab

if lab == 1
    inputmeasurments = 'Cond + LOI.xlsx' ;
    [T17,T18] = readmeasurments(inputmeasurments);
    inputname = 'Output_XRF_T17N.xlsx';
    [standardized1,standardized2,standardized3,standardized4,standardizedrange,fig,mfig,figstand,Elements2,Depth] = readxrf(inputname);
    T17Height = cell2mat(data.core{1,(find(cores == 17))}(5,6))*100;
end

%% Plotting

if lab == 1
    
    figure,subplot(1,20,[1:14])   
    for l = 1:numberofcores
    plot(coreplot{2,l}(:,1),coreplot{2,l}(:,2),'k');
    if labels(l) == 'Core '+string(labcores(2)) | labels(l) == 'Core '+string(labcores(1))
    else
        text(coreplot{2,l}(1,1),coreplot{2,l}(1,2)+30,labels(l),'FontSize',8,'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
    
    if sandplotting == 0
        hold on, scatter(coreplot{2,l}(:,1),coreplot{2,l}(:,2),30,coreplot{2,l}(:,4),'square','filled');
    for z = 1:length(coreplot{2,l})
        if z < length(char{1,(l+2)})
        if isnan(char{25,(l+2)}{z}) == 0 
            text(double((coreplot{2,l}(z,1))-(1*(dist/140))),(coreplot{2,l}(z,2)),char{25,(l+2)}(z),'FontSize',6,'HorizontalAlignment','Center','VerticalAlignment','Middle')
        end
        end
    end  
    end
    
    for z = 1:length(coreplot{2,l})
        if sandplotting == 1
            if coreplot{2,l}(z,4) < 9
                hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),30,coreplot{2,l}(z,4),'square','filled');
            else
                if z < length(char{1,(l+2)})
                    if isnan(char{8,(l+2)}{z}) == 0 
                        text(double((coreplot{2,l}(z,1))+(1*(dist/140))),(coreplot{2,l}(z,2)),char{8,(l+2)}(z),'FontSize',8,'HorizontalAlignment','Center','VerticalAlignment','Middle')
                    end
                    if isnan(char{7,(l+2)}{z}) == 0
                        hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),(cell2mat(char{7,(l+2)}(z)))/3,coreplot{2,l}(z,4),'square','filled');
                    else
                        hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),30,coreplot{2,l}(z,4),'square','filled');
                    end
                end
            end
        end
        
        if coreplot{2,l}(z,6) == 1.0
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),7,'O','filled','r');
        elseif coreplot{2,l}(z,6) == 2.0
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),7,'O','filled','b');
        end
        
        if coreplot{2,l}(z,5) == 2 | coreplot{2,l}(z,5) == 1
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),20,'+','k');
        elseif coreplot{2,l}(z,5) == 3 | coreplot{2,l}(z,5) == 4
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),20,'*','k');
        end
       
        if coreplot{2,l}(z,7) == 2.0
            hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),40,'v','k');
        end
    end
    end

    hold on, plot(horizons(:,1),horizons(:,2)*10^2,'k')
    if horizonsplot == 'p'
        hold on, scatter(horizons(:,1),horizons(:,3)*10^2,30,'square','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,4)*10^2,30,'diamond','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,5)*10^2,30,'v','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,6)*10^2,30,'^','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,7)*10^2,30,'O','filled','k')
    elseif horizonsplot == 'l'
        hold on, plot(horizons(:,1),horizons(:,3)*10^2,'r-.x')
        hold on, plot(horizons(:,1),horizons(:,4)*10^2,'k:x')
        hold on, plot(horizons(:,1),horizons(:,5)*10^2,'k--x')
        hold on, plot(horizons(:,1),horizons(:,6)*10^2,'k-.x')
        hold on, plot(horizons(:,1),horizons(:,7)*10^2,'k:x')
    end
    grid; grid minor;
    xlim([0,dist]); ylim([bot,top]);title(plotname)
    ylabel('Height (cm)'); xlabel('Cross-sectional distance (m)')
    c = colorbar;
    c.Label.String = 'Grainsize Fine --> Course';
    map = [0.400000005960464 0.200000002980232 0;0.400000005960464 0.200000002980232 0;0.400000005960464 0.200000002980232 0;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0 1 0;0 1 0;0 1 0;0 1 0;0 1 0;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;1 1 0;1 1 0;1 1 0;1 1 0;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.501960813999176 0;1 0.501960813999176 0;1 0.501960813999176 0;1 0.501960813999176 0;1 0.215126067399979 0;1 0.179271727800369 0;1 0.143417373299599 0;1 0.107563033699989 0;0.600000023841858 0 0;0.600000023841858 0 0;0.600000023841858 0 0];
        colormap(map); 
    colorbar('Ticks',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],...
             'TickLabels',{'P','CP','PC','C','SiC','(Si)CL','Si','SiL','L','SL','LS','vfS','fS','mS','cS'})
    caxis([1 15]);

    subplot(1,20,15),plot(T17(:,8),(-T17(:,4)+T17Height),'x-k');
    hold on,plot(T18(:,8),(-T18(:,4)+T17Height),'x-r');
    grid;ylim([bot,top]);title('Conductivity');xlabel('(mS/cm)')

    subplot(1,20,16),plot(T17(:,7),(-T17(:,4)+T17Height),'-xk');
    hold on,plot(T18(:,7),(-T18(:,4)+T17Height),'-xr');
    legend('17','18');grid;ylim([bot,top]);title('LOI');xlabel('(%)')

    subplot(1,20,17),plot(standardized1((2:23),1),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized1((2:23),2),-cell2mat(Depth(2:23))+T17Height,'x-')
    grid;ylim([bot,top]);legend(Elements2{1,[3,6]})

    subplot(1,20,18),plot(standardized2((2:23),1),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized2((2:23),2),-cell2mat(Depth(2:23))+T17Height,'x-')
    grid;ylim([bot,top]);legend(Elements2{1,[4,2]})

    subplot(1,20,19),plot(standardized4((2:23),1),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized4((2:23),2),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized4((2:23),3),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized4((2:23),4),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized4((2:23),5),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized4((2:23),6),-cell2mat(Depth(2:23))+T17Height,'x-')
    grid;ylim([bot,top]);legend(Elements2{1,[1,12,7,8,9,10]});xlabel('standardized concentrations');title('XRF standardized')

    subplot(1,20,20),plot(standardized3((2:23),1),-cell2mat(Depth(2:23))+T17Height,'x-')
    hold on,plot(standardized3((2:23),2),-cell2mat(Depth(2:23))+T17Height,'x-')
    grid;ylim([bot,top]);legend(Elements2{1,[5,11]})

elseif lab == 0
    figure,subplot(1,20,[1:20])
    
    for l = 1:numberofcores
    plot(coreplot{2,l}(:,1),coreplot{2,l}(:,2),'k');
    if labels(l) == 'Core '+string(labcores(2)) | labels(l) == 'Core '+string(labcores(1))
    else
        text(coreplot{2,l}(1,1),coreplot{2,l}(1,2)+30,labels(l),'FontSize',8,'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    end
    
    if sandplotting == 0
        hold on, scatter(coreplot{2,l}(:,1),coreplot{2,l}(:,2),30,coreplot{2,l}(:,4),'square','filled');
    for z = 1:length(coreplot{2,l})
        if z < length(char{1,(l+2)})
        if isnan(char{25,(l+2)}{z}) == 0 
            text(double((coreplot{2,l}(z,1))-(1*(dist/140))),(coreplot{2,l}(z,2)),char{25,(l+2)}(z),'FontSize',6,'HorizontalAlignment','Center','VerticalAlignment','Middle')
        end
        end
    end
       
    end
    
    for z = 1:length(coreplot{2,l})
        if sandplotting == 1
            if coreplot{2,l}(z,4) < 9
                hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),30,coreplot{2,l}(z,4),'square','filled');
            else
                if z < length(char{1,(l+2)})
                    if isnan(char{8,(l+2)}{z}) == 0 
                        text(double((coreplot{2,l}(z,1))+(1*(dist/140))),(coreplot{2,l}(z,2)),char{8,(l+2)}(z),'FontSize',8,'HorizontalAlignment','Center','VerticalAlignment','Middle')
                    end
                    if isnan(char{7,(l+2)}{z}) == 0
                        hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),(cell2mat(char{7,(l+2)}(z)))/3,coreplot{2,l}(z,4),'square','filled');
                    else
                        hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),30,coreplot{2,l}(z,4),'square','filled');
                    end
                end
            end
        end
        
        if coreplot{2,l}(z,6) == 1.0
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),7,'O','filled','r');
        elseif coreplot{2,l}(z,6) == 2.0
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),7,'O','filled','b');
        end
        
        if coreplot{2,l}(z,5) == 2 | coreplot{2,l}(z,5) == 1
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),20,'+','k');
        elseif coreplot{2,l}(z,5) == 3 | coreplot{2,l}(z,5) == 4
            hold on, scatter(coreplot{2,l}(z,1)+(dist/140),coreplot{2,l}(z,2),20,'*','k');
        end
       
        if coreplot{2,l}(z,7) == 2.0
            hold on, scatter(coreplot{2,l}(z,1),coreplot{2,l}(z,2),40,'v','k');
        end
    end
    end

    hold on, plot(horizons(:,1),horizons(:,2)*10^2,'k')
    if horizonsplot == 'p'
        hold on, scatter(horizons(:,1),horizons(:,3)*10^2,30,'square','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,4)*10^2,30,'diamond','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,5)*10^2,30,'v','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,6)*10^2,30,'^','filled','k')
        hold on, scatter(horizons(:,1),horizons(:,7)*10^2,30,'O','filled','k')
    elseif horizonsplot == 'l'
        hold on, plot(horizons(:,1),horizons(:,3)*10^2,'r-.x')
        hold on, plot(horizons(:,1),horizons(:,4)*10^2,'k:x')
        hold on, plot(horizons(:,1),horizons(:,5)*10^2,'k--x')
        hold on, plot(horizons(:,1),horizons(:,6)*10^2,'k-.x')
        hold on, plot(horizons(:,1),horizons(:,7)*10^2,'k:x')
    end
    
    grid; grid minor;
    xlim([0,dist]); ylim([bot,top]);title(plotname)
    ylabel('Height (cm)'); xlabel('Cross-sectional distance (m)')
    c = colorbar;
    c.Label.String = 'Grainsize Fine --> Course';
    map = [0.400000005960464 0.200000002980232 0;0.400000005960464 0.200000002980232 0;0.400000005960464 0.200000002980232 0;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0 1 0;0 1 0;0 1 0;0 1 0;0 1 0;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;1 1 0;1 1 0;1 1 0;1 1 0;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.501960813999176 0;1 0.501960813999176 0;1 0.501960813999176 0;1 0.501960813999176 0;1 0.215126067399979 0;1 0.179271727800369 0;1 0.143417373299599 0;1 0.107563033699989 0;0.600000023841858 0 0;0.600000023841858 0 0;0.600000023841858 0 0];
        colormap(map); 
    colorbar('Ticks',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],...
             'TickLabels',{'P','CP','PC','C','SiC','(Si)CL','Si','SiL','L','SL','LS','vfS','fS','mS','cS'})
    caxis([1 15]);  
end

%% 3D-ploting

gridplotdat(:,3:7) = horizons(:,2:6)
gridplotdat(:,1:2) = info(:,2:3)

x = gridplotdat(:,1);
y = gridplotdat(:,2);
z = gridplotdat(:,3)*100;
z1 = gridplotdat(:,4)*100;
z2 = gridplotdat(:,5)*100;
z3 = gridplotdat(:,6)*100;
z4 = gridplotdat(:,6)*100;

[xq,yq] = meshgrid((min(x)):10:(max(x)), (min(y):10:(max(y))));
vq = griddata(x,y,z,xq,yq,'natural');
vq1 = griddata(x,y,z1,xq,yq,'natural');
vq2 = griddata(x,y,z2,xq,yq,'natural');
vq3 = griddata(x,y,z3,xq,yq,'natural');
vq4 = griddata(x,y,z4,xq,yq,'natural');

figure,
mesh(xq,yq,vq); hold on, plot3(x,y,z,'x');
hold on, mesh(xq,yq,vq1);hold on,plot3(x,y,z1,'x');
hold on, mesh(xq,yq,vq2);hold on,plot3(x,y,z2,'x');
hold on, mesh(xq,yq,vq3);hold on,plot3(x,y,z3,'x');
hold on, mesh(xq,yq,vq4);hold on,plot3(x,y,z4,'x');
for b = 1:numberofcores
    plotcoord([1:length(coreplot{2,b})],1) = gridplotdat(b,1);
    plotcoord([1:length(coreplot{2,b})],2) = gridplotdat(b,2);
    hold on, plot3(plotcoord(:,1),plotcoord(:,2),coreplot{2,b}(:,2),'-','Color','k');
    hold on, scatter3(plotcoord(:,1),plotcoord(:,2),coreplot{2,b}(:,2),20,coreplot{2,b}(:,4));
        c = colorbar;
    c.Label.String = 'Grainsize Fine --> Course';
    map = [0.400000005960464 0.200000002980232 0;0.400000005960464 0.200000002980232 0;0.400000005960464 0.200000002980232 0;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.600000023841858 0.200000002980232 0.10196078568697;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.560784339904785 0.600000023841858 0.0784313753247261;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.10196078568697 0.301960796117783 0.10196078568697;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0.200000002980232 0.600000023841858 0.0509803928434849;0 1 0;0 1 0;0 1 0;0 1 0;0 1 0;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.490196079015732 0.180392161011696 0.560784339904785;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0 1;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.600000023841858 0.501960813999176 0.200000002980232;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;0.800000011920929 0.800000011920929 0.10196078568697;1 1 0;1 1 0;1 1 0;1 1 0;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.701960802078247 0.200000002980232;1 0.501960813999176 0;1 0.501960813999176 0;1 0.501960813999176 0;1 0.501960813999176 0;1 0.215126067399979 0;1 0.179271727800369 0;1 0.143417373299599 0;1 0.107563033699989 0;0.600000023841858 0 0;0.600000023841858 0 0;0.600000023841858 0 0];
        colormap(map); 
    colorbar('Ticks',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],...
             'TickLabels',{'P','CP','PC','C','SiC','(Si)CL','Si','SiL','L','SL','LS','vfS','fS','mS','cS'})
    caxis([1 15]);
end
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
zlim([95350 95800]);

%% end of this script
