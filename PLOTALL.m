% Manyara allplots combi
clear all
close all

% Rainfall data plotting

rmeanwind = 5;
ss = xlsread('sunspots.xlsx','sel');
Rainfallfunct(rmeanwind,ss)

% Coring & lab data plotting
inputfile = 'Corings_example_crosssection_stan';        % specify file name, has to be in same folder as script
cores = (1:18);         % Specify which cores to use, WORKS ONLY IF THE SHEETS ARE DEFINED (in excel) AS: 'Core 1', 'Core 2' etc.

% figure settings
sortdir = 'Ver';       % sorting of cores, Ver == N to S, Hor == E to W
heightset = 0;         % Set an absolute height, if the base level is higher then Sea level
horizonsplot = 'l';      % 0 = no horizons, p = points, l = lines, Plotting interpretation horizons in some way (insert as string)
nrhor = 5;             % nr of interpretation horizons
sandplotting = 0;      % 0 = off, 1 = on; (plotting the grainsize (dotsize in plot) and sorting (as text next to core)
top = 95850; bot = 95350;                                  % specify core plot top and bottom
plotname = 'Corings tanzania lake manyara NP 2019';      % specify plot title

% laboratory settings
lab = 0;               % 0 = off, 1 = on, Plotting labdata asside to cross section
labcores = [17,18];         % specify lab core numbers

Coreplotfunct(inputfile,cores,sortdir,heightset,horizonsplot,nrhor,sandplotting,top,bot,plotname,lab,labcores)

% Paleo record plotting
inputmeasurments = 'Cond + LOI.xlsx' 
Africafunct(inputmeasurments)

