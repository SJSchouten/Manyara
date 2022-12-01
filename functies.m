% functions aanroepen
clear
close all

inputfile = 'Corings_example_crosssection.xlsx';
top = 300; bot = -400;
infolength = 11
nrchar = 12

lab = 'on'
horizonsplot = 'p'
heightset = 0
sortdir = 'Ver'
cores = [7]
labcores = [17,18]

coreplot(inputfile,top,bot,infolength,nrchar,lab,horizonsplot,heightset,sortdir,cores)

inputmeasurments = 'Cond + LOI.xlsx' ;
[T17,T18] = readmeasurments(inputmeasurments);

inputname = 'Output_XRF_T17N.xlsx';
[standardized1,standardized2,standardized3,standardized4,standardizedrange,fig,mfig,figstand,Elements2,Depth] = readxrf(inputname);

% Africafunct(inputmeasurments)