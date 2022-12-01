function [standardized1,standardized2,standardized3,standardized4,standardizedrange,fig,mfig,figstand,Elements2,Depth]=readxrf(inputname)
    % Reading XRF sheet

    [~,~, data] = xlsread(inputname);
    % 'Output_XRF_T17N.xlsx'

    Elements(1,:) = data(1,(20:71));
    Elements(2,:) = data(3,(20:71));
    Elements(3,:) = data(7,(20:71));
    Elements((4:7),:) = data((9:12),(20:71));
    Elements(8,:) = data(40,(20:71));
    Elements((9:23),:) = data((17:31),(20:71));

    Depth(1,:) = data(1,19);
    Depth(2,:) = data(3,19);
    Depth(3,:) = data(7,19);
    Depth((4:7),:) = data((9:12),19);
    Depth((8:23),:) = data((16:31),19);

    Elementsrange = Elements(:,[2,4,6,8,10,12,16,18,26,32,34,52]);
    Elements2 = Elements(:,[1,3,5,7,9,11,15,17,25,31,33,51]);

    for j = 1:12
        standardized((2:23),j) = cell2mat(Elements2((2:23),j))/(max(cell2mat(Elements2((2:23),j))));
        standardizedrange((2:23),j) = cell2mat(Elementsrange((2:23),j))/(max(cell2mat(Elements2((2:23),j))));
    end

    standardized1 = standardized(:,[3,6]);
    standardized2 = standardized(:,[4,2]);
    standardized3 = standardized(:,[5,11]);
    standardized4 = standardized(:,[1,12,7,8,9,10]);

    CaTiratio = cell2mat(Elements2((2:23),5))/cell2mat(Elements2((2:23),6));
    CaTiratio = CaTiratio(:,2);

    inputmeasurments = 'Cond + LOI.xlsx' ;

    [T17,T18,f] = readmeasurments(inputmeasurments);

    figstand = figure,
    for l = 1:12
    subplot(5,3,l),plot(cell2mat(Depth(2:23)),(standardized((2:23),l)),'-xk')
    hold on,plot(cell2mat(Depth(2:23)),((standardized((2:23),l))+(standardizedrange((2:23),l))),':r')
    hold on,plot(cell2mat(Depth(2:23)),((standardized((2:23),l))-(standardizedrange((2:23),l))),':r')
    xlim([0,120]);grid;
    title(Elements2{1,l})
    end
    subplot(5,3,13),plot(T17(:,3),T17(:,8),'-xk');grid;title('Cond')
    subplot(5,3,14),plot(T17(:,3),T17(:,7),'-xk');grid;title('LOI')

    fig = figure,
    for l = 1:12
    subplot(5,3,l),plot(cell2mat(Depth(2:23)),cell2mat(Elements2((2:23),l)),'-xk')
    hold on,plot(cell2mat(Depth(2:23)),(cell2mat(Elements2((2:23),l))+cell2mat(Elementsrange((2:23),l))),':r')
    hold on,plot(cell2mat(Depth(2:23)),(cell2mat(Elements2((2:23),l))-cell2mat(Elementsrange((2:23),l))),':r')
    xlim([0,120]);grid;
    title(Elements2{1,l})
    end
    subplot(5,3,13),plot(T17(:,3),T17(:,8),'-xk');grid;title('Cond')
    subplot(5,3,14),plot(T17(:,3),T17(:,7),'-xk');grid;title('LOI')

    mfig = figure,subplot(2,2,2),plot(cell2mat(Depth(2:23)),standardized1((2:23),1),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized1((2:23),2),'x-')
    grid;xlim([0,120]);legend(Elements2{1,[3,6]});xlabel('Depth(cm)');ylabel('Standardized concentrations');
    subplot(2,2,4),plot(cell2mat(Depth(2:23)),standardized2((2:23),1),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized2((2:23),2),'x-')
    grid;xlim([0,120]);legend(Elements2{1,[4,2]});xlabel('Depth(cm)');ylabel('Standardized concentrations');
    subplot(2,2,3),plot(cell2mat(Depth(2:23)),standardized4((2:23),1),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized4((2:23),2),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized4((2:23),3),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized4((2:23),4),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized4((2:23),5),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized4((2:23),6),'x-')
    grid;xlim([0,120]);legend(Elements2{1,[1,12,7,8,9,10]});xlabel('Depth(cm)');ylabel('Standardized concentrations');
    subplot(2,2,1),plot(cell2mat(Depth(2:23)),standardized3((2:23),1),'x-')
    hold on,plot(cell2mat(Depth(2:23)),standardized3((2:23),2),'x-')
    grid;xlim([0,120]);legend(Elements2{1,[5,11]});xlabel('Depth(cm)');ylabel('Standardized concentrations');

end