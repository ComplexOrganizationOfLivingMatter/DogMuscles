%%Extraer Mejor Descriptor PCA, cc de ese descriptor, las 10 mejores cc

clear all
close all



% Asignamos notación para posteriormente cargar cada una de las matrices de
% descriptores

names={'Health','GRMD','GRMD_MIB','GRMD_MIB_38_39_41'};

r=0;
for i=2:length(names)
    
    load(['..\Datos_PCA\PCA_' names{1,1} '_' names{1,i} '_selection_cc_69.mat'])

    matrix{i,1}= [names{1,1} '_' names{1,i}];
    matrix{i,2}=Mejor_pca;
    %matrix{i,3}=indice_cc_seleccionadas;
    
    
    %Nos quedamos con las cc de los 10 mejores descriptores
    matriz_diez_descrip=Mejores(1:10,2:end);
    
    %%Pasamos la matríz a forma de vector y calculamos la frecuencia de
    %%cada cc
    vect_cc=reshape(matriz_diez_descrip,1,size(matriz_diez_descrip,1)*size(matriz_diez_descrip,2));
    [a,b]=hist(vect_cc,1:69);
    
    %%Ordenamos las frecuencias de mayor a menor y nos quedamos con 10
    %%mejores frecuencias
    freq_ord=sort(a,'descend');
    freq_max=freq_ord(1:10);
    freq_max=unique(freq_max);
    freq_max=sort(freq_max,'descend');
    
    
    %Identificamos el número de las cc a las que corresponden las 10 máximas
    %frecuencias
    diez_mejores_cc=[];
    j=1;
    while length(diez_mejores_cc)<10
        cc=find(a==freq_max(j));
        diez_mejores_cc=[diez_mejores_cc,cc];
        j=j+1;
    end

    %matrix{i,10}=sort(diez_mejores_cc(1:10));
    
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],indice_cc_seleccionadas,'Hoja1',['D' num2str(i+3) ':J' num2str(i+3)] );
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],sort(diez_mejores_cc(1:10)),'Hoja1',['K' num2str(i+3)]);

    
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],Mejores_des(1,[2:end])','Hoja1',['C' num2str(11+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{matrix{i,1}},'Hoja1',['B' num2str(10+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'Selection order'},'Hoja1',['C' num2str(10+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],eigenvectors,'Hoja1',['D' num2str(11+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'Projection'},'Hoja1',['D' num2str(9+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'P1','P2'},'Hoja1',['D' num2str(10+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'|P1|+|P2|'},'Hoja1',['F' num2str(10+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'|P1| or |P2|'},'Hoja1',['H' num2str(10+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'ORDER 1'},'Hoja1',['G' num2str(10+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'ORDER 2'},'Hoja1',['I' num2str(10+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{'PCA descriptor'},'Hoja1',['J' num2str(9+(10*r))]);
    xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],{matrix{i,2}},'Hoja1',['J' num2str(10+(10*r))]);

    r=r+1;
end

xlswrite(['..\Datos_PCA\PCA_descriptors_69_cc_' date],matrix,'Hoja1','B4');


