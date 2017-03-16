%%Llamada extracción de característcas no dapi
clear all
close all

cd ..\..\Matrices_cc

matrices = load ('Matrices_cc.mat');
matrices=struct2cell(matrices);



names={'Health','GRMD','GRMD_MIB'};

cd ..\Codigo
parfor i=2:3

    PCA_2_cc(matrices{1},matrices{i},names{1},names{i})

end

PCA_2_cc(matrices{7},matrices{3},'Health_without_Dubai','GRMD_MIB_similar_age')
PCA_2_cc(matrices{4},[matrices{5};matrices{6}],'Itchy_against','Izmir & Indigo')
PCA_2_cc(matrices{5},[matrices{4};matrices{6}],'Izmir_against','Itchy & Indigo')
PCA_2_cc(matrices{6},[matrices{4};matrices{5}],'Indigo_against','Itchy & Izmir')