cd ..\Matrices_cc

%% Cargamos todos los índices de las posiciones de los perros y todas las matrices
load('Indices_total.mat');
load('Indices_musc.mat');
load('Indices_musc_edad.mat');
load('Matrices_cc.mat')
n_img_tipo1=size(Health_cc,1);


%% Representacion del PCA Health VS GRMD global

cd ..\Datos_PCA

load ('PCA_Health_GRMD_selection_cc_69.mat')

Proyecc=Proy{1,1};
h=figure; p1=plot(Proyecc(1,1:n_img_tipo1),Proyecc(2,1:n_img_tipo1),'ok','MarkerSize',10)

% hold on, plot(Proyecc(1,n_img_tipo1+indices_Dubai),Proyecc(2,n_img_tipo1+indices_Dubai),'ok','MarkerSize',10)
hold on, p2=plot(Proyecc(1,n_img_tipo1+indices_ebouge_msk),Proyecc(2,n_img_tipo1+indices_ebouge_msk),'.','Color',[0.8 0 0],'MarkerSize',30)
hold on, p3=plot(Proyecc(1,n_img_tipo1+indices_ebouge_msq),Proyecc(2,n_img_tipo1+indices_ebouge_msq),'.','Color',[1 0.2 0],'MarkerSize',30)
hold on, p4=plot(Proyecc(1,n_img_tipo1+indices_felix_msk),Proyecc(2,n_img_tipo1+indices_felix_msk),'.','Color',[1 0 1],'MarkerSize',30)
hold on, p5=plot(Proyecc(1,n_img_tipo1+indices_felix_msq),Proyecc(2,n_img_tipo1+indices_felix_msq),'.','Color',[1 0.6 1],'MarkerSize',30)
hold on, p6=plot(Proyecc(1,n_img_tipo1+indices_Fiasko_msk),Proyecc(2,n_img_tipo1+indices_Fiasko_msk),'.','Color',[1 1 0],'MarkerSize',30)
hold on, p7=plot(Proyecc(1,n_img_tipo1+indices_Fiasko_msq),Proyecc(2,n_img_tipo1+indices_Fiasko_msq),'.','Color',[1 1 0.6],'MarkerSize',30)
hold on, p8=plot(Proyecc(1,n_img_tipo1+indices_indigo_mib),Proyecc(2,n_img_tipo1+indices_indigo_mib),'.','Color',[0 0 0],'MarkerSize',30)
hold on, p9=plot(Proyecc(1,n_img_tipo1+indices_itchy_mib),Proyecc(2,n_img_tipo1+indices_itchy_mib),'.','Color',[0 0 1],'MarkerSize',30)
hold on, p10=plot(Proyecc(1,n_img_tipo1+indices_itchy_mst),Proyecc(2,n_img_tipo1+indices_itchy_mst),'.','Color',[0 0.6 1],'MarkerSize',30)
hold on, p11=plot(Proyecc(1,n_img_tipo1+indices_izmir_mib),Proyecc(2,n_img_tipo1+indices_izmir_mib),'.','Color',[0.6 0.4 0],'MarkerSize',30)

legend([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11],{'Control','Ebouge MSK','Ebouge MSQ','Felix MSK','Felix MSQ','Fiasko MSK','Fiasko MSQ','Indigo MIB','Itchy MIB','Itchy MST','Izmir MIB'})
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)

if isdir(['PCA coloreado\' date])~=1
    mkdir(['PCA coloreado\' date]);
end

saveas(h,['PCA coloreado\' date '\PCA_Health_GRMD_full_colours_' date '.jpg'])

close all

%%Representamos todo sin fiasko,eboune y felix

h=figure; plot(Proyecc(1,1:n_img_tipo1),Proyecc(2,1:n_img_tipo1),'ok','MarkerSize',10)

% hold on, plot(Proyecc(1,n_img_tipo1+indices_Dubai),Proyecc(2,n_img_tipo1+indices_Dubai),'.g','MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_indigo_mib),Proyecc(2,n_img_tipo1+indices_indigo_mib),'.','Color',[0 0 0],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_itchy_mib),Proyecc(2,n_img_tipo1+indices_itchy_mib),'.','Color',[0 0 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_itchy_mst),Proyecc(2,n_img_tipo1+indices_itchy_mst),'.','Color',[0 0.6 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_izmir_mib),Proyecc(2,n_img_tipo1+indices_izmir_mib),'.','Color',[0.6 0.4 0],'MarkerSize',30)

legend('Control','Indigo MIB','Itchy MIB','Itchy MST','Izmir MIB')
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)
saveas(h,['PCA coloreado\' date '\PCA_Health_GRMD_without_colours_fiasko_felix_eboune_' date '.jpg'])

close all

%%Representamos control con fiasko,eboune y felix

h=figure; plot(Proyecc(1,1:n_img_tipo1),Proyecc(2,1:n_img_tipo1),'ok','MarkerSize',10)

hold on, plot(Proyecc(1,n_img_tipo1+indices_ebouge_msk),Proyecc(2,n_img_tipo1+indices_ebouge_msk),'.','Color',[0.8 0 0],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_ebouge_msq),Proyecc(2,n_img_tipo1+indices_ebouge_msq),'.','Color',[1 0.2 0],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_felix_msk),Proyecc(2,n_img_tipo1+indices_felix_msk),'.','Color',[1 0 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_felix_msq),Proyecc(2,n_img_tipo1+indices_felix_msq),'.','Color',[1 0.6 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_Fiasko_msk),Proyecc(2,n_img_tipo1+indices_Fiasko_msk),'.','Color',[1 1 0],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_Fiasko_msq),Proyecc(2,n_img_tipo1+indices_Fiasko_msq),'.','Color',[1 1 0.6],'MarkerSize',30)

legend('Control','Ebouge MSK','Ebouge MSQ','Felix MSK','Felix MSQ','Fiasko MSK','Fiasko MSQ')
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)
saveas(h,['PCA coloreado\' date '\PCA_Health_GRMD_without_colours_izmir_itchy_indigo_dubai_' date '.jpg'])
close all

%% Representacion del PCA Health VS GRMD MIB

load ('PCA_Health_GRMD_MIB_selection_cc_69.mat')

Proyecc=Proy{1,1};
h=figure; plot(Proyecc(1,1:n_img_tipo1),Proyecc(2,1:n_img_tipo1),'ok','MarkerSize',10)

% hold on, plot(Proyecc(1,n_img_tipo1+indices_Dubai_musc),Proyecc(2,n_img_tipo1+indices_Dubai_musc),'.g','MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_indigo_musc),Proyecc(2,n_img_tipo1+indices_indigo_musc),'.','Color',[0 0 0],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_itchy_musc),Proyecc(2,n_img_tipo1+indices_itchy_musc),'.','Color',[0 0 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+indices_izmir_musc),Proyecc(2,n_img_tipo1+indices_izmir_musc),'.','Color',[0.6 0.4 0],'MarkerSize',30)

legend('Control','Indigo MIB','Itchy MIB','Izmir MIB')
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)
saveas(h,['PCA coloreado\' date '\PCA_Health_GRMD_MIB_colours_' date '.jpg'])

close all
%% Representacion del PCA Health VS GRMD MIB similar age

load ('PCA_Health_without_Dubai_GRMD_MIB_similar_age_selection_cc_69.mat')

Proyecc=Proy{1,1};
h=figure; plot(Proyecc(1,1:n_img_tipo1-max(indices_Dubai)),Proyecc(2,1:n_img_tipo1-max(indices_Dubai)),'ok','MarkerSize',10)

hold on, plot(Proyecc(1,n_img_tipo1-max(indices_Dubai)+indices_indigo_musc_edad),Proyecc(2,n_img_tipo1-max(indices_Dubai)+indices_indigo_musc_edad),'.','Color',[0 0 0],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1-max(indices_Dubai)+indices_itchy_musc_edad),Proyecc(2,n_img_tipo1-max(indices_Dubai)+indices_itchy_musc_edad),'.','Color',[0 0 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1-max(indices_Dubai)+indices_izmir_musc_edad),Proyecc(2,n_img_tipo1-max(indices_Dubai)+indices_izmir_musc_edad),'.','Color',[0.6 0.4 0],'MarkerSize',30)

legend('Control','Indigo MIB','Itchy MIB','Izmir MIB')
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)
saveas(h,['PCA coloreado\' date '\PCA_Health_without_Dubai_GRMD_MIB_similar_age_colours_' date '.jpg'])

close all

%% Representacion del PCA Itchy VS Izmir & Indigo (MIB)

load('PCA_Itchy_against_Izmir & Indigo_selection_cc_69.mat')

n_img_tipo1=size(Itchy_cc,1);

Proyecc=Proy{1,1};
h=figure; plot(Proyecc(1,1:n_img_tipo1),Proyecc(2,1:n_img_tipo1),'.','Color',[0 0 1],'MarkerSize',30)

hold on, plot(Proyecc(1,n_img_tipo1+[1:size(Izmir_cc,1)]),Proyecc(2,n_img_tipo1+[1:size(Izmir_cc,1)]),'.','Color',[0.6 0.4 0],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+size(Izmir_cc,1)+[1:size(Indigo_cc,1)]),Proyecc(2,n_img_tipo1+size(Izmir_cc,1)+[1:size(Indigo_cc,1)]),'.','Color',[0 0 0],'MarkerSize',30)

legend('Itchy MIB','Izmir MIB','Indigo MIB')
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)
saveas(h,['PCA coloreado\' date '\PCA_Itchy_against_Izmir & Indigo_MIB_colours_' date '.jpg'])

close all

%% Representacion del PCA Izmir VS Itchy & Indigo (MIB)

load('PCA_Izmir_against_Itchy & Indigo_selection_cc_69.mat')

n_img_tipo1=size(Izmir_cc,1);

Proyecc=Proy{1,1};
h=figure; plot(Proyecc(1,1:n_img_tipo1),Proyecc(2,1:n_img_tipo1),'.','Color',[0.6 0.4 0],'MarkerSize',30)

hold on, plot(Proyecc(1,n_img_tipo1+[1:size(Itchy_cc,1)]),Proyecc(2,n_img_tipo1+[1:size(Itchy_cc,1)]),'.','Color',[0 0 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+size(Itchy_cc,1)+[1:size(Indigo_cc,1)]),Proyecc(2,n_img_tipo1+size(Itchy_cc,1)+[1:size(Indigo_cc,1)]),'.','Color',[0 0 0],'MarkerSize',30)

legend('Izmir MIB','Itchy MIB','Indigo MIB')
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)
saveas(h,['PCA coloreado\' date '\PCA_Izmir_against_Itchy & Indigo_MIB_colours_' date '.jpg'])

close all

%% Representacion del PCA Indigo VS Itchy & Izmir (MIB)

load('PCA_Indigo_against_Itchy & Izmir_selection_cc_69.mat')

n_img_tipo1=size(Indigo_cc,1);

Proyecc=Proy{1,1};
h=figure; plot(Proyecc(1,1:n_img_tipo1),Proyecc(2,1:n_img_tipo1),'.','Color',[0 0 0],'MarkerSize',30)

hold on, plot(Proyecc(1,n_img_tipo1+[1:size(Itchy_cc,1)]),Proyecc(2,n_img_tipo1+[1:size(Itchy_cc,1)]),'.','Color',[0 0 1],'MarkerSize',30)
hold on, plot(Proyecc(1,n_img_tipo1+size(Itchy_cc,1)+[1:size(Izmir_cc,1)]),Proyecc(2,n_img_tipo1+size(Itchy_cc,1)+[1:size(Izmir_cc,1)]),'.','Color',[0.6 0.4 0],'MarkerSize',30)

legend('Indigo MIB','Itchy MIB','Izmir MIB')
legend('Location','eastoutside')

stringres=strcat(num2str(indice_cc_seleccionadas));
title(stringres)
saveas(h,['PCA coloreado\' date '\PCA_Indigo_against_Itchy & Izmir_MIB_colours_' date '.jpg'])

close all

cd ..\Codigo
