%% Organizar_grupos_PCA

clear all
close all

cd ..
cd Caracteristicas_extraidas

%Leer todas las fechas de la ruta
%Obtenemos la ruta actual
[stat,struc] = fileattrib;
PathCurrent = struc.Name;

lee_fechas = dir(PathCurrent);
lee_fechas(1:2)=[];

Health_cc=[];
Health_without_Dubai_cc=[];
GRMD_cc=[];
GRMD_MIB_cc=[];
GRMD_MIB_38_39_41_cc=[];
Itchy_cc=[];

%%Indices para pca global

    indices_itchy=[];
    indices_itchy_mib=[];
    indices_itchy_mst=[];
    
    indices_izmir=[];
    indices_izmir_mib=[];
    indices_izmir_mst=[];

    indices_indigo=[];
    indices_indigo_mib=[];
    indices_indigo_mst=[];
    
    indices_ebouge=[];
    indices_ebouge_msk=[];
    indices_ebouge_msq=[];
    
    indices_felix=[];
    indices_felix_msk=[];
    indices_felix_msq=[];
    
    indices_Fiasko=[];
    indices_Fiasko_msk=[];
    indices_Fiasko_msq=[];
    
    indices_Dubai=[];
    
    
    
    
%%Indices para pca por musculo mib   
    
    indices_itchy_musc=[];
    indices_izmir_musc=[];
    indices_indigo_musc=[];
%     indices_Dubai_musc=[];

%%Indices para pca por musculo mib y edad ~39
%     
%     indices_itchy_musc_edad=[];
%     indices_izmir_musc_edad=[];
%     indices_indigo_musc_edad=[];
    
    indice=1;
    indice_musc=1;
    indice_musc_edad=1;
    indice_health=1;


    
for i=1:length(lee_fechas)
    
    fecha=lee_fechas(i).name;
    cd (fecha)
    
    [stat,struc] = fileattrib;
    PathCurrent = struc.Name;
    
    lee_fotos = dir(PathCurrent);
    lee_fotos(1:2)=[];
    
    for j=1:length(lee_fotos)
        foto=lee_fotos(j).name;
        cd (foto)

        [stat,struc] = fileattrib;
        PathCurrent = struc.Name;
        lee_rois = dir([PathCurrent '\*.mat']);
        
        for k=1:length(lee_rois)
        
       
            load (lee_rois(k).name)
            cc=[Mean_Area, Std_Area, Mean_Area_rojas,Std_Area_rojas, Mean_Area_negras, Std_Area_negras, Mean_mayor_eje, Mean_menor_eje,Mean_Relacion_ejes,Std_Relacion_ejes,Mean_Pix_region_convexa,Std_Pix_region_convexa,Mean_relacion_areas,Std_relacion_areas,Mean_vecinos,Std_vecinos,Std_vecinos_de_rojos,Std_vecinos_de_negros,Mean_vecinos_rojos_de_rojos,Mean_vecinos_negros_de_rojos,Mean_vecinos_rojos_de_negros,Mean_vecinos_negros_de_negros,Mean_Relacion_areas_vecindad,Std_Relacion_areas_vecindad, Mean_relacion_eje_mayor_vecinos,Std_relacion_eje_mayor_vecinos,Mean_relacion_eje_menor_vecinos,Std_relacion_eje_menor_vecinos,Mean_relacion_relacion_ejes_vecinos,Std_relacion_relacion_ejes_vecinos,Mean_relacion_Pix_region_convexa_vecinos,Std_relacion_Pix_region_convexa_vecinos,Mean_relacion_Angulos_vecinos,Std_relacion_Angulos_vecinos,Mean_suma_pesos,Desv_suma_pesos,Mean_pesos_celulas_negras,Desv_pesos_celulas_negras,Mean_pesos_celulas_rojas,Desv_pesos_celulas_rojas,Mean_Coef_cluster,Desv_Coef_cluster,Mean_Coef_cluster_negras,Desv_Coef_cluster_negras,Mean_Coef_cluster_rojas,Desv_Coef_cluster_rojas,Mean_excentricidad,Desv_excentricidad,Mean_excentricidad_negras,Desv_excentricidad_negras,Mean_excentricidad_rojas,Desv_excentricidad_rojas,Mean_BC,Desv_BC,Mean_BC_negras,Desv_BC_negras,Mean_BC_rojas,Desv_BC_rojas,Mean_dist,Desv_dist,Mean_dist_negras_negras,Desv_dist_negras_negras,Mean_dist_negras_rojas,Desv_dist_negras_rojas,Mean_dist_rojas_rojas,Desv_dist_rojas_rojas,Mean_dist_rojas_negras,Desv_dist_rojas_negras,Average_celulas_rojas];

    %         load 'Nucleos_cc'
    %         cc=[cc,media_picos_en_cel_water_val,media_picos_en_colageno_val,media_picos_celula_total_val,desv_n_picos_cwater,desv_n_picos_colageno,desv_n_picos_cel_total,Porcentaje_Area_Objeto_water_cell,Porcentaje_Area_Objeto_colageno,Porcentaje_Area_Objeto_cell,desv_area_obj_cwater,desv_area_obj_colageno,desv_area_obj_cel_total];


            switch foto(1:5)

                %'09065' se trata de Dubai
                case {'12007','12008','12011','09065'}   %Numeraciones de perros sanos 

                    Health_cc=[Health_cc; cc];   %Guardamos las cc de todos los sanos (todos tienen 39 semanas y las muestras son tomadas del musculo MIB)
                    
                    if strcmp(foto(1:5),'09065')==1
                        indices_Dubai=[indices_Dubai, indice_health];
                    else
                        Health_without_Dubai_cc=[Health_without_Dubai_cc; cc];
                    end
                    
                    indice_health=indice_health+1;

                case {'14074','14095','10107','11031','11043'}  % Numeraciones de los perros GRMD
                    
                    GRMD_cc=[GRMD_cc; cc];                              %Tomamos todos los perros con mutación
                    
                                        
                    switch foto(1:5)
                        
                        %%%Perro Itchy
                        case {'14074'}
                            
                            indices_itchy=[indices_itchy, indice];
                                                        
                            switch foto(7:9)
                                
                                case {'MIB'}
                                    indices_itchy_mib=[indices_itchy_mib,indice];
                                    
                                case {'MST'}    
                                    indices_itchy_mst=[indices_itchy_mst,indice];
                            end
                            
                            
                        %%%Perro Izmir    
                        case {'14095'}
                            if foto(6)=='6'
                                indices_izmir=[indices_izmir,indice];
                                
                                switch foto(7:9)
                                    case {'MIB'}
                                        indices_izmir_mib=[indices_izmir_mib,indice];

                                    case {'MST'}    
                                        indices_izmir_mst=[indices_izmir_mst,indice];
                                end
                                
                        %%%Perro Indigo
                            else
                                indices_indigo=[indices_indigo,indice];
                                
                                switch foto(7:9)
                                    case {'MIB'}
                                        indices_indigo_mib=[indices_indigo_mib,indice];

                                    case {'MST'}    
                                        indices_indigo_mst=[indices_indigo_mst,indice];
                                end                                                          
                               
                            end
                            
                            
                        %%%Perro ebouge    
                        case {'10107'}
                            indices_ebouge=[indices_ebouge,indice];
                            
                            switch foto(7:9)
                                case {'MSK'}
                                    indices_ebouge_msk=[indices_ebouge_msk,indice];

                                case {'MSQ'}    
                                    indices_ebouge_msq=[indices_ebouge_msq,indice];
                            end    
                            
                        %%%%Perro Felix    
                        case {'11031'}
                            indices_felix=[indices_felix,indice];
                            
                            switch foto(8:9)
                                case {'SK'}
                                    indices_felix_msk=[indices_felix_msk,indice];

                                case {'SQ','Qg'}    
                                    indices_felix_msq=[indices_felix_msq,indice];
                            end
                            
                            
                        %%%Perro fiasko    
                        case {'11043'}
                            indices_Fiasko=[indices_Fiasko,indice];
                            
                            switch foto(8:9)
                                case {'SK'}
                                    indices_Fiasko_msk=[indices_Fiasko_msk,indice];

                                case {'SQ','Qg'}    
                                    indices_Fiasko_msq=[indices_Fiasko_msq,indice];
                            end
                            
                        %%%Perro dubai    
%                         case {'09065'}
%                             indices_Dubai=[indices_Dubai,indice];
                                                      
                            
                    end
                    indice=indice+1; 
                        
                    %%Aqui guardamos las matrices que se especializan
                    %%en músculo mib (itchy,izmir,indigo y dubai) 
                    
                    if strcmp(foto(7:8),'MI')==1 || strcmp(foto(7:8),'IB')==1

                        GRMD_MIB_cc=[GRMD_MIB_cc;cc];
                        
                        switch foto(1:5)
                            case {'14074'}
                                indices_itchy_musc=[indices_itchy_musc, indice_musc];
                            case {'14095'}
                                if foto(6)=='6'
                                    indices_izmir_musc=[indices_izmir_musc,indice_musc];
                                else
                                    indices_indigo_musc=[indices_indigo_musc,indice_musc];
                                end
%                             case {'09065'}
%                                 indices_Dubai_musc=[indices_Dubai_musc,indice_musc];
                        end
                        indice_musc=indice_musc+1; 

%                         %%Aqui guardamos las matrices que se especializan en edad y musculo mib(itchy,izmir e indigo)
%                         
%                         if strcmp(foto(1:2),'14')==1
% 
%                             GRMD_MIB_38_39_41_cc=[GRMD_MIB_38_39_41_cc;cc];
%                             
%                             switch foto(1:5)
%                                 case {'14074'}
%                                     indices_itchy_musc_edad=[indices_itchy_musc_edad, indice_musc_edad];
%                                 case {'14095'}
%                                     if foto(6)=='6'
%                                         indices_izmir_musc_edad=[indices_izmir_musc_edad,indice_musc_edad];
%                                     else
%                                         indices_indigo_musc_edad=[indices_indigo_musc_edad,indice_musc_edad];
%                                     end                               
%                             end
%                             indice_musc_edad=indice_musc_edad+1; 
% 
%                         end

                    end

            end
            
        end
                                               
           cd ..       
 
    end
    cd ..
end

Itchy_cc=GRMD_MIB_cc(indices_itchy_musc,:);
Izmir_cc=GRMD_MIB_cc(indices_izmir_musc,:);
Indigo_cc=GRMD_MIB_cc(indices_indigo_musc,:);


cd ..\Matrices_cc

save ('Matrices_cc.mat','Health_cc','Health_without_Dubai_cc','GRMD_cc','GRMD_MIB_cc','Itchy_cc','Izmir_cc','Indigo_cc')
save('Indices_total.mat','indices_itchy','indices_izmir','indices_indigo','indices_ebouge','indices_felix','indices_Fiasko','indices_Dubai','indices_itchy_mib','indices_itchy_mst','indices_izmir_mib','indices_izmir_mst','indices_indigo_mib','indices_indigo_mst','indices_ebouge_msk','indices_ebouge_msq','indices_felix_msk','indices_felix_msq','indices_Fiasko_msk','indices_Fiasko_msq')
save('Indices_musc.mat','indices_itchy_musc','indices_izmir_musc','indices_indigo_musc')%,'indices_Dubai_musc')
% save('Indices_musc_edad.mat','indices_itchy_musc_edad','indices_izmir_musc_edad','indices_indigo_musc_edad')
cd ..\Codigo

