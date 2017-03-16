% Llamada segmentación de núcleos mediante el uso de ROIS

close all
clear all

cd ('..\..\Segmentadas')

%Leer todas las fechas de la ruta
%Obtenemos la ruta actual
[stat,struc] = fileattrib;
PathCurrent = struc.Name;

%Obtenemos la lista de carpetas fecha

lee_fechas = dir(PathCurrent);
lee_fechas(1:2)=[];

for i=7:length(lee_fechas)
    
    fecha=lee_fechas(i).name;
    cd (fecha)
    
    %Obtenemos la lista de fotos de cada carpeta
    [stat,struc] = fileattrib;
    PathCurrent = struc.Name;
    
    lee_fotos = dir(PathCurrent);
    lee_fotos(1:2)=[];
    lee_fotos(end-2:end)=[];
    
   
    for j=1:length(lee_fotos)
       
        foto=lee_fotos(j).name;
        cd (foto)
        Img=imread([foto ' P.jpg']);
                  
        cd ('..\..\..\Seleccion_Roi\7 Seleccion_Roi_Pedro\Imagenes')
        cd (fecha(end-7:end))
        cd (foto)

        %Cargamos los datos de las rois
        
        load ('Datos_Rois.mat')

        cd ..\..\..\..\..\Deteccion_dapi
        
        %Creamos el directorio de la carpeta fecha
        
        if isdir(['6 Deteccion Dapi ' fecha(end-7:end)]) == 0
            mkdir(['6 Deteccion Dapi ' fecha(end-7:end)]);
        end
        
        cd (['6 Deteccion Dapi ' fecha(end-7:end)])
              
        %%Creamos el directorio de la foto
        
        if isdir(foto) == 0
            mkdir(foto);
        end
        
        %%Guardamos la imagen dentro de la carpeta destino de detección
        %%dapi
        imwrite(Img,[foto '\' foto ' P.jpg']);
        
        %%Hallamos el número de ROIS a partir de las caracteristicas
        %%extraidas
        
        cd (['..\..\Caracteristicas_extraidas\8 Caracteristicas_extraidas ' fecha(end-7:end)])
        cd (foto)
        num_rois=ls('*.mat');
        num_rois=size(num_rois,1);
        
        %%Volvemos a codigo para llamar la funcion de extraccion de nucleos
        cd ..\..\..\Codigo
        
        
        %%Llamamos a la funcion de extracción de nucleos tantas veces
        %%como rois existan.
        parfor n_roi=1:num_rois
            Extraccion_nucleos_roi(fecha,foto,n_roi)
            close all
        end
            
        cd ..\Segmentadas
        cd (fecha)
        
    end
    cd ..
end