%%Llamada extracción de característcas no dapi
clear all
close all

cd ..\..
cd Segmentadas

%Leer todas las fechas de la ruta
%Obtenemos la ruta actual
[stat,struc] = fileattrib;
PathCurrent = struc.Name;


lee_fechas = dir(PathCurrent);
lee_fechas(1:2)=[];



for i=1:length(lee_fechas)
    
    fecha=lee_fechas(i).name;
    cd (fecha)
    
    [stat,struc] = fileattrib;
    PathCurrent = struc.Name;
    
    lee_fotos = dir(PathCurrent);
    lee_fotos(1:2)=[];
    lee_fotos(end-2:end)=[];
    
    cd ..\..
    cd Codigo
    for j=1:length(lee_fotos)
        cd ('..\Seleccion_Roi\7 Seleccion_Roi_Pedro\Imagenes')
        cd (fecha(end-7:end))
        cd (lee_fotos(j).name)
        
        load ('Datos_Rois.mat')
        
        cd ..\..\..\..\..\Codigo
        for n_roi=1:size(mask_ROIS,1)
            foto=lee_fotos(j).name;
            Extraccion_cc_no_DAPI_roi(fecha,foto,n_roi)
            close all
        end
    end
    cd ..
    cd Segmentadas
    
 
end

cd ..
cd ('Codigo\Llamadas')

