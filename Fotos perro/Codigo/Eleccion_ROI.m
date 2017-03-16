%Eleccion_ROI
function [f,centroides_ROI, mask_ROIS, area_rois]=Eleccion_ROI(Img,f)


     
%%Elegimos el numero de ROIS que queremos, con unas medidas determinadas. Colocamos las ROIS en la imagen con un api
%%proporcionado por matlab.

%%Habrá varias interactuaciones con pantalla: -elegir el nº de ROIs,
%%-colocar las ROIs y -verificar su corrección

bandera=0;
    while bandera~=1
        
        centroides_ROI=[];
        mask_ROIS=[];
        area_rois=[];
        
        NUMEROROI = input('How many ROIs would you like to get? : ');

        for i=1:NUMEROROI
            %limite=100;
            limite=0;
            diametro=1400;


            %%Definimos el circulo que va a delimitar el ROI
            h = imellipse(gca, [limite limite diametro diametro]); 
            api = iptgetapi(h);

            fcn = getPositionConstraintFcn(h);

            api.setPositionConstraintFcn(fcn);
            pause

            BW = createMask(h);       

            area_ROI=regionprops(BW,'Area');
            area_ROI=cat(1,area_ROI.Area);
            area_rois(i)=area_ROI;

            mask_ROIS{i,1}=BW;

            centroide=regionprops(BW,'Centroid');
            centroide_ROI=cat(1,centroide.Centroid);
            centroides_ROI{i,1}=centroide_ROI;

        end 
        
        bandera = input('If ROIs are correct click ---> 1 , else click any different number ');
        pause;
        
        if bandera ~=1
            close all
            f=figure;imshow(Img)
                   
        end
    end
end