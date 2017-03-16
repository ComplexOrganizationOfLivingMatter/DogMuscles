
function Extraccion_nucleos(fecha,nombre)
%%Cargamos las imágenes

    carpeta_fecha=fecha;


    cd ..
    cd 'Segmentadas'
    
    cd (fecha)
    cd (nombre)
    


    nombre_r=strcat(nombre,'r');
    nombre_g=strcat(nombre,'g');
    nombre_b=strcat(nombre,'b');



    Img_r=imread(nombre_r,'jpeg');
    Img_g=imread(strcat(nombre_g,'_edited'),'jpeg');
    Img_g_original=imread(nombre_g,'jpeg');

    Img_b=imread(nombre_b,'jpeg');


    Img=(Img_r+Img_g+Img_b);
    
    R=Img(:,:,1);
    G=Img(:,:,2);
    B=im2double(Img(:,:,3));    
    [H,W,c]=size(Img);
    
    cd ..\..\..
    cd 'Imagenes_procesadas'
    cd (carpeta_fecha)
    cd (nombre)
    cd 'Imagenes_necesarias'
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%Extracion plano azul %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% datos_celulas=cell(1,6);
% datos_celulas_water=cell(1,6);
% datos_celulas_menos_water=cell(1,6);
% datos_objetos=cell(1,6);
% datos_objetos_pico=cell(1,6);
% datos_objetos_pico_centroide=cell(1,6);
% ruta=strcat(nombre3,'-images\imagenes_necesarias');
% load(ruta,'contorno','contorno_water','mascara_mejorada')
load('imagenes_necesarias.mat','contorno','contorno_water','mascara_mejorada','celulas')    




%% Obtenemos objetos
    
%Suavizamos contornos
h1=fspecial('average') ;
B=imfilter(B,h1);

%Binarizamos imagen
% level=0.05; % Nivel de intensidad al menos de 5%
BW = im2bw(B);
%figure, imshow(BW)

%Eliminamos objetos menores de 300 pix
BWmin= bwareaopen(BW,400);


%Redondeamos objetos
% se = strel('disk',2);
% BWmin = imdilate(BWmin,se);
% figure, imshow(BWmin)
%imwrite(BW,'Blue segmentation.bmp')

% etiqueto objetos
L = bwlabel(BWmin,8);
% figure, imshow(L)

% Clasifico objetos
BWmax=zeros(H,W);

Area_ob = regionprops(L, 'area');
Area_ob = cat(1, Area_ob.Area);
area_media=mean(Area_ob);
for i=1:max(max(L))

    BWmax(L==i)=i;
    a=BWmax(L==i);
    rest=max(B(L==i))-min(B(L==i));
    
    if Area_ob(i)>(area_media*1.2)
            int_umbral=min(B(L==i))+rest*0.85;       
    else
        int_umbral=min(B(L==i))+rest*0.40;
    end
    
    b=find(B(L==i)<int_umbral);
    a(b)=0;
    BWmax(L==i)=a;
    BWmax(BWmax>0)=1;
    
end
    
    
% level=0.65*max(max(B)); % Nivel de intensidad al menos de 75%
% 
% 
% 
% BWmax = im2bw(B, level);
BWmax= bwareaopen(BWmax,80);
% figure, imshow(BWmax)



Lmax = bwlabel(BWmax,8);
% figure, imshow(Lmax)

% celulas=bwlabel(1-contorno,8);
% % figure, imshow(celulas)





%% La idea aqui es concatenar


aux=zeros(H,W);
s  = regionprops(Lmax, 'PixelList');

for i=1:max(max(Lmax))
    data=[];
    data=s(i,1).PixelList;
    aux(L==L(data(1,2),data(1,1)))=1;
end
% figure, imshow(aux)
BWmin=aux;
L = bwlabel(BWmin,8);
% 
% %% Seleccion de picos(nodos)
% masc_nucleo_TOTAL=zeros(H,W);
% 
% for i=1:max(max(L))
%     %% Creamos imagen mascara del objeto i
%     masc_obj=zeros(H,W);
%     masc_obj(L==i)=1;
%     
%     %% Obtenemos informacion del plano azul del objeto i
%     obj=masc_obj.*B;
%     
%     %% Umbralizamos para dividir el objeto i en unidades mas pequeñas
%     level=graythresh(obj(masc_obj==1));
%     masc_sub_obj_TOTAL = im2bw(obj, level);
%     L_masc_sub_obj_TOTAL = bwlabel(masc_sub_obj_TOTAL,8);
%     % Obtenemos area de cada subobjeto j del objeto i
%     AREA = regionprops(L_masc_sub_obj_TOTAL, 'area');
%     AREA = cat(1, AREA.Area);
%     
%     %% Recorremos todos los subobjetos j del objeto i
%     for j=1:max(max(L_masc_sub_obj_TOTAL))
%         
%         % Creamos imagen mascara del subobjeto j
%         masc_sub_obj=zeros(H,W);
%         masc_sub_obj(L_masc_sub_obj_TOTAL==j)=1;
%         sub_obj=masc_sub_obj.*B;
%         % Si subobjeto no presenta una intensidad maxima mayor de 0.8 no es
%         % objeto pico y lo eliminamos
%         if max(sub_obj(L_masc_sub_obj_TOTAL==j)) <=0.80
%             masc_sub_obj(L_masc_sub_obj_TOTAL==j)=0;
%         else % En caso contrario volvemos a umbralizar subobjeto para obtener objetos picos definitivos
%             level=graythresh(sub_obj(L_masc_sub_obj_TOTAL==j));
%             masc_nucleo = im2bw(sub_obj, level);
%             L_masc_nucleo = bwlabel(masc_nucleo,8);
%             
%             % Obtenemos area de los picos perteneciente a subobjeto j del
%             % objeto i
%             AREAP = regionprops(L_masc_nucleo, 'area');
%             AREAP = cat(1, AREAP.Area);
%             %Recorremos todos los picos k en potencia eliminos todos
%             %aquellos que tenga area inferior a 1% del area del subobjeto j
%             %al que pertenece , ademas comprobamos que picos sean mayor que
%             %0.99
%             for k=1:max(max(L_masc_nucleo))
%                 sub_nucleo=masc_nucleo.*B;
%                 proporcion=AREAP(k)/AREA(j);
%                 if proporcion <= 0.02
%                     masc_nucleo(L_masc_nucleo==k)=0;
%                 end
%                 
%                 if max(sub_nucleo(L_masc_nucleo==k)) <=0.9
%                     masc_nucleo(L_masc_nucleo==k)=0;
%                 end
%             end
%             masc_nucleo_TOTAL=masc_nucleo_TOTAL+masc_nucleo;
%         end
%     end
% end

masc_nucleo_TOTAL=Lmax;

% figure, imshow(masc_nucleo_TOTAL)






%%Etiquetamos los núcleos 
L_masc_nucleo_TOTAL = bwlabel(masc_nucleo_TOTAL,8);
s  = regionprops(L_masc_nucleo_TOTAL, 'centroid');
centroids = cat(1, s.Centroid);
centroids=fliplr(round(centroids));
aux2=zeros(H,W);
for i=1:size(centroids,1)
    aux2(centroids(i,1),centroids(i,2))=i;
end




%% Relaciono datos con respecto a celulas
L_contorno_water = bwlabel(contorno_water,8);
area_water = regionprops(L_contorno_water, 'Area');
area_water = cat(1, area_water.Area);
[value,ix]=max(area_water);
contorno_water(L_contorno_water==ix)=0;
contorno_water(contorno_water~=0)=1;
L_contorno_water = bwlabel(contorno_water,8);
area_water = regionprops(L_contorno_water, 'Area');
area_water = cat(1, area_water.Area);
mascara_mejorada=im2bw(mascara_mejorada,0.5);
%%Recorro celulas para terminar limpiar conexiones restantes
for i=1:max(max(L_contorno_water))
    v=find(mascara_mejorada(L_contorno_water==i)==0);
    probabilidad_eliminar=(length(v)/area_water(i))*100;
    if probabilidad_eliminar>80
        contorno_water(L_contorno_water==i)=0;
    end
end
L_contorno_water = bwlabel(contorno_water,8);
L_contorno_water(L_contorno_water>0)=1;
L_contorno_water=L_contorno_water.*celulas;

contorno_total_menos_water=1-contorno-contorno_water;
L_contorno_total_menos_water = bwlabel(contorno_total_menos_water,8);

s  = regionprops(celulas, 'Area'); 
Area_celula = cat(1, s.Area);
s  = regionprops(L_contorno_water, 'Area');
Area_water = cat(1, s.Area);
s  = regionprops(L_contorno_total_menos_water, 'Area');
Area_total_menos_water = cat(1, s.Area);
s  = regionprops(L, 'Area');
Area_objeto = cat(1, s.Area);
s  = regionprops(L_masc_nucleo_TOTAL, 'Area');
Area_Nodo = cat(1, s.Area);



















for i=1:max(max(celulas))
    
    aux=zeros(H,W);
    aux(celulas==i)=1;
    %%Celula
    datos_celulas{i,1}=i;
    
    %%Celula Water
    objeto_azul=aux.*L_contorno_water; %%concatenacion de celulas borde 1, con celulas originales
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas{i,2}=objetos; 
    
    %%Espacio sobrante
    objeto_azul=aux.*L_contorno_total_menos_water; %%concatenacion de celulas borde 1, con el espacio intercelular
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas{i,3}=objetos; %Espacio sobrante que pertenece a cada célula water
    
    %%Objeto
    objeto_azul=aux.*L; %%concatenacion de núcleos con células (Célula X contiene L núcleos)
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas{i,4}=objetos; %nucleos que pertenecen a cada célula
    
    %%Nodo
    objeto_azul=aux.*L_masc_nucleo_TOTAL; %concatenacion de picos de nucleo con celulas borde 1.
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas{i,5}=objetos;
    
    %%Centroides
    objetos=aux2.*aux; %% concatenacion de centroides de picos de nucleo con las celulas borde 1.
    objetos=unique(objetos);
    objetos=objetos(2:end,1);
    datos_celulas{i,6}=centroids(objetos,:);
    
    %% Area Celula
    datos_celulas{i,7}=Area_celula(i);
    
    %% Area pico que solapan con la ´celula i
    s  = regionprops(objeto_azul, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_celulas{i,8}=Area_nodo_celula(datos_celulas{i,5});
    
    %% Porcentaje de solapamiento de nodo respecto al area total del nodo
    s  = regionprops(objeto_azul, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_celulas{i,9}=(datos_celulas{i,8}./Area_Nodo(datos_celulas{i,5}))*100;
    
    %% Area objeto que solapan con la ´celula i
    objeto_azul=aux.*L;
    s  = regionprops(objeto_azul, 'Area');
    Area_objeto_celula = cat(1, s.Area);
    datos_celulas{i,10}=Area_objeto_celula(datos_celulas{i,4});
    
    %% Porcentaje de solapamiento de objeto respecto al area total del objeto
    datos_celulas{i,11}=(datos_celulas{i,10}./Area_objeto(datos_celulas{i,4}))*100;
end

for i=1:max(max(L_contorno_water))
    aux=zeros(H,W);
    aux(L_contorno_water==i)=1;
    
    %%Celula Water
    datos_celulas_water{i,1}=i;
    
    %%Celula 
    objeto_celulas=aux.*celulas;
    objetos=unique(objeto_celulas);
    objetos=objetos(2:end,1);
    datos_celulas_water{i,2}=objetos;
    
    %%Objeto
    objeto_azul=aux.*L; %%concatenación de nucleos con celula original (No colágeno)
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas_water{i,4}=objetos;
    
    %%Nodo
    objeto_azul=aux.*L_masc_nucleo_TOTAL; %%concatenacion de celulas originales con picos de nucleo
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas_water{i,5}=objetos;
    
    %%Centroides
    objetos=aux2.*aux; 
    objetos=unique(objetos);
    objetos=objetos(2:end,1);
    datos_celulas_water{i,6}=centroids(objetos,:);
    
    %% Area Water
    datos_celulas_water{i,7}=Area_water(i);
    
     %% Area pico que solapan con la celula water i
    s  = regionprops(objeto_azul, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_celulas_water{i,8}=Area_nodo_celula(datos_celulas_water{i,5});
    
    %% Porcentaje de solapamiento de nodo respecto al area total del nodo
    s  = regionprops(objeto_azul, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_celulas_water{i,9}=(datos_celulas_water{i,8}./Area_Nodo(datos_celulas_water{i,5}))*100;
    
    %% Area objeto que solapan con la ´celula water i
    objeto_azul=aux.*L;
    s  = regionprops(objeto_azul, 'Area');
    Area_objeto_celula = cat(1, s.Area);
    datos_celulas_water{i,10}=Area_objeto_celula(datos_celulas_water{i,4});
    
    %% Porcentaje de solapamiento de objeto respecto al area water del objeto
    datos_celulas_water{i,11}=(datos_celulas_water{i,10}./Area_objeto(datos_celulas_water{i,4}))*100;
    
end

for i=1:max(max(L_contorno_total_menos_water))
    aux=zeros(H,W);
    aux(L_contorno_total_menos_water==i)=1;
    
    %%Espacio sobrante
    datos_celulas_menos_water{i,1}=i;
    
    %%Celula 
    objeto_celulas=aux.*celulas; %concatenacion de espacio sobrante con celulas borde 1.
    objetos=unique(objeto_celulas);
    objetos=objetos(2:end,1);
    datos_celulas_menos_water{i,2}=objetos;
    
    %%Objeto
    objeto_azul=aux.*L; %concatenacion de espacio sobrante con los nucleos
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas_menos_water{i,4}=objetos;
    
    %%Nodo
    objeto_azul=aux.*L_masc_nucleo_TOTAL; %concatenacion de nucleos pico con espacio sobrante
    objetos=unique(objeto_azul);
    objetos=objetos(2:end,1);
    datos_celulas_menos_water{i,5}=objetos;
    
    %%Centroides
    objetos=aux2.*aux; %concatenacion de centroides de celulas pico con espacio sobrante
    objetos=unique(objetos);
    objetos=objetos(2:end,1);
    datos_celulas_menos_water{i,6}=centroids(objetos,:);
     
     %% Area total_menos_water
    datos_celulas_menos_water{i,7}=Area_total_menos_water(i);
    
     %% Area pico que solapan con la celula water i
    s  = regionprops(objeto_azul, 'Area');
    Area_nodo_celula = cat(1, s.Area);
   datos_celulas_menos_water{i,8}=Area_nodo_celula(datos_celulas_menos_water{i,5});
    
    %% Porcentaje de solapamiento de nodo respecto al area total del nodo
    s  = regionprops(objeto_azul, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_celulas_menos_water{i,9}=(datos_celulas_menos_water{i,8}./Area_Nodo(datos_celulas_menos_water{i,5}))*100;
    
    %% Area objeto que solapan con la ´celula total menos water i
    objeto_azul=aux.*L;
    s  = regionprops(objeto_azul, 'Area');
    Area_objeto_celula = cat(1, s.Area);
    datos_celulas_menos_water{i,10}=Area_objeto_celula(datos_celulas_menos_water{i,4});
    
    %% Porcentaje de solapamiento de objeto respecto al area water del objeto
    datos_celulas_menos_water{i,11}=(datos_celulas_menos_water{i,10}./Area_objeto(datos_celulas_menos_water{i,4}))*100;
end

% Relaciono datos con respecto a objetos
for i=1:max(max(L))
    aux=zeros(H,W);
    aux(L==i)=1;
    
    %%Objeto
    datos_objetos{i,1}=i;
    
    %%Celula 
    objeto_celulas=aux.*celulas;
    objetos=unique(objeto_celulas);
    objetos=objetos(2:end,1);
    datos_objetos{i,2}=objetos;
    
    %%Celula Water
    objeto_celulas=aux.*L_contorno_water;
    objetos=unique(objeto_celulas);
    objetos=objetos(2:end,1);
    datos_objetos{i,3}=objetos;
    
    %%Espacio sobrante
    objeto_celulas=aux.*L_contorno_total_menos_water;
    objetos=unique(objeto_celulas);
    objetos=objetos(2:end,1);
    datos_objetos{i,4}=objetos;
    
    %%Nodo
    objeto_celulas=aux.*L_masc_nucleo_TOTAL;
    objetos=unique(objeto_celulas);
    objetos=objetos(2:end,1);
    datos_objetos{i,5}=objetos;
    
    %%Centroides
     objetos=aux2.*aux;
    objetos=unique(objetos);
    objetos=objetos(2:end,1);
    datos_objetos{i,6}=centroids(objetos,:);
    
    %% Area Objeto
    s  = regionprops(L, 'Area');
    Area_objeto = cat(1, s.Area);
    datos_objetos{i,7}=Area_objeto(i);
end

% Relaciono datos con respecto a picos total
for i=1:max(max(L_masc_nucleo_TOTAL))
    aux=zeros(H,W);
    aux(L_masc_nucleo_TOTAL==i)=1;
    
    %%Nodo
    datos_objetos_pico{i,1}=i;
    
    %%Celula
    picoc=aux.*celulas;  %%diferentes celulas (borde 1px) con las que solapan los nucles
    objetos=unique(picoc);
    objetos=objetos(2:end,1);
    datos_objetos_pico{i,2}=objetos;
    
    %%Celula Water
    picow=aux.*L_contorno_water; %%celulas originales con las que solapan los nucleos
    objetos=unique(picow);
    objetos=objetos(2:end,1);
    datos_objetos_pico{i,3}=objetos;
    
    %%Espacio sobrante
    picos=aux.*L_contorno_total_menos_water; %%colageno perteneciente a la celula X con el que solapan los nucleos
    objetos=unique(picos);
    objetos=objetos(2:end,1);
    datos_objetos_pico{i,4}=objetos;
    
    %%Objeto
    pico=aux.*L; %%objeto al que pertenecen los nucleos pico
    objetos=unique(pico);
    objetos=objetos(2:end,1);
    datos_objetos_pico{i,5}=objetos;
    
    %%Centroides
    objetos=aux2.*aux;
    objetos=unique(objetos);
    objetos=objetos(2:end,1);
    datos_objetos_pico{i,6}=centroids(objetos,:);
    
    %% Area Nodo
    s  = regionprops(L_masc_nucleo_TOTAL, 'Area'); %area nucleos pico
    Area_Nodo = cat(1, s.Area);
    datos_objetos_pico{i,7}=Area_Nodo(i);
    
    %% Area pico que solapan con la celula  i
    s  = regionprops(picoc, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_objetos_pico{i,8}=Area_nodo_celula(datos_objetos_pico{i,2});
    
    %% Area pico que solapan con la celula water  i
    s  = regionprops(picow, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_objetos_pico{i,9}=Area_nodo_celula(datos_objetos_pico{i,3});
    
    %% Area pico que solapan con espacio celula menos water  i
    s  = regionprops(picos, 'Area');
    Area_nodo_celula = cat(1, s.Area);
    datos_objetos_pico{i,10}=Area_nodo_celula(datos_objetos_pico{i,4});
    
    %% Porcentaje de solapamiento de nodo solapante con celula respecto al area total del nodo
    datos_objetos_pico{i,11}=(datos_objetos_pico{i,8}./Area_Nodo(i))*100;
    
    %% Porcentaje de solapamiento de nodo solapante con celula water respecto al area total del nodo
    datos_objetos_pico{i,12}=(datos_objetos_pico{i,9}./Area_Nodo(i))*100;
    
    %% Porcentaje de solapamiento de nodo solapante con celula water- celula total respecto al area total del nodo
    datos_objetos_pico{i,13}=(datos_objetos_pico{i,10}./Area_Nodo(i))*100;
end

% Relaciono datos con respecto a picos centroide
for i=1:size(centroids,1)
    
    
    %%centroide
    datos_objetos_pico_centroide{i,1}=[centroids(i,1),centroids(i,2)];
    
    %%Celula 
    objetos=celulas(centroids(i,1),centroids(i,2));
    datos_objetos_pico_centroide{i,2}=objetos;
    
    %%Celula Water    
    objetos=L_contorno_water(centroids(i,1),centroids(i,2));
    datos_objetos_pico_centroide{i,3}=objetos;
    
    %%Espacio sobrante
    objetos=L_contorno_total_menos_water(centroids(i,1),centroids(i,2));
    datos_objetos_pico_centroide{i,4}=objetos;
    
    %%Objeto    
    objetos=L(centroids(i,1),centroids(i,2));
    datos_objetos_pico_centroide{i,5}=objetos;
    
    %%Nodo
    objetos=L_masc_nucleo_TOTAL(centroids(i,1),centroids(i,2));
    datos_objetos_pico_centroide{i,6}=objetos;
end



imwrite(B,'Blue channel.bmp')
imwrite(L,'Blue object segmentation.bmp')
imwrite(L_masc_nucleo_TOTAL,'Blue node segmentation.bmp')
save('Datos_plano_azul','L','L_masc_nucleo_TOTAL','L_contorno_water','L_contorno_total_menos_water','centroids','datos_celulas','datos_celulas_water','datos_celulas_menos_water','datos_objetos','datos_objetos_pico','datos_objetos_pico_centroide')

cd ..\..\..\..
cd 'código'

% Extraccion_nucleos
end