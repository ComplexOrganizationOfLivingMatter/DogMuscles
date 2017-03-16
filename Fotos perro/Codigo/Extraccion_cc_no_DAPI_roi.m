%% EXTRACCIÓN DE CARACTERÍSTICAS DE IMAGENES SIN DAPI (ROIS)

function Extraccion_cc_no_DAPI_roi(fecha,nombre,n_roi)

%Cargamos las imagenes necesarias y leemos la foto
cd ..
cd Segmentadas
cd (fecha)
cd (nombre)

Img=imread([nombre ' P.jpg']);

cd ([nombre ' P-images'])
load ('imagenes_necesarias.mat')

%Cargamos las ROIs
cd ..\..\..\..

cd 'Seleccion_Roi\7 Seleccion_Roi_Pedro\Imagenes'
cd (fecha(end-7:end))
cd (nombre)
load Datos_Rois
Img2=imread('Rois.bmp');

%Definimos la imagen célula, y le asignamos la region de la roi.
celulas=1-contorno;

roi=mask_ROIS{n_roi};
celulas_roi=bwlabel(celulas.*roi);




%% OBTENEMOS VECINOS

[H,W]=size(celulas_roi);
aux=celulas_roi;   %%En 'celulas' tenemos las células organizadas como en un voronoi, es decir, con bordes de 1 o 2 píxeles. De aquí obtenemos el nº de vecinos
radio=4;
num_celulas=sort(unique(aux));
num_celulas=num_celulas(2:end);
n_vecinos=zeros(length(num_celulas),1);
vecinos_real=cell(1);
for i=1 : length(num_celulas)
    BW2 = bwperim(aux==i);
    [pi,pj]=find(BW2==1);
    
    se = strel('disk',radio);
    BW2_dilate=imdilate(BW2,se);
    pixeles_vecinos=find(BW2_dilate==1);
    vecinos=unique(aux(pixeles_vecinos));
    vecinos_real{i}=vecinos(vecinos ~= 0 & vecinos ~= i);
    n_vecinos(i)=length(vecinos_real{i});
end

%% Deteccion de células validas (aquellas no pertenecientes al borde ni sus vecinos)

%Obtenemos los bordes de la roi para calcular las celulas no validas del
%borde
se=strel('disk',4);
roi_erode= imerode(roi,se);
per_roi=roi-roi_erode;



% Señalamos como celulas no validas las que estan en el borde
aux3=per_roi.*aux;
celulas_no_validas=sort(unique(aux3));
celulas_no_validas=celulas_no_validas(2:end);
celulas_validas=(1:length(num_celulas));
for i=1:length(celulas_no_validas)
    celulas_validas(celulas_validas==celulas_no_validas(i))=[];
end
celulas_validas=celulas_validas';
celulas_validas_no_borde=celulas_validas;

% Señalamos como celulas no validas las vecinas de las células del borde

celulas_no_validas_2=[];
for i=1:size(celulas_no_validas,1)
    celulas_no_validas_2=[celulas_no_validas_2;vecinos_real{celulas_no_validas(i)}];
end
celulas_no_validas_2=sort(unique(celulas_no_validas_2));

for i=1:length(celulas_no_validas_2)
    celulas_validas(celulas_validas==celulas_no_validas_2(i))=[];
end
celulas_validas=celulas_validas';
celulas_validas_no_borde_no_vecinos=celulas_validas;


% Eliminamos celulas aisladas del grupo de celulas validas sin borde
bandera=0;
celulas_validas_previa=celulas_validas_no_borde;
celulas_no_validas_previa=celulas_no_validas;
for j=1:length(celulas_validas_previa) % Bucle que recorre todas celulas validas para comprobar
    vec_cel_ind_j=vecinos_real{celulas_validas_previa(j)};
    no_coincidencia=[];
    for i=1:length(vec_cel_ind_j)% Bucle que recorre todos los vecinos de la celula bajo estudio
        no_coincidencia(i)=isempty(find(vec_cel_ind_j(i)==celulas_no_validas_previa, 1)); %Variable que vale 1 si no existe coincidencia de la celula j con alguna de las celulas no validas
    end
    if sum(no_coincidencia)==0
        celulas_validas_no_borde=celulas_validas_previa(celulas_validas_previa~=celulas_validas_previa(j));
        celulas_no_validas_no_borde=sort([celulas_validas_previa(j);celulas_no_validas_previa]);
        bandera=1;
    end
end
celulas_validas_no_borde=celulas_validas_no_borde';
if (bandera==0)
    celulas_validas_no_borde=celulas_validas_previa;
    celulas_no_validas_no_borde=celulas_no_validas_previa;
end


% Eliminamos celulas aisladas del grupo de celulas validas sin borde ni
% sus vecinos
bandera=0;
celulas_validas_previa=celulas_validas_no_borde_no_vecinos;
celulas_no_validas_previa=celulas_no_validas_2;
for j=1:length(celulas_validas_previa) % Bucle que recorre todas celulas validas para comprobar
    vec_cel_ind_j=vecinos_real{celulas_validas_previa(j)};
    no_coincidencia=[];
    for i=1:length(vec_cel_ind_j)% Bucle que recorre todos los vecinos de la celula bajo estudio
        no_coincidencia(i)=isempty(find(vec_cel_ind_j(i)==celulas_no_validas_previa, 1)); %Variable que vale 1 si no existe coincidencia de la celula j con alguna de las celulas no validas
    end
    if sum(no_coincidencia)==0
        celulas_validas_no_borde_no_vecinos=celulas_validas_previa(celulas_validas_previa~=celulas_validas_previa(j));
        celulas_no_validas=sort([celulas_validas_previa(j);celulas_no_validas_previa]);
        bandera=1;
    end
end
celulas_validas_no_borde_no_vecinos=celulas_validas_no_borde_no_vecinos';
if (bandera==0)
    celulas_validas_no_borde_no_vecinos=celulas_validas_previa;
    celulas_no_validas=celulas_no_validas_previa;
end


% %%%%%%%%%%%      CELULAS_VÁLIDAS      %%%%%%%%%%%%
% 
% %celulas_validas -> (Ni celulas del borde ni sus vecinos) == celulas_validas_no_borde_no_vecinos
% %celulas_validas_no_borde -> (Todas las células excepto las del borde)
% %celulas_no_validas -> (Células del borde y vecinas)
% %celulas_no_validas_no_borde -> (Células vecinas de las del borde)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Señalamos como celulas validas celulas que formen un grafo independiente
%pequeño
grafo_bin=zeros(length(num_celulas),length(num_celulas));
for i=1:size(celulas_validas_no_borde,1)
    numero_celula=celulas_validas_no_borde(i);
    vecinos_network=vecinos_real{numero_celula};

    %pintamos grafica
    for ii=1:(length(vecinos_network))
        grafo_bin(numero_celula,vecinos_network(ii))=1;
        grafo_bin(vecinos_network(ii),numero_celula)=1;
    end
end
%Eliminamos conexiones con celulas no validas
for i=1:(length(num_celulas))
    if isempty(find(i==celulas_validas_no_borde, 1))==1
        grafo_bin(i,:)=0;
        grafo_bin(:,i)=0;
    end    
end

%Cuento numero de grafos existentes
for i=1:(length(num_celulas))
   v=find(grafo_bin(i,:)==1);
   for j=1:length(v)
       grafo_bin(i,:)=(grafo_bin(v(j),:)|grafo_bin(i,:));
       grafo_bin(v(j),:)=(grafo_bin(v(j),:)|grafo_bin(i,:));
   end
end
Suma=sum(grafo_bin,1);
dat=unique(Suma);
grafo_may=dat(1,end);
%celulas_validas_no_borde=find(Suma==grafo_may);%comentar si error
celulas_validas_no_borde=celulas_validas_no_borde';
%celulas_no_validas_no_borde=(1:(length(num_celulas)));%comentar si error
%celulas_no_validas_no_borde(celulas_validas_no_borde)=[];%comentar si error
celulas_no_validas_no_borde=celulas_no_validas_no_borde';
 


%% CARACTERISTICAS 

%%%Extraemos en canal rojo para obtener propiedades de las celulas
%%%musculares

R=Img(:,:,1);
R=imadjust(R);
R(roi~=1)=0;



%% Media de areas y centros, n_celulas, areas de cel rojas y negras
%Area 

s  = regionprops(celulas_roi, 'Area'); 
Area_celula = cat(1, s.Area);

Area_celulas_validas=Area_celula(celulas_validas_no_borde);
Mean_Area=mean(Area_celulas_validas);
Std_Area=std(Area_celulas_validas);


%Multiplicamos las celulas etiquetadas de la roi por la mascara mejorada,
%para que tengan las mismas etiquetas

mascara_mejorada(roi~=1)=0;
masc_mejorada_logica=logical(mascara_mejorada);
BW_mascara_mejorada=masc_mejorada_logica.*celulas_roi;


Mean_R = regionprops(celulas_roi, R, 'MeanIntensity'); %%He cambiado aqui BW_mascara_mejorada por celulas_roi en las 3y4  linea siguiente
Mean_R = cat(1, Mean_R.MeanIntensity);
for i=1:max(max(celulas_roi))
    Intensidades=R(celulas_roi==i);
    n_int=length(find(Intensidades>=100));
    porc_rojo_intenso(i)=n_int/length(Intensidades);
end
Mean_R_validas=Mean_R(celulas_validas_no_borde);
porc_rojo_intenso_validas=porc_rojo_intenso(celulas_validas_no_borde);

% Numero celulas rojas y negras
n_celulas_rojas=length(unique([find(Mean_R_validas>50)' find(porc_rojo_intenso_validas>=0.15)]));
n_celulas_negras=length(Mean_R_validas)-n_celulas_rojas;

Average_celulas_rojas=n_celulas_rojas/(n_celulas_negras+n_celulas_rojas);

%Area roja, Area negra
Mean_Area_rojas=mean(Area_celula(unique([find(Mean_R_validas>50)' find(porc_rojo_intenso_validas>=0.15)])));
Std_Area_rojas=std(Area_celula(Mean_R_validas>50));
Mean_Area_negras=mean(Area_celula(Mean_R_validas<=50));
Std_Area_negras=std(Area_celula(Mean_R_validas<=50));

%elongacion
mayor_eje = regionprops(BW_mascara_mejorada,'MajorAxisLength');
mayor_eje = cat(1, mayor_eje.MajorAxisLength);
Mean_mayor_eje = mean(mayor_eje(celulas_validas_no_borde));

menor_eje = regionprops(BW_mascara_mejorada,'MinorAxisLength');
menor_eje = cat(1, menor_eje.MinorAxisLength);
Mean_menor_eje = mean(menor_eje(celulas_validas_no_borde));

Relacion_ejes=mayor_eje./menor_eje;
Mean_Relacion_ejes=mean(Relacion_ejes(celulas_validas_no_borde));
Std_Relacion_ejes=std(Relacion_ejes(celulas_validas_no_borde));

%convex hull
Pix_region_convexa=regionprops(BW_mascara_mejorada,'Solidity');
Pix_region_convexa = cat(1, Pix_region_convexa.Solidity);
Mean_Pix_region_convexa=mean(Pix_region_convexa(celulas_validas_no_borde));
Std_Pix_region_convexa=std(Pix_region_convexa(celulas_validas_no_borde));

% angulo con respecto al eje x
angulos = regionprops(BW_mascara_mejorada,'Orientation');
angulos = cat(1, angulos.Orientation);
Mean_angulos = mean(angulos(celulas_validas_no_borde));
Std_angulos = std(angulos(celulas_validas_no_borde));

%centros de la imagen entera
centro = regionprops(celulas_roi,'Centroid'); %he modficado celulas_roi por BW_mascara_mejorada
centros = cat(1, centro.Centroid);

%Area water
contorno_water(roi~=1)=0;

L_contorno_water = bwlabel(contorno_water,8);
area_water = regionprops(L_contorno_water, 'Area');
area_water = cat(1, area_water.Area);
[value,ix]=max(area_water);
contorno_water(L_contorno_water==ix)=0;
contorno_water(contorno_water~=0)=1;

L_contorno_water = contorno_water.*celulas_roi;
area_water = regionprops(L_contorno_water, 'Area');
area_water = cat(1, area_water.Area);
mascara_mejorada=im2bw(mascara_mejorada,0.5);

%%Recorro celulas para terminar de limpiar conexiones restantes
for i=1:max(max(L_contorno_water))
    v=find(mascara_mejorada(L_contorno_water==i)==0);
    probabilidad_eliminar=(length(v)/area_water(i))*100;
    if probabilidad_eliminar>80
        contorno_water(L_contorno_water==i)=0;
    end
end
L_contorno_water = bwlabel(contorno_water,8);
L_contorno_water(L_contorno_water>0)=1;
L_contorno_water=L_contorno_water.*celulas_roi;

s  = regionprops(L_contorno_water, 'Area');
Area_water = cat(1, s.Area);


%relacion entre areas celulas y las reales del water
% Area_water_validas = Area_water(celulas_validas_no_borde);
Relacion_areas_validas=Area_water(celulas_validas_no_borde)./Area_celula(celulas_validas_no_borde);

Mean_relacion_areas=mean(Relacion_areas_validas);
Std_relacion_areas=std(Relacion_areas_validas);

%Media de vecinos
n_vecinos_validos=n_vecinos(celulas_validas_no_borde);
Mean_vecinos=mean(n_vecinos_validos);
Std_vecinos=std(n_vecinos_validos);

% represento la network en el contorno en rojo
aux=mascara_mejorada;
aux(aux~=0)=1;
contorno2=contorno.*roi;
aux2=aux+contorno2;


pi_contorno=find(aux2==0);
pi_contorno2=find(contorno2==1);
R_c=ones(H,W);
G_c=ones(H,W);
B_c=ones(H,W);

R_c(pi_contorno)=1;
G_c(pi_contorno)=0;
B_c(pi_contorno)=0;

R_c(pi_contorno2)=0.8;
G_c(pi_contorno2)=0.8;
B_c(pi_contorno2)=0.8;

con_rojo=cat(3,R_c,G_c,B_c);
cont=0;
ind_centros=2;
imshow(con_rojo)

% Creacion de network
i_rojos=0;
i_negros=0;
relacion_area_vecinos=zeros(1,1);
relacion_eje_mayor_vecinos=zeros(1,1);
relacion_eje_menor_vecinos=zeros(1,1);
relacion_relacion_ejes_vecinos=zeros(1,1);
relacion_Pix_region_convexa_vecinos=zeros(1,1);
relacion_Angulos_vecinos=zeros(1,1);
relacion_Relacion_areas_vecinos=zeros(1,1);
grafo_binario=zeros(size(Area_celula,1),size(Area_celula,1));
grafo=zeros(size(Area_celula,1),size(Area_celula,1));
celulas_rojas=zeros(1);
celulas_negras=zeros(1);
n_vecinos_de_rojos=zeros(1);
n_vecinos_rojos_de_rojos=zeros(1);
n_vecinos_negros_de_rojos=zeros(1);
n_vecinos_de_negros=zeros(1);
n_vecinos_rojos_de_negros=zeros(1);
n_vecinos_negros_de_negros=zeros(1);

for i=1:length(celulas_validas_no_borde_no_vecinos)
    numero_celula=celulas_validas_no_borde_no_vecinos(i);
    vecinos_network=vecinos_real{numero_celula};
    
    % calculo el area de cada celula /[ media(area de las celulas vecina) + std(area celulas vecina)]
    relacion_area_vecinos(i)=Area_celula(numero_celula)/(mean(Area_celula(vecinos_network))+std(Area_celula(vecinos_network)));
    relacion_eje_mayor_vecinos(i)=mayor_eje(numero_celula)/(mean(mayor_eje(vecinos_network))+std(mayor_eje(vecinos_network)));
    relacion_eje_menor_vecinos(i)=menor_eje(numero_celula)/(mean(menor_eje(vecinos_network))+std(menor_eje(vecinos_network)));
    relacion_relacion_ejes_vecinos(i)=(mayor_eje(numero_celula)./menor_eje(numero_celula))/(mean((mayor_eje(vecinos_network))./menor_eje(vecinos_network)));
    relacion_Pix_region_convexa_vecinos(i)=Pix_region_convexa(numero_celula)/(mean(Pix_region_convexa(vecinos_network))+std(Pix_region_convexa(vecinos_network)));
    relacion_Angulos_vecinos(i)=angulos(numero_celula)/(mean(angulos(vecinos_network))+std(angulos(vecinos_network)));
    relacion_Relacion_areas_vecinos(i)=Relacion_areas_validas(i)/(mean(Relacion_areas_validas(i))+std(Relacion_areas_validas(i)));
    
end

for i=1:length(celulas_validas_no_borde)
    numero_celula=celulas_validas_no_borde(i);
    vecinos_network=vecinos_real{numero_celula};

        
    % Vecinos de celulas rojas y negras
    if Mean_R(numero_celula)>100 || porc_rojo_intenso(numero_celula)>=0.15
        hold on,
        plot(round(centros(numero_celula,1)),round(centros(numero_celula,2)),'.r','MarkerSize',20)
        
        i_rojos=i_rojos+1;
        celulas_rojas(i_rojos)=numero_celula;
        n_vecinos_de_rojos(i_rojos)=length(vecinos_real{numero_celula}); % vecinos de celulas rojas
        if ~isempty(vecinos_real{numero_celula})
            n_vecinos_rojos_de_rojos(i_rojos)=length(union(find(Mean_R(vecinos_real{numero_celula})>100), find(porc_rojo_intenso(vecinos_real{numero_celula})>=0.15)));
            n_vecinos_negros_de_rojos(i_rojos)=length(intersect(find(Mean_R(vecinos_real{numero_celula})<=100), find(porc_rojo_intenso(vecinos_real{numero_celula})<0.15)));
        end
    else
        
        hold on,
        plot(round(centros(numero_celula,1)),round(centros(numero_celula,2)),'.g','MarkerSize',20)
        
        i_negros=i_negros+1;
        celulas_negras(i_negros)=numero_celula;
        n_vecinos_de_negros(i_negros)=length(vecinos_real{numero_celula});
        if ~isempty(vecinos_real{numero_celula})
            n_vecinos_rojos_de_negros(i_negros)=length(union(find(Mean_R(vecinos_real{numero_celula})>100), find(porc_rojo_intenso(vecinos_real{numero_celula})>=0.15)));
            n_vecinos_negros_de_negros(i_negros)=length(intersect(find(Mean_R(vecinos_real{numero_celula})<=100), find(porc_rojo_intenso(vecinos_real{numero_celula})<0.15)));
        end
    end
        
    %pintamos grafica
    for ii=1:(length(vecinos_network))
        grafo_binario(numero_celula,vecinos_network(ii))=1;
        grafo_binario(vecinos_network(ii),numero_celula)=1;
    end
    
    for ii=1:(length(vecinos_network))
        grafo(numero_celula,vecinos_network(ii))=sqrt((centros(numero_celula,1)-centros(vecinos_network(ii),1))^2+(centros(numero_celula,2)-centros(vecinos_network(ii),2))^2);
        grafo(vecinos_network(ii),numero_celula)=sqrt((centros(numero_celula,1)-centros(vecinos_network(ii),1))^2+(centros(numero_celula,2)-centros(vecinos_network(ii),2))^2);
    end

    
    for k=1:length(vecinos_network)
         hold on,
        plot([round(centros(numero_celula,1)),round(centros(vecinos_network(k),1))], [round(centros(numero_celula,2)),round(centros(vecinos_network(k),2))],'Color','black')
    end
    
    
end

cd ..\..\..\..\..
cd Caracteristicas_extraidas

if isdir(['8 Caracteristicas_extraidas ' fecha(end-7:end)])~=1
    mkdir(['8 Caracteristicas_extraidas ' fecha(end-7:end)])
end

cd (['8 Caracteristicas_extraidas ' fecha(end-7:end)])

if isdir([nombre])~=1
    mkdir([nombre])
end

cd (nombre)

stringres=strcat(nombre,'_network_',num2str(n_roi),'.bmp');
print('-f1','-dbmp',stringres)

%Caracteristicas vecinos rojas\negras
if i_rojos>0
    Mean_vecinos_de_rojos=mean(n_vecinos_de_rojos);
    Std_vecinos_de_rojos=std(n_vecinos_de_rojos);
    Mean_vecinos_rojos_de_rojos=mean(n_vecinos_rojos_de_rojos);
    Std_vecinos_rojos_de_rojos=std(n_vecinos_rojos_de_rojos);
    Mean_vecinos_negros_de_rojos=mean(n_vecinos_negros_de_rojos);
    Std_vecinos_negros_de_rojos=std(n_vecinos_negros_de_rojos);
else
    Mean_vecinos_de_rojos=0;
    Std_vecinos_de_rojos=0;
    Mean_vecinos_rojos_de_rojos=0;
    Std_vecinos_rojos_de_rojos=0;
    Mean_vecinos_negros_de_rojos=0;
    Std_vecinos_negros_de_rojos=0;
end

if i_negros>0
    Mean_vecinos_de_negros=mean(n_vecinos_de_negros);
    Std_vecinos_de_negros=std(n_vecinos_de_negros);
    Mean_vecinos_rojos_de_negros=mean(n_vecinos_rojos_de_negros);
    Std_vecinos_rojos_de_negros=std(n_vecinos_rojos_de_negros);
    Mean_vecinos_negros_de_negros=mean(n_vecinos_negros_de_negros);
    Std_vecinos_negros_de_negros=std(n_vecinos_negros_de_negros);
else
    Mean_vecinos_de_negros=0;
    Std_vecinos_de_negros=0;
    Mean_vecinos_rojos_de_negros=0;
    Std_vecinos_rojos_de_negros=0;
    Mean_vecinos_negros_de_negros=0;
    Std_vecinos_negros_de_negros=0;
end

% caracteristicas de network
grafo_rect=grafo(celulas_validas_no_borde,celulas_validas_no_borde);

cd ..\..\..
cd 'Codigo\Codigo_BCT'

[n_conexiones_cada_nodos, Suma_pesos_cada_nodo, Correlacion_entre_grados_nodos, Densidad_conexiones, Coef_cluster, T, estructura_optima,modularidad_maximizada, Matriz_distancias_mas_cortas_de_todos_nodos,lambda,efficiency,ecc,radius,diameter, BC]=Prueba_brain(grafo,grafo_binario,grafo_rect,celulas_validas_no_borde);
cd ..


Suma_pesos=Suma_pesos_cada_nodo(celulas_validas_no_borde);

Mean_Area;
Std_Area;
Mean_Area_rojas;
Std_Area_rojas;
Mean_Area_negras;
Std_Area_negras;
Mean_mayor_eje;
Mean_menor_eje;
Mean_Relacion_ejes;
Std_Relacion_ejes;
Mean_Pix_region_convexa;
Std_Pix_region_convexa;
Mean_angulos;
Std_angulos;
Mean_relacion_areas;
Std_relacion_areas;
Mean_vecinos;
Std_vecinos;
Std_vecinos_de_rojos;
Std_vecinos_de_negros;
Mean_vecinos_rojos_de_rojos;
Mean_vecinos_negros_de_rojos;
Mean_vecinos_rojos_de_negros;
Mean_vecinos_negros_de_negros;


Mean_Relacion_areas_vecindad=mean(relacion_area_vecinos);
Std_Relacion_areas_vecindad=std(relacion_area_vecinos);
Mean_relacion_eje_mayor_vecinos=mean(relacion_eje_mayor_vecinos);
Std_relacion_eje_mayor_vecinos=std(relacion_eje_mayor_vecinos);
Mean_relacion_eje_menor_vecinos=mean(relacion_eje_menor_vecinos);
Std_relacion_eje_menor_vecinos=std(relacion_eje_menor_vecinos);
Mean_relacion_relacion_ejes_vecinos=mean(relacion_relacion_ejes_vecinos);
Std_relacion_relacion_ejes_vecinos=std(relacion_relacion_ejes_vecinos);
Mean_relacion_Pix_region_convexa_vecinos=mean(relacion_Pix_region_convexa_vecinos);
Std_relacion_Pix_region_convexa_vecinos=std(relacion_Pix_region_convexa_vecinos);
Mean_relacion_Angulos_vecinos=mean(relacion_Angulos_vecinos);
Std_relacion_Angulos_vecinos=std(relacion_Angulos_vecinos);
Mean_relacion_Relacion_areas_vecinos=mean(relacion_Relacion_areas_vecinos);
Std_relacion_Relacion_areas_vecinos=std(relacion_Relacion_areas_vecinos);


Mean_suma_pesos=mean(Suma_pesos);
Desv_suma_pesos=std(Suma_pesos);

if i_negros==0
    
    Mean_pesos_celulas_negras=0;
    Desv_pesos_celulas_negras=0;
    Mean_Coef_cluster_negras=0;
    Desv_Coef_cluster_negras=0;
    Mean_excentricidad_negras=0;
    Desv_excentricidad_negras=0;
    Mean_BC_negras=0;
    Desv_BC_negras=0;
    
else
    
    Mean_pesos_celulas_negras=mean(Suma_pesos_cada_nodo((celulas_negras)));
    Desv_pesos_celulas_negras=std(Suma_pesos_cada_nodo((celulas_negras)));
    Mean_Coef_cluster_negras=mean(Coef_cluster(celulas_negras));
    Desv_Coef_cluster_negras=std(Coef_cluster(celulas_negras));
    Mean_excentricidad_negras=mean(ecc(celulas_negras));
    Desv_excentricidad_negras=std(ecc(celulas_negras));
    Mean_BC_negras=mean(BC(celulas_negras));
    Desv_BC_negras=std(BC(celulas_negras));
    
    
end



if i_rojos==0
    Mean_pesos_celulas_rojas=0;
    Desv_pesos_celulas_rojas=0;
    Mean_Coef_cluster_rojas=0;
    Desv_Coef_cluster_rojas=0;
    Mean_excentricidad_rojas=0;
    Desv_excentricidad_rojas=0;
    Mean_BC_rojas=0;
    Desv_BC_rojas=0;
else
    
    Mean_pesos_celulas_rojas=mean(Suma_pesos_cada_nodo((celulas_rojas)));
    Desv_pesos_celulas_rojas=std(Suma_pesos_cada_nodo((celulas_rojas)));
    Mean_Coef_cluster_rojas=mean(Coef_cluster(celulas_rojas));
    Desv_Coef_cluster_rojas=std(Coef_cluster(celulas_rojas));
    Mean_excentricidad_rojas=mean(ecc(celulas_rojas));
    Desv_excentricidad_rojas=std(ecc(celulas_rojas));
    Mean_BC_rojas=mean(BC(celulas_rojas));
    Desv_BC_rojas=std(BC(celulas_rojas));
    
end



Mean_Coef_cluster=mean(Coef_cluster);
Desv_Coef_cluster=std(Coef_cluster);
Mean_excentricidad=mean(ecc(celulas_validas_no_borde));
Desv_excentricidad=std(ecc(celulas_validas_no_borde));
Mean_BC=mean(BC(celulas_validas_no_borde));
Desv_BC=std(BC(celulas_validas_no_borde));



M=triu(Matriz_distancias_mas_cortas_de_todos_nodos(celulas_validas_no_borde,celulas_validas_no_borde));% me quedo solo con los valores de una celula a otra y no viceversa, ya que es una matriz simetriica
Mean_dist=mean(M(M~=0));
Desv_dist=std(M(M~=0));


if i_negros~=0 
    M1=triu(Matriz_distancias_mas_cortas_de_todos_nodos(celulas_negras,celulas_negras));% me quedo solo con los valores de una celula a otra y no viceversa, ya que es una matriz simetriica
    Mean_dist_negras_negras=mean(M1(M1~=0));
    Desv_dist_negras_negras=std(M1(M1~=0));
else
    Mean_dist_negras_negras=0;
    Desv_dist_negras_negras=0;
end

if i_rojos~=0
    M3=triu(Matriz_distancias_mas_cortas_de_todos_nodos(celulas_rojas,celulas_rojas));% me quedo solo con los valores de una celula a otra y no viceversa, ya que es una matriz simetriica
    Mean_dist_rojas_rojas=mean(M3(M3~=0));
    Desv_dist_rojas_rojas=std(M3(M3~=0));
else
    Mean_dist_rojas_rojas=0;
    Desv_dist_rojas_rojas=0;
end

if i_rojos~=0 && i_negros~=0
    M2=triu(Matriz_distancias_mas_cortas_de_todos_nodos(celulas_negras,celulas_rojas));% me quedo solo con los valores de una celula a otra y no viceversa, ya que es una matriz simetriica
    Mean_dist_negras_rojas=mean(M2(M2~=0));
    Desv_dist_negras_rojas=std(M2(M2~=0));
    
    M4=triu(Matriz_distancias_mas_cortas_de_todos_nodos(celulas_rojas,celulas_negras));% me quedo solo con los valores de una celula a otra y no viceversa, ya que es una matriz simetriica
    Mean_dist_rojas_negras=mean(M4(M4~=0));
    Desv_dist_rojas_negras=std(M4(M4~=0));
    
else 
    Mean_dist_negras_rojas=0;
    Desv_dist_negras_rojas=0;
    
    Mean_dist_rojas_negras=0;
    Desv_dist_rojas_negras=0;
end



radius; %minimo ecc
diameter; %maximo ecc
efficiency;  %

cd Gergana_Bounova
[prs, a, s]=Prueba_Bounova(grafo,grafo_binario,grafo_rect,celulas_validas_no_borde);
cd ..\..
Coef_pearson=prs;
conectividad_algebraica=a;
Metrica_s=s;

Assortattivity=Correlacion_entre_grados_nodos;
Densidad_conexiones;
Transitivity=T;
modularidad_maximizada;





% % %% Extraccion de caracteristicas del plano azul
% % %Numero de nucleos en celula water %cc 1
% % 
% % %Para seleccionar celulas validas ordenamos cells_water segun cells
% % picos_centroide_en_cell_water=datos_celulas_water(cell2mat(datos_celulas_water(:,2)),6); %%nº de centroides de picos dentro de las celulas water (no colágeno), incluidas no validas
% % picos_centroide_en_cell_water_validos=picos_centroide_en_cell_water(celulas_validas_no_borde,1); %nº de nucleos dentro de celulas validas
% % cantidad_picos_cell_water=zeros(size(picos_centroide_en_cell_water_validos,1),1);
% % for i=1:size(picos_centroide_en_cell_water_validos,1)
% %     cantidad_picos_cell_water(i)=size(picos_centroide_en_cell_water_validos{i,1},1);
% % end
% % 
% % Picos_centr_cell_water=sum(cantidad_picos_cell_water);
% % 
% % 
% % media_picos_en_cel_water_val=Picos_centr_cell_water/length(celulas_validas_no_borde);
% % 
% % %Numero de nucleos en el colágeno perteneciente a la celula %cc 2
% % picos_centroide_en_colageno=datos_celulas_menos_water(cell2mat(datos_celulas_menos_water(:,2)),6);
% % picos_centroide_en_colageno_valido=picos_centroide_en_colageno(celulas_validas_no_borde,1);
% % cantidad_picos_colageno =zeros(size(picos_centroide_en_colageno_valido,1),1);
% % for i=1:size(picos_centroide_en_colageno_valido,1)
% %     cantidad_picos_colageno(i)=size(picos_centroide_en_colageno_valido{i,1},1);
% % end
% % 
% % Picos_centr_colageno=sum(cantidad_picos_colageno);
% % media_picos_en_colageno_val=Picos_centr_colageno/length(celulas_validas_no_borde);
% % 
% % 
% % %Numero de nucleos en celula total (c.water + colágeno) %cc 3
% % picos_centroide_en_cell_total=datos_celulas(:,6);
% % picos_centroide_en_cell_total_validas=picos_centroide_en_cell_total(celulas_validas_no_borde,1);
% % cantidad_picos_validos=zeros(size(picos_centroide_en_cell_total_validas,1),1);
% % for i=1:size(picos_centroide_en_cell_total_validas,1)
% %     cantidad_picos_validos(i)=size(picos_centroide_en_cell_total_validas{i,1},1);
% % end
% % 
% % Picos_centr_cell_total=sum(cantidad_picos_validos);
% % media_picos_celula_total_val=Picos_centr_cell_total/length(celulas_validas_no_borde);
% % 
% % 
% % 
% % %Desviación estándar del numero de picos en celulas water %%cc4
% % 
% % desv_n_picos_cwater=std(cantidad_picos_cell_water);
% % 
% % %Desviación estándar del numero de picos en colageno %%cc5
% % 
% % desv_n_picos_colageno=std(cantidad_picos_colageno);
% % 
% % %Desviación estándar del numero de picos en celulas+colágeno %%cc6
% % 
% % desv_n_picos_cel_total=std(cantidad_picos_validos);
% % 
% % 
% % 
% % % %porcentaje nodo en area water %cc 4
% % % %Para seleccionar celulas validas ordenamos cells_water segun cells
% % % 
% % % area_nucleo_cell_water=datos_celulas_water(cell2mat(datos_celulas_water(:,2)),8);
% % % area_nucleo_cell_water_validas=cell2mat(area_nucleo_cell_water(celulas_validas_no_borde,1));
% % % area_cell_water=cell2mat(datos_celulas_water(cell2mat(datos_celulas_water(:,2)),7));
% % % area_cell_water_validas=area_cell_water(celulas_validas_no_borde,1);
% % % area_total_nucleo_water=sum(area_nucleo_cell_water_validas);
% % % area_total_water=sum(area_cell_water_validas);
% % % 
% % % Porcentaje_Area_Nucleos_water_cell=(area_total_nucleo_water/area_total_water)*100;
% % % 
% % % %porcentaje nodo en area cell menos water %cc 5
% % % area_nucleo_cell_menos_water=datos_celulas_menos_water(cell2mat(datos_celulas_menos_water(:,2)),8);
% % % area_nucleo_cell_menos_water_validas=cell2mat(area_nucleo_cell_menos_water(celulas_validas_no_borde,1));
% % % area_cell_menos_water=cell2mat(datos_celulas_menos_water(cell2mat(datos_celulas_menos_water(:,2)),7));
% % % area_cell_menos_water_validas=area_cell_menos_water(celulas_validas_no_borde,1);
% % % area_total_nucleo_menos_water=sum(area_nucleo_cell_menos_water_validas);
% % % area_total_menos_water=sum(area_cell_menos_water_validas);
% % % 
% % % Porcentaje_Area_Nucleos_menos_water_cell=(area_total_nucleo_menos_water/area_total_menos_water)*100;
% % % 
% % % %porcentaje nodo en area total %cc 6
% % % area_nucleo_cell=datos_celulas(:,8);
% % % area_nucleo_cell_validas=cell2mat(area_nucleo_cell(celulas_validas_no_borde,1));
% % % area_cell=cell2mat(datos_celulas(:,7));
% % % area_cell_validas=area_cell(celulas_validas_no_borde,1);
% % % area_total_nucleo=sum(area_nucleo_cell_validas);
% % % area_total=sum(area_cell_validas);
% % % 
% % % Porcentaje_Area_Nucleos_cell=(area_total_nucleo/area_total)*100;
% % 
% % %porcentaje objeto en area water %cc 7
% % area_objeto_cell_water=datos_celulas_water(cell2mat(datos_celulas_water(:,2)),10);
% % area_objeto_cell_water_validas=cell2mat(area_objeto_cell_water(celulas_validas_no_borde,1));
% % area_cell_water=datos_celulas_water(cell2mat(datos_celulas_water(:,2)),7);
% % area_cell_water_validas=cell2mat(area_cell_water(celulas_validas_no_borde,1));
% % area_total_objeto_water=sum(area_objeto_cell_water_validas);
% % area_total_water=sum(area_cell_water_validas);
% % 
% % Porcentaje_Area_Objeto_water_cell=(area_total_objeto_water/area_total_water)*100;
% % 
% % %porcentaje objeto en colagenos %cc 8
% % area_objeto_colageno=datos_celulas_menos_water(cell2mat(datos_celulas_menos_water(:,2)),10);
% % area_objeto_colageno_valido=cell2mat(area_objeto_colageno(celulas_validas_no_borde,1));
% % area_colageno=datos_celulas_menos_water(cell2mat(datos_celulas_menos_water(:,2)),7);
% % area_colageno_valido=cell2mat(area_colageno(celulas_validas_no_borde,1));
% % area_total_objeto_colageno=sum(area_objeto_colageno_valido);
% % area_total_colageno=sum(area_colageno_valido);
% % 
% % Porcentaje_Area_Objeto_colageno=(area_total_objeto_colageno/area_total_colageno)*100;
% % 
% % %porcentaje objeto en area cell %cc 9
% % area_objeto_cell=datos_celulas(:,10);
% % area_objeto_cell_validas=cell2mat(area_objeto_cell(celulas_validas_no_borde,1));
% % area_cell=datos_celulas(:,7);
% % area_cell_validas=cell2mat(area_cell(celulas_validas_no_borde,1));
% % area_total_objeto=sum(area_objeto_cell_validas);
% % area_total=sum(area_cell_validas);
% % 
% % Porcentaje_Area_Objeto_cell=(area_total_objeto/area_total)*100;
% % 
% % 
% % 
% % %Desviación estándar del area de objetos en celulas water %%cc10
% % 
% % desv_area_obj_cwater=std(area_objeto_cell_water_validas);
% % 
% % %Desviación estándar del area de objetos en colageno %%cc11
% % 
% % desv_area_obj_colageno=std(area_objeto_colageno_valido);
% % 
% % %Desviación estándar del area de objetos en celulas+colágeno %%cc12
% % 
% % desv_area_obj_cel_total=std(area_objeto_cell_validas);
% 
% 
% 
% % %Nucleos cuyo porcentaje caiga 90% fuera de AREA WATER %cc 10
% % PORCENTAJE_INCURSION_NODOS_EN_WATER_CELL=datos_objetos_pico(:,12); %porcentaje de n_pico en celula original
% % NUM_CELL_PICO=datos_objetos_pico(:,2); %nº de nucleos pico existentes
% % Nucleo_90_out_water=0;
% % for i=1:size(PORCENTAJE_INCURSION_NODOS_EN_WATER_CELL,1)
% %     if i==33
% %         k=0;
% %     end
% %     pic=NUM_CELL_PICO{i,1};
% %     Dato=PORCENTAJE_INCURSION_NODOS_EN_WATER_CELL{i,1};
% %     aux=pic;
% %     iii=1;
% %     v=zeros(1);
% %     for j=1:length(pic)
% %        if isempty(find(pic(j)==celulas_validas_no_borde, 1))==1
% %            cantidad_borrada=j-iii+1;
% %             aux(cantidad_borrada)=[];
% %             v(iii)=pic(j);
% %             iii=iii+1;
% %        end
% %     end
% %     NUM_CELL_PICO{i,1}=aux;
% %     for k=1:length(v)
% %         regulacion=cell2mat(datos_celulas_water(:,2));
% %         picos_ordenados=regulacion(datos_objetos_pico{i,3});
% %        if isempty(v(k)==picos_ordenados)==0
% %            [~,ix]=find(v(k)==picos_ordenados);
% %            Dato(ix)=[];
% %        end
% %     end
% %     num_pico_cell=length(NUM_CELL_PICO{i,1});
% %     cell_fuera=length(find(Dato>=10));
% %     Nucleo_90_out_water=Nucleo_90_out_water+(num_pico_cell-cell_fuera);
% % end
% 
% 
% % %Numero de núcleos contenidos en celula water por célula válida % cc.11
% % Average_nucleos_water_cell=Picos_centr_cell_water/length(celulas_validas_no_borde);
% % 
% % %Numero de núcleos contenido en colágeno por célula válida % cc.12
% % Average_nucleos_colageno_cell=Nucleos_centr_cell_menos_water/length(celulas_validas_no_borde);
% % 
% % %Numero de núcleos por célula válida % cc.13
% % Average_nucleos_cell=Nucleos_centr_cell/length(celulas_validas_no_borde);
% % 
% % 
% % 
% % 
% % imshow(con_rojo)
% % grafo_binario_2=zeros(size(Area_celula,1),size(Area_celula,1));
% % grafo_2=zeros(size(Area_celula,1),size(Area_celula,1));
% % for i=1:size(celulas_validas_no_borde,1)
% %     numero_celula=celulas_validas_no_borde(i);
% %     vecinos_network=vecinos_real{numero_celula};
% %     cuenta=length(vecinos_network);
% %     k=1;
% %    while k<=cuenta
% %         if isempty(find(vecinos_network(k)==(celulas_validas_no_borde), 1))==1
% %             vecinos_network(k)=[];
% %             k=k-1;
% %             cuenta=cuenta-1;
% %         end
% %         k=k+1;
% %     end
% % 
% %     %pintamos grafica
% %     for ii=1:(length(vecinos_network))
% %         grafo_binario_2(numero_celula,vecinos_network(ii))=1;
% %         grafo_binario_2(vecinos_network(ii),numero_celula)=1;
% %     end
% %     
% %     for ii=1:(length(vecinos_network))
% %         grafo_2(numero_celula,vecinos_network(ii))=(1-Relacion_areas(numero_celula))+(1-Relacion_areas(vecinos_network(ii)));
% %         grafo_2(vecinos_network(ii),numero_celula)=(1-Relacion_areas(numero_celula))+(1-Relacion_areas(vecinos_network(ii)));
% %     end
% % 
% %     hold on,
% %     plot(round(centros(numero_celula,1)),round(centros(numero_celula,2)),'.g','MarkerSize',15)
% %     for k=1:length(vecinos_network)
% %          hold on,
% %         plot([round(centros(numero_celula,1)),round(centros(vecinos_network(k),1))], [round(centros(numero_celula,2)),round(centros(vecinos_network(k),2))],'Color','black')
% %     end    
% % end
% % 
% % 
% % stringres=strcat('Imagenes_procesadas\',fecha,'\',nombre,'\',nombre,'-network_colageno.bmp');
% % print('-f1','-dbmp',stringres)
% % 
% % 
% % grafo_rect_2=grafo_2(celulas_validas_no_borde,celulas_validas_no_borde);
% % cd código
% % cd Codigo_BCT
% % [n_conexiones_cada_nodos_2, Suma_pesos_cada_nodo_2, Correlacion_entre_grados_nodos_2, Densidad_conexiones_2, Coef_cluster_2, T_2, estructura_optima_2,modularidad_maximizada_2, Matriz_distancias_mas_cortas_de_todos_nodos_2,lambda_2,efficiency_2,ecc_2,radius_2,diameter_2, BC_2]=Prueba_brain(grafo_2,grafo_binario_2,grafo_rect_2,celulas_validas_no_borde);
% % cd ..;
% % 
% % 
% % Suma_pesos_grafo_2=Suma_pesos_cada_nodo_2(celulas_validas_no_borde);
% % Mean_suma_pesos_grafo_2=mean(Suma_pesos_grafo_2);
% % Desv_suma_pesos_grafo_2=std(Suma_pesos_grafo_2);
% % 
% % Mean_Coef_cluster_grafo_2=mean(Coef_cluster_2);
% % Desv_Coef_cluster_grafo_2=std(Coef_cluster_2);
% % 
% % Mean_excentricidad_grafo_2=mean(ecc_2(celulas_validas_no_borde));
% % Desv_excentricidad_grafo_2=std(ecc_2(celulas_validas_no_borde));
% % 
% % Mean_BC_grafo_2=mean(BC_2(celulas_validas_no_borde));
% % Desv_BC_grafo_2=std(BC_2(celulas_validas_no_borde));
% % 
% % M2=triu(Matriz_distancias_mas_cortas_de_todos_nodos_2(celulas_validas_no_borde,celulas_validas_no_borde));% me quedo solo con los valores de una celula a otra y no viceversa, ya que es una matriz simetriica
% % Mean_dist_grafo_2=mean(M2(M2~=0));
% % Desv_dist_grafo_2=std(M2(M2~=0));
% % 
% % radius_grafo_2=radius_2; %minimo ecc
% % diameter_grafo_2=diameter_2; %maximo ecc
% % efficiency_grafo_2=efficiency_2;  %
% % 
% % cd Gergana_Bounova
% % [prs_2, a_2, s_2]=Prueba_Bounova(grafo_2,grafo_binario_2,grafo_rect_2,celulas_validas_no_borde);
% % cd ..
% % Coef_pearson_grafo_2=prs_2;
% % conectividad_algebraica_grafo_2=a_2;
% % Metrica_s_grafo_2=s_2;
% % 
% % Assortattivity_grafo_2=Correlacion_entre_grados_nodos_2;
% % Densidad_conexiones_grafo_2=Densidad_conexiones_2;
% % Transitivity_grafo_2=T_2;
% % modularidad_maximizada_grafo_2=modularidad_maximizada_2;


cd Caracteristicas_extraidas
cd (['8 Caracteristicas_extraidas ' fecha(end-7:end)])
cd (nombre)

%%Salvamos las 69 cc de la teselación

stringres=['Results_69_cc_roi_',num2str(n_roi),'.mat'];
save (stringres, 'Mean_Area', 'Std_Area', 'Mean_Area_rojas','Std_Area_rojas', 'Mean_Area_negras', 'Std_Area_negras', 'Mean_mayor_eje', 'Mean_menor_eje','Mean_Relacion_ejes','Std_Relacion_ejes','Mean_Pix_region_convexa','Std_Pix_region_convexa','Mean_relacion_areas','Std_relacion_areas','Mean_vecinos','Std_vecinos','Std_vecinos_de_rojos','Std_vecinos_de_negros','Mean_vecinos_rojos_de_rojos','Mean_vecinos_negros_de_rojos','Mean_vecinos_rojos_de_negros','Mean_vecinos_negros_de_negros','Mean_Relacion_areas_vecindad','Std_Relacion_areas_vecindad', 'Mean_relacion_eje_mayor_vecinos','Std_relacion_eje_mayor_vecinos','Mean_relacion_eje_menor_vecinos','Std_relacion_eje_menor_vecinos','Mean_relacion_relacion_ejes_vecinos','Std_relacion_relacion_ejes_vecinos','Mean_relacion_Pix_region_convexa_vecinos','Std_relacion_Pix_region_convexa_vecinos','Mean_relacion_Angulos_vecinos','Std_relacion_Angulos_vecinos','Mean_suma_pesos','Desv_suma_pesos','Mean_pesos_celulas_negras','Desv_pesos_celulas_negras','Mean_pesos_celulas_rojas','Desv_pesos_celulas_rojas','Mean_Coef_cluster','Desv_Coef_cluster','Mean_Coef_cluster_negras','Desv_Coef_cluster_negras','Mean_Coef_cluster_rojas','Desv_Coef_cluster_rojas','Mean_excentricidad','Desv_excentricidad','Mean_excentricidad_negras','Desv_excentricidad_negras','Mean_excentricidad_rojas','Desv_excentricidad_rojas','Mean_BC','Desv_BC','Mean_BC_negras','Desv_BC_negras','Mean_BC_rojas','Desv_BC_rojas','Mean_dist','Desv_dist','Mean_dist_negras_negras','Desv_dist_negras_negras','Mean_dist_negras_rojas','Desv_dist_negras_rojas','Mean_dist_rojas_rojas','Desv_dist_rojas_rojas','Mean_dist_rojas_negras','Desv_dist_rojas_negras','Average_celulas_rojas')
imwrite(Img,[nombre,'.bmp'])
imwrite(Img2,[nombre,'_rois.bmp'])

%%%%% CC NO GUARDADAS EN ESTE CASO %%%%%

%%'Mean_relacion_Relacion_areas_vecinos',
%%'Std_relacion_Relacion_areas_vecinos',                         Relaciona las relaciones de areas solapadas por las celulas water y celulas
% 'Mean_angulos','Std_angulos'
%%'radius','diameter','efficiency','Coef_pearson','conectividad_algebraica',
% 'Metrica_s','Assortattivity','Densidad_conexiones','Transitivity','modularidad_maximizada'

%%%%%------------------------------%%%%%
% 
% 
% %%Guardamos las 12 características de los nucleos
% 
% % stringres=strcat('Imagenes_procesadas\',fecha,'\',nombre,'\Imagenes_necesarias\', 'Nucleos_cc.mat');
% % save (stringres,'media_picos_en_cel_water_val','media_picos_en_colageno_val','media_picos_celula_total_val','desv_n_picos_cwater','desv_n_picos_colageno','desv_n_picos_cel_total','Porcentaje_Area_Objeto_water_cell','Porcentaje_Area_Objeto_colageno','Porcentaje_Area_Objeto_cell','desv_area_obj_cwater','desv_area_obj_colageno','desv_area_obj_cel_total')
% % 
% 
% 
% %,'Mean_suma_pesos_grafo_2','Desv_suma_pesos_grafo_2','Mean_Coef_cluster_grafo_2','Desv_Coef_cluster_grafo_2','Mean_excentricidad_grafo_2','Desv_excentricidad_grafo_2','Mean_BC_grafo_2','Desv_BC_grafo_2','Mean_dist_grafo_2','Desv_dist_grafo_2','radius_grafo_2','diameter_grafo_2','efficiency_grafo_2','Coef_pearson_grafo_2','conectividad_algebraica_grafo_2','Metrica_s_grafo_2','Assortattivity_grafo_2','Densidad_conexiones_grafo_2','Transitivity_grafo_2','modularidad_maximizada_grafo_2'
% 
% stringres=strcat('Imagenes_procesadas\SOD1 sin DAPI\',fecha,'\',nombre,'\Imagenes_necesarias\','Datos_cc.mat');
% save (stringres, 'celulas','celulas_validas_no_borde_no_vecinos','celulas_validas_no_borde','celulas_negras','celulas_rojas','Area_celula', 'mayor_eje', 'Relacion_ejes','Pix_region_convexa','angulos','Relacion_areas','n_vecinos','vecinos_real','Suma_pesos','Coef_cluster','BC','ecc','grafo','grafo_binario')
% 
% 
% 
% %'Suma_pesos_grafo_2','Coef_cluster_2','BC_2','ecc_2','grafo_2','grafo_binario_2'



cd ..\..\..
cd 'Codigo'

end