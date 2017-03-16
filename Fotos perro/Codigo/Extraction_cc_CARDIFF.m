function Extraction_cc(date,name) %We have classify photos by date and name, it is only to load segmented images and save features.

%% PREPARING NECESSARIES IMAGES TO EXTRACT FEATURES
%load perfectly segmented image (You should specify images path)
cellular_mask=load('Img_segmented.jpg');

%get image with cells outline
contour=load('Img_contour.jpg');

%load image with segmented border cells
contour_watershed=load('Img_contour_watershed.jpg');

%define invested contour as cells and we label each cell with a number
cells=1-contour;
cells=bwlabel(cells);


%% GETTING NEIGHBORS

%define image dimensions and save cells in an auxiliar parameter to operate with him 
[H,W]=size(cells);
aux=cells;

%define a radius, always bigger than your contour cell width
radius=4;
%get the number of labelled cells
num_cells=sort(unique(aux));
%discart label=0 (this label represent outlines)
num_cells=num_cells(2:end);

%define a vector to save the number of neighbour for each cell
n_neighbors=zeros(length(num_cells),1);
%define a type cell to save labelled neighbor cells for each cell
neighbors_real=cell(1);

%With this loop for we are going to get the number and the specifics
%neighbor for each labelled cell.


for i=1 : length(num_cells)
    
    %Find 'i' cell region positions
    BW2 = bwperim(aux==i);
    [pi,pj]=find(BW2==1);
    
    %We dilate 'i' cell region in 'radius' pixels
    se = strel('disk',radius);
    BW2_dilate=imdilate(BW2,se);
    
    %Compare labelled than occupy dilated region
    pixels_neighbors=find(BW2_dilate==1);
    neighbors=unique(aux(pixels_neighbors));
    
    %Labels than are different to 0 and 'i' are the neighbors of 'i'.
    neighbors_real{i}=neighbors(neighbors ~= 0 & neighbors ~= i);
    n_neighbors(i)=length(neighbors_real{i});
end

%% VALID CELLS DETECTION (no consider border image cells (to calculate all features) and their neighbour (to calculate some characteristics))

%Mark border with 1 to connect all border cells.
aux(1:end,1)=1;
aux(1:end,W)=1;
aux(1,1:end)=1;
aux(H,1:end)=1;

%aux2 has all border cells labelled as 1, and we save this positions in mask 
aux2=bwlabel(aux,8);
mask=zeros(H,W);
mask(aux2==1)=1;

%Mark no valid cells -> border cells. When we multiply mask and aux, we
%only obtain labels of border cells. Despising valid cells.
aux3=mask.*aux;
no_valid_cells=sort(unique(aux3));
no_valid_cells=no_valid_cells(2:end);

%Calculate valid cells no border, despising no_valid_cells
valid_cells=(1:length(num_cells));
for i=1:length(valid_cells)
    valid_cells(valid_cells==no_valid_cells(i))=[];
end
valid_cells=valid_cells';
valid_cells_no_border=valid_cells;


%Mark no valid cells -> neighbors border cells.
no_valid_cells_2=[];
for i=1:size(no_valid_cells,1)
    no_valid_cells_2=[no_valid_cells_2;neighbors_real{no_valid_cells(i)}];
end
no_valid_cells_2=sort(unique(no_valid_cells_2));

for i=1:length(no_valid_cells_2)
    valid_cells(valid_cells==no_valid_cells_2(i))=[];
end
valid_cells=valid_cells';
valid_cells_no_border_no_neighbors=valid_cells;


%Delete isolated cells (valid cell surrounded by no valid cells) from valid_cells_no_border
    flag=0;
    valid_cells_previous=valid_cells_no_border;
    no_valid_cells_previous=no_valid_cells;
    %analyze all valid cells
    for j=1:length(valid_cells_previous)
        neighbors_j=neighbors_real{valid_cells_previous(j)};
        no_coincidence=[];
        %analyze if all neighbors are all no valid cells or not.
        for i=1:length(neighbors_j)
            no_coincidence(i)=isempty(find(neighbors_j(i)==no_valid_cells_previous, 1));
        end
        %if all neighbors are no valid cells, discard j as valid cell.
        if sum(no_coincidence)==0
            valid_cells_no_border=valid_cells_previous(valid_cells_previous~=valid_cells_previous(j));
            no_valid_cells_no_border=sort([valid_cells_previous(j);no_valid_cells_previous]);
            flag=1;
        end
    end
    valid_cells_no_border=valid_cells_no_border';

    %if there aren't any isolated cell, no change.
    if (flag==0)
        valid_cells_no_border=valid_cells_previous;
        no_valid_cells_no_border=no_valid_cells_previous;
    end
    
%Delete isolated cells (valid cell surrounded by no valid cells) from
%valid_cells_no_border_no_neighbors (same process than before)
    
    flag=0;
    valid_cells_previous=valid_cells_no_border_no_neighbors;
    no_valid_cells_previous=no_valid_cells_2;
    for j=1:length(valid_cells_previous) 
        neighbors_j=neighbors_real{valid_cells_previous(j)};
        no_coincidence=[];
        for i=1:length(neighbors_j)
            no_coincidence(i)=isempty(find(neighbors_j(i)==no_valid_cells_previous, 1));
        end
        if sum(no_coincidence)==0
            valid_cells_no_border_no_neighbors=valid_cells_previous(valid_cells_previous~=valid_cells_previous(j));
            no_valid_cells=sort([valid_cells_previous(j);no_valid_cells_previous]);
            flag=1;
        end
    end
    valid_cells_no_border_no_neighbors=valid_cells_no_border_no_neighbors';
    if (flag==0)
        valid_cells_no_border_no_neighbors=valid_cells_previous;
        no_valid_cells=no_valid_cells_previous;
    end
    


%%%%%%%%%%%      CALCULATE VALID CELLS WITH A GRAPH      %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mark as valid cells, all them that make a small independent graph
graph_bin=zeros(length(num_cells),length(num_cells));
for i=1:size(valid_cells_no_border,1)
    n_cell=valid_cells_no_border(i);
    neighbor_network=neighbors_real{n_cell};

    %graph development
    for ii=1:(length(neighbor_network))
        graph_bin(n_cell,neighbor_network(ii))=1;
        graph_bin(neighbor_network(ii),n_cell)=1;
    end
end
%Delete connections with no valid cells
for i=1:(length(num_cells))
    if isempty(find(i==valid_cells_no_border, 1))==1
        graph_bin(i,:)=0;
        graph_bin(:,i)=0;
    end    
end

%Count the number of graphs
for i=1:(length(num_cells))
   v=find(graph_bin(i,:)==1);
   for j=1:length(v)
       graph_bin(i,:)=(graph_bin(v(j),:)|graph_bin(i,:));
       graph_bin(v(j),:)=(graph_bin(v(j),:)|graph_bin(i,:));
   end
end
Sum=sum(graph_bin,1);
dat=unique(Sum);
grafo_may=dat(1,end);
valid_cells_no_border=find(Sum==grafo_may); %comment this line if error in data extraction due to the number of valid cells (due to presence of isolated valid cells)
valid_cells_no_border=valid_cells_no_border';
no_valid_cells_no_border=(1:(length(num_celulas)));%comment this line if error in data extraction due to the number of valid cells (due to presence of isolated valid cells)
no_valid_cells_no_border(valid_cells_no_border)=[];%comment this line if error in data extraction due to the number of valid cells (due to presence of isolated valid cells)
no_valid_cells_no_border=no_valid_cells_no_border';
    
%% FEATURES EXTRACTION

Img=imread('Original_image.jpg'); %load rgb original image (r:slow cells ; g: collagen ; b:dapi)
%We want adquire information about channel R intensity in cellular regions
%to classify cells as slow or fast cell.
R=Img(:,:,1);
R=imadjust(R); %enhance

%%%% Average areas, centroids,num_cells, slow and fast cells areas %%%%

s=regionprops(cells,'Area');
cell_area=cat(1,s.Area);

Area_valid_cells=cell_area(valid_cells_no_border);
Mean_Area=mean(Area_valid_cells);
Std_Area=std(Area_valid_cells);

%Get slow intensity in for each cellular region (cells slow = slow cells)
Mean_R = regionprops(cellular_mask, R, 'MeanIntensity'); 
Mean_R = cat(1, Mean_R.MeanIntensity);

%A cell will be a slow cell if the 15% of slow intensity of cellular region
%is >=100, or his average is > 50.(this values could be modified depends of
%the intensity image), else cell will be a fast cell.
for i=1:max(max(cellular_mask))
    Intensities=R(cellular_mask==i);
    n_int=length(find(Intensities>=100));
    porc_intense_slow(i)=n_int/length(Intensities);
end
Mean_R_valid_cells=Mean_R(valid_cells_no_border);
porc_intense_slow_valid_cells=porc_intense_slow(valid_cells_no_border);

% Number of slow and fast cells
n_slow_cells=length(unique([find(Mean_R_valid_cells>50)' find(porc_intense_slow_valid_cells>=0.15)]));
n_fast_cells=length(Mean_R_valid_cells)-n_slow_cells;

Average_slow_cells=n_slow_cells/(n_fast_cells+n_slow_cells);

%Slow and fast area
Mean_slow_cells_area=mean(cell_area(unique([find(Mean_R_valid_cells>50)' find(porc_intense_slow_valid_cells>=0.15)])));
Std_slow_cells_area=std(cell_area(Mean_R_valid_cells>50));
Mean_fast_cells_area=mean(cell_area(Mean_R_valid_cells<=50));
Std_fast_cells_area=std(cell_area(Mean_R_valid_cells<=50));

%elongation (major and minor axis cells)
major_axis = regionprops(cellular_mask,'MajorAxisLength');
major_axis = cat(1, major_axis.MajorAxisLength);
Mean_major_axis = mean(major_axis(valid_cells_no_border));

minor_axis = regionprops(cellular_mask,'MinorAxisLength');
minor_axis = cat(1, minor_axis.MinorAxisLength);
Mean_minor_axis = mean(minor_axis(valid_cells_no_border));

Relation_axis=major_axis./minor_axis;
Mean_relation_axis(Relation_axis(valid_cells_no_border));
Std_relation_axis=std(Relation_axis(valid_cells_no_border));

%convex hull
Pix_convex_region=regionprops(cellular_mask,'Solidity');
Pix_convex_region = cat(1, Pix_convex_region.Solidity);
Mean_Pix_convex_region=mean(Pix_convex_region(valid_cells_no_border));
Std_Pix_convex_region=std(Pix_convex_region(valid_cells_no_border));

% angle in relation with X axis
angles = regionprops(cellular_mask,'Orientation');
angles = cat(1, angles.Orientation);
Mean_angles = mean(angles(valid_cells_no_border));
Std_angles = std(angles(valid_cells_no_border));

%centers
center = regionprops(cellular_mask,'Centroid'); 
centers = cat(1, center.Centroid);


%Area water
L_contour_water = bwlabel(contour_watershed,8);
area_water = regionprops(L_contour_water, 'Area');
area_water = cat(1, area_water.Area);
[value,ix]=max(area_water);
contour_watershed(L_contour_water==ix)=0;
contour_watershed(contour_watershed~=0)=1;

L_contour_water = contour_watershed.*cells;
area_water = regionprops(L_contour_water, 'Area');
area_water = cat(1, area_water.Area);
cellular_mask=im2bw(cellular_mask,0.5);



%Go around the cells to clean the remaining connections
for i=1:max(max(L_contour_water))
    v=find(cellular_mask(L_contour_water==i)==0);
    prob_deleting=(length(v)/area_water(i))*100; %% If a cell in cellular_mask is really small in comparison with a cell of cells
    if prob_deleting>80
        contour_watershed(L_contour_water==i)=0;
    end
end
L_contour_water = bwlabel(contour_watershed,8);
L_contour_water(L_contour_water>0)=1;
L_contour_water=L_contour_water.*cells;

s  = regionprops(L_contour_water, 'Area');
Area_water = cat(1, s.Area);

%Relation between cells areas (cells) and real cells areas (cellular_mask)
Relation_valid_areas=Area_water(valid_cells_no_border)./cell_area(valid_cells_no_border);

Mean_relation_areas=mean(Relation_valid_areas);
Std_relation_areas=std(Relation_valid_areas);

%Mean of neighbors
n_valid_neighbors=n_neighbors(valid_cells_no_border);
Mean_neighbors=mean(n_valid_neighbors);
Std_neighbors=std(n_valid_neighbors);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Represent cellular network with a slow background%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux=cellular_mask;
aux(aux~=0)=1;

aux2=aux+contour;


pi_contour=find(aux2==0);
pi_contour2=find(contour==1);
R_c=ones(H,W);
G_c=ones(H,W);
B_c=ones(H,W);

R_c(pi_contour)=1;
G_c(pi_contour)=0;
B_c(pi_contour)=0;

R_c(pi_contour2)=0.8;
G_c(pi_contour2)=0.8;
B_c(pi_contour2)=0.8;

with_slow=cat(3,R_c,G_c,B_c);
count=0;
ind_centers=2;
imshow(with_slow)

% Building network, we are going to add links between valid cells centroids
% overwriting a watershed image with slow background.

i_slow=0;
i_fast=0;
relation_area_neighbors=zeros(1,1);
relation_major_axis_neighbors=zeros(1,1);
relation_minor_axis_neighbors=zeros(1,1);
relation_relation_axis_neighbors=zeros(1,1);
relation_Pix_region_convex_neighbors=zeros(1,1);
relation_Angles_neighbors=zeros(1,1);
relation_Relation_areas_neighbors=zeros(1,1);
graph_binary=zeros(size(cell_area,1),size(cell_area,1));
graph=zeros(size(cell_area,1),size(cell_area,1));
slow_cells=zeros(1);
fast_cells=zeros(1);
n_neighbors_of_slow=zeros(1);
n_slow_neighbors_of_slow=zeros(1);
n_fast_neighbors_of_slow=zeros(1);
n_neighbors_of_fast=zeros(1);
n_slow_neighbors_of_fast=zeros(1);
n_fast_neighbors_of_fast=zeros(1);

for i=1:length(valid_cells_no_border_no_neighbors)
    n_cell=valid_cells_no_border_no_neighbors(i);
    neighbors_network=neighbors_real{n_cell};
    
    % calculate: X of each cell/[ mean(X neighbors cells) + std(X neighbors cells)]
    relation_area_neighbors(i)=cell_area(n_cell)/(mean(cell_area(neighbors_network))+std(cell_area(neighbors_network)));
    relation_major_axis_neighbors(i)=major_axis(n_cell)/(mean(major_axis(neighbors_network))+std(major_axis(neighbors_network)));
    relation_minor_axis_neighbors(i)=minor_axis(n_cell)/(mean(minor_axis(neighbors_network))+std(minor_axis(neighbors_network)));
    relation_relation_axis_neighbors(i)=(major_axis(n_cell)./minor_axis(n_cell))/(mean((major_axis(neighbors_network))./minor_axis(neighbors_network)));
    relation_Pix_region_convex_neighbors(i)=Pix_convex_region(n_cell)/(mean(Pix_convex_region(neighbors_network))+std(Pix_convex_region(neighbors_network)));
    relation_Angles_neighbors(i)=angles(n_cell)/(mean(angles(neighbors_network))+std(angles(neighbors_network)));
    relation_Relation_areas_neighbors(i)=Relation_valid_areas(i)/(mean(Relation_valid_areas(i))+std(Relation_valid_areas(i)));
    
end


for i=1:length(valid_cells_no_border)
    n_cell=valid_cells_no_border(i);
    neighbors_network=neighbors_real{n_cell};

        
    % Neighbors of slow and fast cells
    
    %Slow cells if conditions are true, and add a slow marker to Figure
    if Mean_R(n_cell)>100 || porc_intense_slow(n_cell)>=0.15
        hold on,
        plot(round(centers(n_cell,1)),round(centers(n_cell,2)),'.r','MarkerSize',20)
        
        i_slow=i_slow+1;
        slow_cells(i_slow)=n_cell;
        n_neighbors_of_slow(i_slow)=length(neighbors_real{n_cell});  
        if ~isempty(neighbors_real{n_cell})
            n_slow_neighbors_of_slow(i_slow)=length(union(find(Mean_R(neighbors_real{n_cell})>100), find(porc_intense_slow(neighbors_real{n_cell})>=0.15)));
            n_fast_neighbors_of_slow(i_slow)=length(intersect(find(Mean_R(neighbors_real{n_cell})<=100), find(porc_intense_slow(neighbors_real{n_cell})<0.15)));
        end
        
    %Fast cells, add fast marker to Figure.    
    else
        hold on,
        plot(round(centers(n_cell,1)),round(centers(n_cell,2)),'.g','MarkerSize',20)
        
        i_fast=i_fast+1;
        fast_cells(i_fast)=n_cell;
        n_neighbors_of_fast(i_fast)=length(neighbors_real{n_cell});
        if ~isempty(neighbors_real{n_cell})
            n_slow_neighbors_of_fast(i_fast)=length(union(find(Mean_R(neighbors_real{n_cell})>100), find(porc_intense_slow(neighbors_real{n_cell})>=0.15)));
            n_fast_neighbors_of_fast(i_fast)=length(intersect(find(Mean_R(neighbors_real{n_cell})<=100), find(porc_intense_slow(neighbors_real{n_cell})<0.15)));
        end
    end
        
    for ii=1:(length(neighbors_network))
        graph_binary(n_cell,neighbors_network(ii))=1;
        graph_binary(neighbors_network(ii),n_cell)=1;
    end
    
    for ii=1:(length(neighbors_network))
        graph(n_cell,neighbors_network(ii))=sqrt((centers(n_cell,1)-centers(neighbors_network(ii),1))^2+(centers(n_cell,2)-centers(neighbors_network(ii),2))^2);
        graph(neighbors_network(ii),n_cell)=sqrt((centers(n_cell,1)-centers(neighbors_network(ii),1))^2+(centers(n_cell,2)-centers(neighbors_network(ii),2))^2);
    end
    
    %Draw graph connections between neighbors
    for k=1:length(neighbors_network)
         hold on,
        plot([round(centers(n_cell,1)),round(centers(neighbors_network(k),1))], [round(centers(n_cell,2)),round(centers(neighbors_network(k),2))],'Color','fast')
    end
    
    
end

%Save the network figure where you want
stringres=strcat(name,'_network.bmp');
print('-f1','-dbmp',stringres)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Features slow cells - neighbors
if i_slow>0
    Mean_neighbors_of_slow=mean(n_neighbors_of_slow);
    Std_neighbors_of_slow=std(n_neighbors_of_slow);
    Mean_slow_neighbors_of_slow=mean(n_slow_neighbors_of_slow);
    Std_slow_neighbors_of_slow=std(n_slow_neighbors_of_slow);
    Mean_fast_neighbors_of_slow=mean(n_fast_neighbors_of_slow);
    Std_fast_neighbors_of_slow=std(n_fast_neighbors_of_slow);
else
    Mean_neighbors_of_slow=0;
    Std_neighbors_of_slow=0;
    Mean_slow_neighbors_of_slow=0;
    Std_slow_neighbors_of_slow=0;
    Mean_fast_neighbors_of_slow=0;
    Std_fast_neighbors_of_slow=0;
end

%Features fast cells - neighbors
if i_fast>0
    Mean_neighbors_of_fast=mean(n_neighbors_of_fast);
    Std_neighbors_of_fast=std(n_neighbors_of_fast);
    Mean_slow_neighbors_of_fast=mean(n_slow_neighbors_of_fast);
    Std_slow_neighbors_of_fast=std(n_slow_neighbors_of_fast);
    Mean_fast_neighbors_of_fast=mean(n_fast_neighbors_of_fast);
    Std_fast_neighbors_of_fast=std(n_fast_neighbors_of_fast);
else
    Mean_neighbors_of_fast=0;
    Std_neighbors_of_fast=0;
    Mean_slow_neighbors_of_fast=0;
    Std_slow_neighbors_of_fast=0;
    Mean_fast_neighbors_of_fast=0;
    Std_fast_neighbors_of_fast=0;
end

% Network features
graph_rect=graph(valid_cells_no_border,valid_cells_no_border);


%CHANGE THE DIRECTORY LOOKING FOR 'CODIGO BCT' to apply Prueba_brain
%function
cd DIRECTORY


[n_connection_each_nodes, Sum_weights_each_node, Correlation_between_grades_nodes, Densidad_conexiones, Coef_cluster, T, optime_structure,max_modularity, Matrix_shortest_distances_all_nodes,lambda,efficiency,ecc,radius,diameter, BC]=Prueba_brain(graph,graph_binary,graph_rect,valid_cells_no_border);
cd ..

Sum_weights=Sum_weights_each_node(valid_cells_no_border);

Mean_Relation_areas_neighborhood=mean(relation_area_neighbors);
Std_Relation_areas_neighborhood=std(relation_area_neighbors);

Mean_relation_major_axis_neighbors=mean(relation_major_axis_neighbors);
Std_relation_major_axis_neighbors=std(relation_major_axis_neighbors);
Mean_relation_minor_axis_neighbors=mean(relation_minor_axis_neighbors);
Std_relation_minor_axis_neighbors=std(relation_minor_axis_neighbors);
Mean_relation_relation_axis_neighbors=mean(relation_relation_axis_neighbors);
Std_relation_relation_axis_neighbors=std(relation_relation_axis_neighbors);
Mean_relation_Pix_convex_region_neighbors=mean(relation_Pix_region_convex_neighbors);
Std_relation_Pix_convex_region_neighbors=std(relation_Pix_region_convex_neighbors);
Mean_relation_Angles_neighbors=mean(relation_Angles_neighbors);
Std_relation_Angles_neighbors=std(relation_Angles_neighbors);
Mean_relation_relation_areas_neighbors=mean(relation_Relation_areas_neighbors);
Std_relation_relation_areas_neighbors=std(relation_Relation_areas_neighbors);


Mean_sum_weights=mean(Sum_weights);
Desv_sum_weights=std(Sum_weights);

if i_fast==0
    
    Mean_weights_fast_cells=0;
    Desv_weights_fast_cells=0;
    Mean_Coef_cluster_fast=0;
    Desv_Coef_cluster_fast=0;
    Mean_excentricity_fast=0;
    Desv_excentricity_fast=0;
    Mean_BC_fast=0;
    Desv_BC_fast=0;
    
else
    
    Mean_weights_fast_cells=mean(Sum_weights_each_node((fast_cells)));
    Desv_weights_fast_cells=std(Sum_weights_each_node((fast_cells)));
    Mean_Coef_cluster_fast=mean(Coef_cluster(fast_cells));
    Desv_Coef_cluster_fast=std(Coef_cluster(fast_cells));
    Mean_excentricity_fast=mean(ecc(fast_cells));
    Desv_excentricity_fast=std(ecc(fast_cells));
    Mean_BC_fast=mean(BC(fast_cells));
    Desv_BC_fast=std(BC(fast_cells));
    
    
end



if i_slow==0
    Mean_weights_slow_cells=0;
    Desv_weights_slow_cells=0;
    Mean_Coef_cluster_slow=0;
    Desv_Coef_cluster_slow=0;
    Mean_excentricity_slow=0;
    Desv_excentricity_slow=0;
    Mean_BC_slow=0;
    Desv_BC_slow=0;
else
    
    Mean_weights_slow_cells=mean(Sum_weights_each_node((slow_cells)));
    Desv_weights_slow_cells=std(Sum_weights_each_node((slow_cells)));
    Mean_Coef_cluster_slow=mean(Coef_cluster(slow_cells));
    Desv_Coef_cluster_slow=std(Coef_cluster(slow_cells));
    Mean_excentricity_slow=mean(ecc(slow_cells));
    Desv_excentricity_slow=std(ecc(slow_cells));
    Mean_BC_slow=mean(BC(slow_cells));
    Desv_BC_slow=std(BC(slow_cells));
    
end



Mean_Coef_cluster=mean(Coef_cluster);
Desv_Coef_cluster=std(Coef_cluster);
Mean_excentricity=mean(ecc(valid_cells_no_border));
Desv_excentricity=std(ecc(valid_cells_no_border));
Mean_BC=mean(BC(valid_cells_no_border));
Desv_BC=std(BC(valid_cells_no_border));



M=triu(Matrix_shortest_distances_all_nodes(valid_cells_no_border,valid_cells_no_border));
Mean_dist=mean(M(M~=0));
Desv_dist=std(M(M~=0));


if i_fast~=0 
    M1=triu(Matrix_shortest_distances_all_nodes(fast_cells,fast_cells));
    Mean_dist_fast_fast=mean(M1(M1~=0));
    Desv_dist_fast_fast=std(M1(M1~=0));
else
    Mean_dist_fast_fast=0;
    Desv_dist_fast_fast=0;
end

if i_slow~=0
    M3=triu(Matrix_shortest_distances_all_nodes(slow_cells,slow_cells));
    Mean_dist_slow_slow=mean(M3(M3~=0));
    Desv_dist_slow_slow=std(M3(M3~=0));
else
    Mean_dist_slow_slow=0;
    Desv_dist_slow_slow=0;
end

if i_slow~=0 && i_fast~=0
    M2=triu(Matrix_shortest_distances_all_nodes(fast_cells,slow_cells));
    Mean_dist_fast_slow=mean(M2(M2~=0));
    Desv_dist_fast_slow=std(M2(M2~=0));
    
    M4=triu(Matrix_shortest_distances_all_nodes(slow_cells,fast_cells));
    Mean_dist_slow_fast=mean(M4(M4~=0));
    Desv_dist_slow_fast=std(M4(M4~=0));
    
else 
    Mean_dist_fast_slow=0;
    Desv_dist_fast_slow=0;
    
    Mean_dist_slow_fast=0;
    Desv_dist_slow_fast=0;
end





%CHANGE DIRECTORY TO GERGANA_BOUNOVA TO EXECUTE 'Prueba_Bounova' FUNCTION
%TO GET A FEW OF NETWORK FEATURES

cd Gergana_Bounova
[prs, a, s]=Prueba_Bounova(graph,graph_binary,graph_rect,valid_cells_no_border);
cd ..\..
Coef_pearson=prs;
algebraic_connectivity=a;
Metrics_s=s;

Assortattivity=Correlation_between_grades_nodes;
Densidad_conexiones;
Transitivity=T;
max_modularity;


%% SAVING

%%Save 69 cc (There are more feature extracted, you could save what you want)

stringres=['Results_69_cc.mat'];
save (stringres, 'Mean_Area', 'Std_Area', 'Mean_slow_cells_area','Std_slow_cells_area', 'Mean_fast_cells_area', 'Std_fast_cells_area', 'Mean_major_axis', 'Mean_minor_axis','Mean_relation_axis','Std_relation_axis','Mean_Pix_convex_region','Std_Pix_convex_region','Mean_relation_areas','Std_relation_areas','Mean_neighbors','Std_neighbors','Std_neighbors_of_slow','Std_neighbors_of_fast','Mean_slow_neighbors_of_slow','Mean_fast_neighbors_slow','Mean_slow_neighbors_of_fast','Mean_fast_neighbors_of_fast','Mean_Relation_areas_neighborhood','Std_Relation_areas_neighborhood', 'Mean_relation_major_axis_neighbors','Std_relation_major_axis_neighbors','Mean_relation_minor_axis_neighbors','Std_relation_minor_axis_neighbors','Mean_relation_relation_axis_neighbors','Std_relation_relation_axis_neighbors','Mean_relation_Pix_convex_region_neighbors','Std_relation_Pix_convex_region_neighbors','Mean_relation_Angles_neighbors','Std_relation_Angles_neighbors','Mean_sum_weights','Desv_sum_weights','Mean_weights_fast_cells','Desv_weights_fast_cells','Mean_weights_slow_cells','Desv_weights_slow_cells','Mean_Coef_cluster','Desv_Coef_cluster','Mean_Coef_cluster_fast','Desv_Coef_cluster_fast','Mean_Coef_cluster_slow','Desv_Coef_cluster_slow','Mean_excentricity','Desv_excentricity','Mean_excentricity_fast','Desv_excentricity_fast','Mean_excentricity_slow','Desv_excentricity_slow','Mean_BC','Desv_BC','Mean_BC_fast','Desv_BC_fast','Mean_BC_slow','Desv_BC_slow','Mean_dist','Desv_dist','Mean_dist_fast_fast','Desv_dist_fast_fast','Mean_dist_fast_slow','Desv_dist_fast_slow','Mean_dist_slow_slow','Desv_dist_slow_slow','Mean_dist_slow_fast','Desv_dist_slow_fast','Average_slow_cells')



end