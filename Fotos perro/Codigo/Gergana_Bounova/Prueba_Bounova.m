% clear all
% construyo grafo

function  [prs, a, s]=Prueba_Bounova(grafo,grafo_binario,grafo_rect,centros_rect)
% grafo=zeros(6);
% 
% grafo(1,2)=1;
% grafo(2,1)=1;
% grafo(2,3)=1;
% grafo(3,2)=1;
% grafo(2,4)=2;
% grafo(4,2)=2;
% grafo(2,5)=3;
% grafo(5,2)=3;
% grafo(4,6)=2;
% grafo(6,4)=2;
% grafo(5,6)=2;
% grafo(6,5)=2;
% grafo(5,7)=1;
% grafo(7,5)=1;
% grafo(4,5)=1;
% grafo(5,4)=1;
% 
% grafo_binario(1,2)=1;
% grafo_binario(2,1)=1;
% grafo_binario(2,3)=1;
% grafo_binario(3,2)=1;
% grafo_binario(2,4)=1;
% grafo_binario(4,2)=1;
% grafo_binario(2,5)=1;
% grafo_binario(5,2)=1;
% grafo_binario(4,6)=1;
% grafo_binario(6,4)=1;
% grafo_binario(5,6)=1;
% grafo_binario(6,5)=1;
% grafo_binario(5,7)=1;
% grafo_binario(7,5)=1;
% grafo_binario(4,5)=1;
% grafo_binario(5,4)=1;


%% Features
prs = pearson(grafo_rect);
a=algebraic_connectivity(grafo_rect);
s=s_metric(grafo_rect);


%%
% Nodos=getNodes(grafo,'adj')
% Ramas=getEdges(grafo,'adj')
% 
% L = adj2adjL(grafo)
% 
% N_nodos=numnodes(grafo)
% N_ramas=numedges(grafo)
% 
% Densidad=link_density(grafo)
% isweighted(Ramas)
% Grado_medio=average_degree(grafo)

% [GC,gc_nodes]=giant_component(grafo)

% C=closeness(grafo)
% C=closeness(grafo_binario)

%% find the (weighted) eccentricity of all vertexes, radius, diameter, center vertexes and the periphery vertexes;
% ec=vertex_eccentricity(grafo)
% diam = diameter(grafo)
% Rg=graph_radius(grafo)
% diam = diameter(grafo)

%%
%  l = ave_path_length(grafo_binario)
%  l = ave_path_length(grafo)
%  
% 
% 
% ddist=distance_distribution(grafo)