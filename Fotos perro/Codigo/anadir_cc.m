function [Mejores_i,Mejores_i_des,Proy, eigenvectorsF]=anadir_cc(Mejores,Mejores_des,vector_todas_caracteristicas,expansion,n_imagenes_tipo1,n_imagenes_tipo2)
cuenta=0;
Niteracion=1;
for paso=1:size(Mejores,1)
    cc=Mejores(paso,2:size(Mejores,2));
    vector_caracteristicas=[];
    for j=1:size(Mejores,2)-1
        vector_caracteristicas=[vector_caracteristicas,vector_todas_caracteristicas(:,cc(1,j))];
        
    end
    for caract=1:size(vector_todas_caracteristicas,2)
        if isempty(find(caract==cc))==1
            p1=sum(vector_todas_caracteristicas(:,caract));
            if p1==0
                Ratio_pca(1,Niteracion)=0;
                Ratio_pca(2,Niteracion)=caract;
                eigenvectors{1,Niteracion} = 0;
                Niteracion=Niteracion+1;
            else
                
                %Voy incluyendo caracteristica
                vector_caracteristicas_defi=[vector_caracteristicas,vector_todas_caracteristicas(:,caract)];
                %normalizamos los vectores
                for car=1:size(vector_caracteristicas_defi,2)
                    vector_caracteristicas_defi(:,car)=vector_caracteristicas_defi(:,car)-min(vector_caracteristicas_defi(:,car));
                    vector_caracteristicas_defi(:,car)=vector_caracteristicas_defi(:,car)/max(vector_caracteristicas_defi(:,car));
                end
                vector_caracteristicas_defi(isnan(vector_caracteristicas_defi))=0;% PARA HACER 0 TODOS LOS NAN
                %% PCA
                X=vector_caracteristicas_defi';
                %Calculo la media
                Media=mean(X,2);
                
                % Le resto la media a cada imagen
                for i = 1:size(X,2)
                    X(:,i) = X(:,i) - Media;
                end
                
                L = X'*X;
                
                % Cálculo de los autovalores/autovectores
                [Vectors,Values] = eig(L);
                [d_ordenados,ind]=sort(diag(Values),'descend');   % se ordenan los autovalores
                V=Vectors(:,ind);
                
                
                
                %Convierto los autovalores de X'X en los autovectores de X*X'
                Vectors = X*V;
                
                % Normalizo los vectores, para conseguir longitud unitaria
                for i=1:size(X,2)
                    Vectors(:,i) = Vectors(:,i)/norm(Vectors(:,i));
                end
                
                
                V=Vectors(:,1:2);
                
                W{1,Niteracion}=V'*X;  %Proyecciones
                
                %%%% Obtencion de numeros a partir de graficas metodo3 (LUCIANO)
                label=[ones(1, n_imagenes_tipo1), 2*ones(1,n_imagenes_tipo2)];
                [T, sintraluc, sinterluc, Sintra, Sinter] = valid_sumsqures(W{1,Niteracion}',label,2);
                C=sinterluc/sintraluc;
                Ratio_pca(1,Niteracion)=abs(trace(C));
                Ratio_pca(2,Niteracion)=caract;
                eigenvectors{1,Niteracion} = V;

                
                Niteracion=Niteracion+1;
            end
        end
    end
    auxiliar=Ratio_pca(1,:);
    for i=1:expansion
        [Mejores_i(cuenta+i,1) num]=max(auxiliar);
        for j=2:size(Mejores,2)
            Mejores_i(cuenta+i,j)=Mejores(paso,j);
            Mejores_i_des(cuenta+i,j)=Mejores_des(paso,j);
        end
        Mejores_i(cuenta+i,size(Mejores,2)+1)=Ratio_pca(2,num);
        Mejores_i_des(cuenta+i,size(Mejores,2)+1)=Ratio_pca(2,num);
        Proy{cuenta+i,1}=W{1,num};
        eigenvectorsF{cuenta+i,1} = eigenvectors{1, num};
        auxiliar(1,num)=0;
    end
    cuenta=cuenta+expansion;
    clear W Ratio_pca
end

