function clusterBimodal(serie_temporal,kbim,clas,cq,cdf_burr,dns,cadena_carac,titulo,k,edges)
%UNTITLED Summary of this function goes here
%  Analisis y caracterizacion de la serie temporal para casos bimodales.
%  1) calculo de los cluster mediante kmeans y calculo de la funcion
%  de error.

%Tipos de clusterizacion:
% - Hierarchical Clustering -> Creacion de un arbol de clustering
% - K-Means -> Particion del conjunto de datos en k-clusters de acuerdo a
%   la distancia respecto de un centroide de cada cluster.

% - Modelos de mixturas gaussianas -> Modela clusters como mixturas de
%   componentes normales multivariantes.

% - Mapas autoorganizados (Self-organizing maps) -> Usan redes neuronales
%   que aprenden la topologia y distribucion de los datos.

% - Modelos de Markov ocultos -> Ejemplo de procesos estocasticos. Procesos
%   que generan secuencias aleatorias de resultados o estados segun ciertas
%   probabilidades. Procesos sin memoria: su estado siguiente depende
%   unicamente de su estado actual.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% K-Means: -https://es.mathworks.com/help/stats/kmeans.html
% 3er argumento (opcional): Name -> 'Distance ('sqeuclidean'
% predeterminado; centroide = media)', 'euclidean', 'cosine (centroide = media de los puntos tras normalizar a longitud euclidiana de la unidad)'
% ('off' determinado), 'correlation (centroide = media de componentes despues de centrar y normalizar esos puntos a media cero y stdev 1)'

% Otros: 'Replicates', 'Options', statset('UseParallel',1)10 -computacion
% paralela mediante la replicacion de clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - K-medoids -> https://es.mathworks.com/help/stats/kmedoids.html


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vecinos mas cercanos -> https://es.mathworks.com/help/stats/nearest-neighbors-1.html?s_tid=CRUX_lftnav

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clustering jerarquico -> https://es.mathworks.com/help/stats/hierarchical-clustering.html

% La agrupaci�n jer�rquica agrupa los datos en una variedad de escalas mediante la creaci�n de un �rbol de cl�steres o 
% Dendrograma El �rbol no es un solo conjunto de cl�steres, sino m�s bien una jerarqu�a multinivel, donde los cl�steres
% de un nivel se unen como cl�steres en el siguiente nivel. Esto le permite decidir el nivel o la escala de agrupaci�n 
% en cl�steres que es m�s adecuado para la aplicaci�n.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    param4 = 'cityblock'
    [idx1,c1] = kmeans(serie_temporal, kbim, 'Distance', param4);
%     calculoKmeans(idx1,c1,kbim,serie_temporal,cq,cdf_burr, param4,dns,cadena_carac,titulo,k,edges)
    try
        calculoKmeans(idx1,c1,kbim,serie_temporal,cq,cdf_burr, param4,dns,cadena_carac,titulo,k,edges)
     catch
        warning('No funciona con cityblock')
    end

    param3 = 'cosine'
    [idx2,c2] = kmeans(serie_temporal, kbim, 'Distance', param3);
    try
        calculoKmeans(idx2,c2,kbim,serie_temporal,cq,cdf_burr, param3,dns,cadena_carac,titulo,k,edges)
     catch
        warning('No funciona con cosine')
     end
  
    param5 = 'sqeuclidean'
    [idx4,c4] = kmeans(serie_temporal, kbim, 'Distance','sqeuclidean');
    try
        calculoKmeans(idx4,c4,kbim,serie_temporal,cq,cdf_burr, param5,dns,cadena_carac,titulo,k,edges)
     catch
        warning('No funciona con sqeuclidean')
    end
     
    
end
