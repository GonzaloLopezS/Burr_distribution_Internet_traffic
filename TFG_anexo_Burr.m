clear all;
close all;
clc;

% Gonzalo Lopez Segovia [TFG: Modelado Estadistico de Trafico de Internet mediante Distribuciones Burr]
prompt = 'Elige de donde quieres obtener las capturas: LAN (default), Servidor UAM (1)';
entrada = input(prompt, 's');

if strcmp(entrada,'1')
    display('Entrando en la BD de capturas desde redes UAM');
    cd capturas_tshark_uam
else
    display('Entrando en la BD de capturas desde red domestica');
    cd capturas_tshark
end
%Panel de seleccion de archivo a procesar:


if length(strfind(pwd,[filesep 'capturas_tshark_uam'])) ~= 0
    %Captura desde host conectado a red UAM
    prompt = 'Seleccione captura:\n(1)elpais\n(2)wikipedia\n(3)live.com\n(4)twitter\n(5)amazon.com\n(6)blogspot\n(7)okdiario\n(8)facebook\n(otherwise)youtube\n';
    val = input(prompt);
    nombre = 'lista_ip.txt';
    path_uam = true;
    pwd_uam = pwd;
    switch val
        case 1
            archivo = 'elpais_15032019.txt';
        case 2
            archivo = 'wikipedia_11032019.txt';
        case 3
            archivo = 'livecom_10032019.txt';
        case 4
            archivo = 'twitter_12032019.txt';
        case 5
            archivo = 'amazon_10032019.txt';
        case 6
            archivo = 'blogspot_13032019.txt';
        case 7
            archivo = 'okdiario_11032019.txt';
        case 8
            archivo = 'facebook_10032019.txt';    
        otherwise
            archivo = 'youtube_uam_09032019.txt';
    end
    else
%     %Capturas desde host conectado a red local
    prompt = 'Seleccione captura:\n(1)elpais\n(2)google\n(3)twitter\n(4)facebook\n(5)youtube\n(6)amazon\n(7)blogspot\n(8)live.com\n(otherwise)wikipedia\n';
    val = input(prompt);
    nombre = 'lista_ip.txt';
    switch val
        case 1
            archivo = 'elpais_04012019.txt';
        case 2
            archivo = 'google_01012019.txt';
        case 3
            archivo = 'twitter_03012019.txt';
        case 4
            archivo = 'facebook_02012019.txt';
        case 5
            archivo = 'youtube_02012019.txt';
        case 6
            archivo = 'amazon_02012019.txt';
        case 7
            archivo = 'blogspot_04012019.txt';
        case 8
            archivo = 'livecom_02012019.txt';
        otherwise
            archivo = 'wikipedia_03012019.txt';

    end
end

if length(strfind(pwd,[filesep 'capturas_tshark_uam'])) ~= 0
    fileID = fopen(nombre);
    for i=1:val-1       
        fgetl(fileID);
    end
    line_ex = fgetl(fileID);
    line_split = split(line_ex," ");
    titulo = 'amazon_6grafs_def'; %M
    dns = strcat(" ",line_ex);
    fclose(fileID);
else
    fileID = fopen(nombre);
    for i=1:val-1       
        fgetl(fileID);
    end
    line_ex = fgetl(fileID);
    line_split = split(line_ex," ");
    file_name = line_split{2};
    a = strfind(file_name,'.');
    titulo = file_name(1:a-1)
    ip = line_split{3};
    dns = strcat(file_name," ",ip);
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str_arch = strcat("sed 's/,/./g' ", archivo, " >local.txt");
system(str_arch);
fileID = fopen('local.txt','r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
fclose(fileID); 

muestra = 15; %minutos

% carpeta = 'capturas_comparar_remuestreo';
% if ~exist(strcat(carpeta,[filesep titulo]),'dir')
%     mkdir(strcat(carpeta,[filesep titulo]));
% end
% 
% cambia = true;
% tshark_uam = length(strfind(pwd,[filesep 'capturas_tshark_uam']));
% matrix_csv = zeros(1,7);
% pdf_estable_alm = zeros

h = 6;
pdf_burr_stable = zeros(4*h,2);
cdf_burr_stable = pdf_burr_stable;
tiempos_burr_stable = zeros(4*h,2);

for k=1:(4*h)
    close all;
    
    x = ((k-1)*muestra*60)+1:k*muestra*60;
    t15min = A(x);
    minuto = k*muestra;
    str = [' Hasta ', num2str(k*15), ' min'];
    
    tpo = [num2str((k-1)*15), ' - ', num2str(k*15), 'min'];
    
    %%% HISTOGRAMA
    max1 = max(t15min);
    min1 = min(t15min);

    clas_N = length(t15min);
    clas = 100;
    t15min_sort = sort(t15min);

    [q, edges] = histcounts(t15min_sort,clas);
    xt = edges(1:end-1)+(edges(2)-edges(1))/2;
    
    q = q';
    for w=1:length(q)
       if q(w) == 0
           q(w) = 10.^-9;
       else
           q(w) = q(w);
       end
    end
    
    %%% Ajuste Burr. Representación PDF. Comparativa con alpha-Stable
    try
        ajuste_burr = tic;
        PD = fitdist(t15min_sort,'burr');
        burr = toc(ajuste_burr)
    catch
        warning('Los datos no ajustan bien a una distribución Burr');
        continue;
    end
    pdf_burr_stable(k,1) = burr;
    
    ajuste_estable = tic;
    PD2 = fitdist(t15min_sort,'stable');
    estable = toc(ajuste_estable);
    pdf_burr_stable(k,2) = estable;
    
    alpha_burr = PD.alpha;
    c_burr = PD.c;
    k_burr = PD.k;

    %Representacion histogramas + ajuste de distribuciones
    %Histograma:
    pdf_draw_burr = tic;
    pdf_burr = pdf('burr',t15min_sort,alpha_burr,c_burr,k_burr);
    pdf_burr_dr = [0; pdf_burr];
    pdf_draw_burr = toc(pdf_draw_burr)
    
    figure;
    h_t15min = histogram(t15min_sort,clas,'Normalization','pdf');
    title({'pdf burr';tpo})
    hold on
    plot([edges(1); t15min_sort],pdf_burr_dr,'r','linewidth',1.5);
    hold off
    
    legend('datos reales','pdf burr')
    ylabel('pdf')
    xlabel('retardos [ms]');
    if tshark_uam
        title(strcat('PDF Burr '+ dns));
    else
        title(strcat('PDF Burr '+ dns))
    end

    %Distribucion acumulada (CDF)
    [cq,xq] = ecdf(t15min_sort);
    ldif = length(t15min)  - length(cq);
    dlt = xq(end) -xq(end-1);
    
    vector = zeros(ldif,1);
    for w=1:ldif
        vector(w) = w*dlt + xq(end);
    end
    
    if ldif > 0
        cq = [cq; ones((length(t15min) - length(cq)),1)];
        xq = [xq; vector];
    elseif ldif < 0
        cq = [cq(1+abs(ldif):end)];
        xq = [xq(1+abs(ldif):end)];
    end
    
    %CDF-Burr Distribution Type XII
    cdf_draw_burr = tic;
    QRScdf=cdf('burr', t15min_sort, alpha_burr,c_burr,k_burr);
    cdf_draw_burr = toc(cdf_draw_burr);
    cdf_burr_stable(k,1) = cdf_draw_burr;
    
    figure;
    plot(t15min_sort,QRScdf,'g','linewidth',1.6);
    if tshark_uam
        title({strcat('cdf relativa '+ dns);tpo})
    else
        title({strcat('cdf relativa '+ dns);tpo})
    end
    
    ylabel('probabilidades acumuladas')
    xlabel('retardo')

    hold on
    plot(xq,cq,'r','LineWidth',2)
    hold off
    legend('empirica','CDF Burr')
    
%%%%% BOOTSTRAP:
    y_bootstrap = rand(size(t15min_sort));
    x_bootstrap_burr = zeros(size(t15min_sort));
    
    bootstrap_burr = tic; bootstrap_estable = tic;
    x_bootstrap_burr = icdf('burr',y_bootstrap,alpha_burr,c_burr,k_burr);
    

%%%% CDF Bootstrap:
    x_bootstrap_burr = sort(x_bootstrap_burr);
    [cq_sim_burr, x_burr_sim] = ecdf(x_bootstrap_burr);
        
    figure;
    plot(x_burr_sim,cq_sim_burr,'r','LineWidth',1.5);
    hold on
    plot(xq,cq,'LineWidth',1.5)
    hold off
    legend('remuestreo burr (ajuste)', 'original (medidas reales)')
    if tshark_uam
        title(strcat('comparacion ECDF '+ dns))
    else
        title(strcat('comparacion ECDF '+ dns))
    end
    xlabel('retardos [ms]')
    ylabel('probabilidades acumuladas')
    
    
    N = 7000;
    diferencia_burr_max = zeros(N,1);
    
    knt = 0;
    qnt = 0;
    if cambia
        cd ..
        cambia = false;
    end
    
    for i=1:N %Burr
        y_bootstrap = rand(1000,1);
        for j=1:length(y_bootstrap)
           y_bootstrap(j) = min(y_bootstrap(j),0.9999); 
        end
       
        x_bootstrap_burr = zeros(1000,1);
        x_bootstrap_burr = icdf('burr',y_bootstrap,alpha_burr,c_burr,k_burr);
        x_bootstrap_burr = sort(x_bootstrap_burr);
        
        if max(x_bootstrap_burr) == Inf
            knt = knt+1;
            continue;
        end
        
        [cq_sim,xq_sim] = ecdf(x_bootstrap_burr);
        valor_maximo = maxima_diferencia(xq,cq,xq_sim,cq_sim);
        diferencia_burr_max(i,1) = valor_maximo;
        
        
    end

    if strcmp(entrada,'1')
        file_captura =  'capturas_tshark_uam';
    else
        file_captura =  'capturas_tshark';
    end
    subplot(3,2,1)
    plot(x,t15min);
    if tshark_uam
        title(dns)
    else
        title(strcat(dns, str))
    end
    xlabel('ICMP')
    ylabel('Tiempo de respuesta [ms]');
    
    %fig 2
    subplot(3,2,3)
    h_t15min = histogram(t15min_sort,clas,'Normalization','pdf');
    title({'pdf burr';['datos de ' tpo]})
   
    hold on
    plot([edges(1); t15min_sort],pdf_burr_dr,'r','linewidth',1.5);
    hold off
    
    legend('datos reales','pdf burr')
    ylabel('pdf')
    xlabel('retardos [ms]');
    if tshark_uam
        title(strcat('PDF Burr '+ dns + ' ' + tpo));
    else
        title(strcat('PDF Burr '+ dns + ' ' + tpo));
    end
    
    %fig 3
    subplot(3,2,5)
    plot(t15min_sort,QRScdf,'g','linewidth',1.5)
    if tshark_uam
        title(strcat('cdf relativa '+ dns));
    else
        title(strcat('cdf relativa '+ dns))
    end
    
    ylabel('probabilidades acumuladas')
    xlabel('retardo')
    
    hold on
    plot(xq,cq,'r','LineWidth',1);
    legend('CDF burr','ECDF (datos reales)')
    hold off
    
    %funcion de error
    subplot(3,2,2)
    plot(diferencia_burr_max,'.')
    title(strcat('Diferencias maximas. Ajuste Burr'," ",dns));
    xlabel('Iteraciones')
    ylabel('error estimado')
    
    ax(1) = subplot(3,2,4)
    histogram(nonzeros(diferencia_burr_max),60,'Normalization','probability') %Nonzeros
    ylims = get(gca,'YLim');
    hold on, line([moda_unimodal moda_unimodal], ylims,'LineStyle','--','LineWidth',1.3,'Color','red'), hold off
    hold on, line([media_unimodal media_unimodal], ylims,'LineStyle','-','LineWidth',1.5,'Color','red'), hold off
    hold on, line([mediana_unimodal mediana_unimodal], ylims, 'LineStyle','-','LineWidth',1.5,'Color', 'black'), hold off
    hold on, line([percentil_95_unimodal percentil_95_unimodal], ylims, 'LineStyle','--','LineWidth',1.3,'Color', 'black'), hold off
    title(strcat('Distribución de Error (Unimodal). '+ dns));
    legend('error estimado',['moda = ' num2str(moda_unimodal)],['media = ' num2str(media_unimodal)],['mediana = ' num2str(mediana_unimodal)],['p95 = ' num2str(percentil_95_unimodal)])
    xlabel('Error estimado')
    ylabel('frecuencia relativa')
    
    try
        [ecdf_error,x_cdf_error] = ecdf(nonzeros(diferencia_burr_max));
    catch
        ('error linea 510')
        continue;
    end
    
    ax(2) = subplot(3,2,6)
    plot(x_cdf_error,ecdf_error)
    hold on, line([percentil_95_unimodal percentil_95_unimodal], [0 .95],'LineStyle','--','LineWidth',1.3,'Color','black'), hold off
    hold on, line([0 percentil_95_unimodal], [.95 .95],'LineStyle','--','LineWidth',1.3,'Color','black'), hold off
    title(strcat('ECDF Funcion de error'," ",dns));
    title('ECDF Error')
    xlabel('error')
    ylabel('probabilidad acumulada')
    set(gcf, 'Position', get(0, 'Screensize'));
    linkaxes(ax,'x','y');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Mixturas. Algoritmo K-Means para clusterizar las bimodales
    
    for kbim=2:6
        display(strcat('Calculo para: ', num2str(kbim), ' clusters'))
        try
            clusterBimodal(t15min,kbim,clas,cq,QRScdf,dns,cadena_carac,titulo,k,edges);
        catch
            warning('No ajusta a bimodal')
            continue;
        end
        pause;
    end  

    pause;
    
end



%Calculos de clusterizacion para k =1, 2 y 3. Y ver cual tiene menor error.
%-> Media y desviacion de la funcion de error (OK)

%Coeficiente silhouette -> Identificar no. clusters.
%Comparativa de tiempos entre ajuste burr y alfa-estable. (OK)
