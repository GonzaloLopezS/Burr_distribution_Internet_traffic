%%% Funcion para el calculo de Mixturas Burr

function calculoKmeans(idx,c,kbim,serie_temporal,cq,cdf_burr, param4,dns,cadena_carac,titulo, k,edges)

    clusters = zeros(length(serie_temporal),kbim);
    cluster_merge = sort(serie_temporal);
    
    cq_sums = clusters;
    PD_c = zeros(3,kbim);
    pdf_c = clusters;
    cdf_c = clusters;
    
    for i=1:kbim
        kluster = sort(serie_temporal(idx==i));
        knt = 0;
        for j=1:i-1
            clusters(knt+1:knt+length(nonzeros(clusters(:,j))),i) = zeros(length(nonzeros(clusters(:,j))),1);
            knt = knt + length(nonzeros(clusters(:,j)));
        end
        clusters(knt+1:knt+length(nonzeros(kluster)),i) = kluster;
        [cq_aux, xq_aux] = ecdf(clusters(:,i));
        
        cq_aux = [zeros(knt-(i-1),1); cq_aux];
        xq_aux = [zeros(knt-(i-1),1); xq_aux];
        
        ldif = length(cq_sums(:,i)) - length(cq_aux);
        if ldif > 0
            cq_sums(:,i) = [cq_aux; ones((length(cq_sums(:,i)) - length(cq_aux)),1)];
            xq_sums(:,i) = [xq_aux; zeros((length(cq_sums(:,i)) - length(xq_aux)),1)];
        elseif ldif < 0
            cq_sums(:,i) = [cq_aux(1+abs(ldif):end)];
            xq_sums(:,i) = [xq_aux(1+abs(ldif):end)];
        end
        
        cq_sums(:,i) = cq_sums(:,i)*length(nonzeros(clusters(:,i)))/(sum(length(clusters)));
        PD_c_aux = fitdist(nonzeros(clusters(:,i)),'burr');
        
        PD_c(1,i) = PD_c_aux.alpha; %Parametro alpha
        PD_c(2,i) = PD_c_aux.c; %Parametro c
        PD_c(3,i) = PD_c_aux.k; %Parametro k
        
        pdf_c(:,i) = (pdf('burr',cluster_merge, PD_c(1,i), PD_c(2,i), PD_c(3,i)))*length(nonzeros(clusters(:,i)))/(sum(length(clusters)));
        cdf_c(:,i) = (cdf('burr',cluster_merge, PD_c(1,i), PD_c(2,i), PD_c(3,i)))*length(nonzeros(clusters(:,i)))/(sum(length(clusters)));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pdf_suma = sum(pdf_c');
    pdf_suma = pdf_suma';
    pdf_suma_dr = [0; pdf_suma];

    cdf_suma = sum(cdf_c');
    cdf_suma = cdf_suma';

    x_empirica = [2*cluster_merge(1)-cluster_merge(2);cluster_merge];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculo funcion de error
    N = 7000
    diferencia_burr_max = zeros(N,1);
    
    y_bootstrap_c = zeros(900,kbim); %length: 2500, 1000, 900...
    x_bootstrap_burr_c = y_bootstrap_c;
    cq_sim_burr_c = zeros(length(y_bootstrap_c)+1,kbim);
    knt = 0;

    %Para simular la bimodal, rehago cada burr por su parte con el tamanyo
    %de sus clusters asociados y despues se aplica la mixtura entre ambas
    for i=1:N
        for iaux=1:kbim
            y_bootstrap_c(:,iaux) = rand(size(cluster_merge));
            for j=1:length(y_bootstrap_c(:,iaux))
                y_bootstrap_c(j,iaux) = min(y_bootstrap_c(j,iaux),0.9999);
            end
            
            x_bootstrap_burr_c(:,iaux) = icdf('burr',y_bootstrap_c(:,iaux),PD_c(1,iaux),PD_c(2,iaux),PD_c(3,iaux)); 
            x_bootstrap_burr_c(:,iaux) = sort(x_bootstrap_burr_c(:,iaux));
            
            if max(x_bootstrap_burr_c(:,iaux)) == Inf
                knt = knt+1;
                continue;
            end        
            cq_sim_aux = ecdf(x_bootstrap_burr_c(:,iaux));

            try
                cq_sim_burr_c(:,iaux) = cq_sim_aux*length(nonzeros(clusters(:,iaux)))/(sum(length(clusters)));
            catch
                cq_sim_burr_c(:,iaux) =  cq_sim_burr_c(:,iaux);
            end  
        end

        [cluster_merge_burr, cq_merge_burr] = sumaClusters_n(x_bootstrap_burr_c,cq_sim_burr_c,kbim);

        %Resolver aqui
        valor_maximo = maxima_diferencia(cluster_merge,cq,cluster_merge_burr,cq_merge_burr); %Hay que hacer bootstrap 'a pedales' de la se�al de las medidas??
        diferencia_burr_max(i,1) = valor_maximo;

    end

    % Estadisticos
    [freqs,xi] = ksdensity(diferencia_burr_max);
    [~,idx] = max(freqs);
    moda_multimodal = xi(idx)
    media_multimodal = mean(diferencia_burr_max)
    varianza_multimodal = var(diferencia_burr_max)
    desv_multimodal = std(diferencia_burr_max)
    asimetria_multimodal = skewness(diferencia_burr_max) %Analisis de asimetria
    mediana_multimodal = median(diferencia_burr_max)
    curtosis_multimodal = kurtosis(diferencia_burr_max) %Analisis de la curtosis
    percentil_95_multimodal = prctile(diferencia_burr_max,95)
    
    [ecdf_error,x_cdf_error] = ecdf(diferencia_burr_max);
    
    figure;
    subplot(3,2,1)
    plot(serie_temporal)
    title(dns)
    xlabel('ICMP')
    ylabel('Tiempo de respuesta [ms]');
    
    subplot(3,2,3)
    histogram(cluster_merge,100,'Normalization','pdf');
    hold on
    plot([edges(1);cluster_merge],pdf_suma_dr,'r','linewidth',1.5)
    hold off
    title({strcat('PDF Burr Multimodal: ','  ',dns), strcat('clusters= ',num2str(kbim), ',  ',param4)})
    xlabel('retardo [ms]')
    ylabel('pdf')
    
    subplot(3,2,5)
    plot(cluster_merge,cdf_burr,'g','linewidth',1.3)
    hold on
    plot(cluster_merge,cdf_suma,'r','linewidth',1.3)
    plot(cluster_merge,cq,'linewidth',1);
    hold off
    title('CDF Burr y Mixturas')
    legend('CDF Unimodal','CDF Multimodal','ECDF datos reales')
    xlabel('retardo [ms]')
    ylabel('probabilidad acumulada')
    
    %Error estimado
    subplot(3,2,2)
    plot(diferencia_burr_max,'.')
    title('Diferencias maximas. Ajuste Burr');
    xlabel('Iteraciones')
    ylabel('error estimado')

    ax(1) = subplot(3,2,4)
    histogram(nonzeros(diferencia_burr_max),60,'Normalization','probability')
    ylims = get(gca,'YLim');
    title(strcat('Error multimodal '+ dns));
    display(strcat('Parametros para ', num2str(kbim),' clusters'))
    hold on, line([moda_multimodal moda_multimodal], ylims,'LineStyle','--','LineWidth',1.3,'Color','red'), hold off
    hold on, line([media_multimodal media_multimodal], ylims,'LineStyle','-','LineWidth',1.5,'Color','red'), hold off
    hold on, line([mediana_multimodal mediana_multimodal], ylims, 'LineStyle','-','LineWidth',1.5,'Color', 'black'), hold off
    hold on, line([percentil_95_multimodal percentil_95_multimodal], ylims, 'LineStyle','--','LineWidth',1.3,'Color', 'black'), hold off
    title(strcat('Distribución de Error (Multimodal). '+ dns));
    legend('error estimado',['moda = ' num2str(moda_multimodal)],['media = ' num2str(media_multimodal)],['mediana = ' num2str(mediana_multimodal)],['p95 = ' num2str(percentil_95_multimodal)])
    xlabel('Error estimado')
    ylabel('frecuencia relativa')
    
    try
        [ecdf_error,x_cdf_error] = ecdf(nonzeros(diferencia_burr_max));
    catch
        warning('error linea 510')
    end
    
    ax(2) = subplot(3,2,6)
    plot(x_cdf_error,ecdf_error);
    hold on, line([percentil_95_multimodal percentil_95_multimodal], [0 .95],'LineStyle','--','LineWidth',1.3,'Color','black'), hold off
    hold on, line([0 percentil_95_multimodal], [.95 .95],'LineStyle','--','LineWidth',1.3,'Color','black'), hold off
    title('ECDF Error')
    xlabel('error')
    ylabel('error acumulado')
    set(gcf, 'Position', get(0, 'Screensize'));
    linkaxes(ax,'x','y')
end
