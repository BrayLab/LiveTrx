function[] = plotViolin(All,Nicknames,ExpLabels,BinWidthHist,BarW,Title,FontSize,FaceAlpha,LineWidth,FontSizeTitle,Ymin, Log)
    Cat = repmat([1:1:length(Nicknames)],2,1)+[-BarW;+BarW];

    Probs = {};
    Bins = {};
    for i = 1:length(Nicknames)
        h = histogram(All(~isnan(All(:,i)),i),'BinWidth',BinWidthHist,'Normalization','probability','FaceAlpha',FaceAlpha,'LineStyle','none','Visible','off'); hold on;
        Probs{i} = [0,h.Values];
        Edges = h.BinEdges; Edges(Edges == 0) = 0.1; % so that it doesnt give error when transforming to log scale
        Bins{i} = [Edges];
        
    end
    Probs = Cell2Mat(Probs);
    Bins = Cell2Mat(Bins);
    %Bins(end-1,:) = [];
    %Bins(Probs==0) = NaN;
    %Probs(Probs==0) = NaN;

    CatX = repmat([1:1:length(Nicknames)],length(All),1);
    %h = gscatter(CatX(:)+(rand(length(CatX(:)),1)*0.5-0.5).*Jitter, All(:), RepAll(:),CMAP,'.',DotSize);hold on
    %patch([Cat;flip(Cat)],[repmat(quantile(All,0.25),2,1);repmat(quantile(All,0.75),2,1)],'red','EdgeColor','none','FaceAlpha',FaceAlpha); hold on
    %plot(Cat,repmat(quantile(BLenAll,0.25)*TimeRes/60,2,1),'-','DisplayName','Q1','Color',[0.6980,0.0941,0.1686,0.8],'LineWidth',LineWidth)
    %plot(Cat,repmat(quantile(BLenAll,0.75)*TimeRes/60,2,1),'-','DisplayName','Q3','Color',[0.6980,0.0941,0.1686,0.8],'LineWidth',LineWidth)
    %errorbar([1:1:length(Nicknames)],nanmean(BLenAll*TimeRes/60), nanstd(BLenAll*TimeRes/60)*0, nanstd(BLenAll*TimeRes/60),'Color',[0.6980,0.0941,0.1686,0.8],'LineStyle','none','LineWidth',0.5)
    Probs = medfilt1(Probs,3);
    Probs(isnan(Bins)) = NaN;
        Bins(isnan(Probs)) = NaN;
%         Bins(end+1,:) = 0;
%         Probs(end+1,:) = 0;


    for i = 1:size(Probs,2)
        %SmoothProbs(:,i) = smooth(Probs(:,i));
                %SmoothProbs(:,i) = Probs(:,i);
        ProbsPoly = [Probs(:,i).*BarW./max(Probs(:,i))+i;flipud(-Probs(:,i).*BarW./max(Probs(:,i))+i)];
        ProbsPoly(isnan(ProbsPoly)) = [];
        BinsPoly = [Bins(:,i);flipud(Bins(:,i))]; 
        BinsPoly(isnan(BinsPoly)) = [];
        try
            Poly = polyshape(ProbsPoly,BinsPoly);
                %Poly = polyshape([Probs(:,i).*BarW./max(Probs(:,i))+i,-Probs(:,i).*BarW./max(Probs(:,i))+i],[Bins(:,i),Bins(:,i)])
            plot(Poly,'LineStyle','none');hold on
        end
    end
    
        plot(Cat,repmat(nanmedian(All),2,1),'-','DisplayName','mean','Color',[0.3,0.3,0.3],'LineWidth',LineWidth)

    legend off
    xlim([BarW,length(Nicknames)+1-BarW])
    xticks([1:1:length(Nicknames)])
    xticklabels(ExpLabels)
    set(gca,'FontSize',FontSize)
    if Log
        set(gca, 'YScale', 'log')
    end
    ylim([Ymin, inf])
    title(Title,'FontSize',FontSizeTitle); 
    box off;
    %xlim([0+Jitter,length(Nicknames)+Jitter]); 
    %set(findobj(gca, 'type', 'line'), 'linew',2)
end