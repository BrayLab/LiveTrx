function[] = plotBoxplot(All,Nicknames,ExpLabels,Jitter,BarW,Title,FontSize,DotSize,CMAP,FaceAlpha,LineWidth,FontSizeTitle,Ylim)
%set(gcf,'defaultAxesColorOrder',CMAP)

Cat = repmat([1:1:length(Nicknames)],2,1)+[-BarW;+BarW];
RepAll = ones(size(All));

CatX = repmat([1:1:length(Nicknames)],size(All,1),1);
gscatter(CatX(:)+(rand(length(CatX(:)),1)-0.5).*Jitter, All(:), CatX(:),CMAP,'.',DotSize); hold on
patch([Cat;flip(Cat)],[repmat(quantile(All,0.25),2,1);repmat(quantile(All,0.75),2,1)],[1:length(Nicknames)]','EdgeColor','none','FaceAlpha',FaceAlpha); hold on
%plot(Cat,repmat(nanmean(All),2,1),'-','DisplayName','mean','Color',[0.3,0.3,0.3],'LineWidth',LineWidth)
%errorbar([1:1:length(Nicknames)],nanmean(All), nanstd(All), nanstd(All),'Color',[0.3,0.3,0.3,0.8],'LineStyle','none','LineWidth',LineWidth./2)
%plot(Cat,repmat(quantile(NumOn,0.25)*100,2,1),'-','DisplayName','Q1','Color',[0.6980,0.0941,0.1686,0.8],'LineWidth',LineWidth)
%plot(Cat,repmat(quantile(NumOn,0.75)*100,2,1),'-','DisplayName','Q3','Color',[0.6980,0.0941,0.1686,0.8],'LineWidth',LineWidth)
plot(Cat,repmat(nanmedian(All),2,1),'-','DisplayName','mean','LineWidth',LineWidth)
errorbar([1:1:length(Nicknames);nan(1,length(Nicknames))],[nanmean(All);nanmean(All)],...
    [nanstd(All);nanstd(All)], [nanstd(All);nanstd(All)],'LineStyle','-','LineWidth',LineWidth./2)


legend off
xticks([1:1:length(Nicknames)])
xticklabels(ExpLabels)
%set(findobj(gca, 'type', 'line'), 'linew',2)
%set(findall(gca,'type','text'),'fontSize',4,'fontWeight','bold')
set(gca,'FontSize',FontSize)
title(Title,'FontSize',FontSizeTitle); box off; 
    xlim([BarW,length(Nicknames)+1-BarW])
    ylim([Ylim])
colormap(CMAP)


end