function[] = PlotScatter(GastTime,TranscTime,DisplayName, XLabel, YLabel,Plot,Fit)


TranscTime(isnan(GastTime)) = [];
GastTime(isnan(GastTime)) = [];  
GastTime(isnan(TranscTime)) = [];    
TranscTime(isnan(TranscTime)) = [];

TranscTime(GastTime==0) = [];
GastTime(GastTime==0) = [];  
GastTime(TranscTime==0) = [];    
TranscTime(TranscTime==0) = [];

if Plot
    plot(GastTime,TranscTime,'.','MarkerSize',20,'DisplayName',DisplayName); 
end
if Fit
    hold on
    p = polyfit(GastTime, TranscTime,1);
    yfit = polyval(p, GastTime);
    yresid = TranscTime - yfit;
    rsq = 1 - sum(yresid.^2)/((length(TranscTime)-1)*var(TranscTime));
    plot(GastTime, yfit, '-', 'Color',[0.3,0.3,0.3],'DisplayName',['R^2 = ', num2str(round(rsq,2))]);
    axis square
    xlim([35,75])
    ylim([35,75])
else
    %plot(NaN, NaN, '-','HandleVisibility','off');
end
xlabel(XLabel);
ylabel(YLabel);

Lgn = legend('show','Location','best');
 %Lgn = legend('boxoff')
Lgn.FontSize = 12;
pbaspect([1,1,1])


end