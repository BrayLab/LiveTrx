function[] = PlotBurstingProp(BurstNum,BurstLength,BurstPeriod,TimeRes,FileOut)
    Fig = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','off');
    subplot(231);histogram([BurstPeriod{:}].*TimeRes./60); 
    title('Burst period'); xlabel('Burst period (min)');
    %semilogy(Hist.Values)
    subplot(232);histogram([BurstLength{:}].*TimeRes./60);
    title('Burst length'); xlabel('Burst lenght (min)');
    subplot(233);histogram([BurstNum{:}]); 
    title('Number of bursts'); xlabel('# per cell');
    subplot(234);boxplot([BurstPeriod{:}].*TimeRes./60,'Jitter',0.5,'Symbol','.'); 
    title('Burst period'); xlabel('Burst period (min)');
    %semilogy(Hist.Values)
    subplot(235);boxplot([BurstLength{:}].*TimeRes./60,'Jitter',0.5,'Symbol','.');
    title('Burst length'); xlabel('Burst lenght (min)');
    subplot(236);boxplot([BurstNum{:}],'Jitter',0.5,'Symbol','.'); 
    title('Number of bursts'); xlabel('# per cell');

    print(Fig,FileOut,'-fillpage', '-dpsc');

end