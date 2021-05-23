function[Fig] = PlotMeans(MaxF,MeanNorm, OnOff,TimeScale,Selected,Properties,Bits,ColorArg,Fig,OnlyMeans,Legend, XLim, YLimits,XLabel)
    set(0, 'DefaulttextInterpreter', 'none')
    PropertiesSelected = Properties(Selected,:);
    %ShortMid = PropertiesSelected.Type=='ShortMidline';
    %LongMid = PropertiesSelected.Type=='LongMidline';
    %ShortMid = strcmp(PropertiesSelected.Type,'ShortMidline');
    %LongMid = strcmp(PropertiesSelected.Type,'LongMidline');
    %Midline = PropertiesSelected.Type~='EarlyOnly' & PropertiesSelected.Type~='LateTrack';
    %Midline = ~strcmp(PropertiesSelected.Type,'EarlyOnly') & ~strcmp(PropertiesSelected.Type,'LateTrack');
    %Midline = ~strcmp(PropertiesSelected.Type,'fkhgfk'); %###############
    MaxFSelected = MaxF(:,Selected);
    OnOffSelected = OnOff(:,Selected);
      

    figure(Fig);
    if OnlyMeans == 0;
    try;subplot(411); hold on
        
        ylim(YLimits); xlabel(XLabel); ylabel('F (AU)')
        xlim([XLim]);title('All selected');
        legend('show');
        t1 = plot(TimeScale,MaxFSelected,'- .','Color',[0.7,0.7,0.7,0.2],...
            'MarkerEdgeColor',ColorArg(1,:),'HandleVisibility','off');
        m1 = plot(TimeScale,nanmean(MaxFSelected,2),'.','MarkerSize',12,...
            'MarkerEdgeColor',ColorArg(5,:),'DisplayName',Legend);
        uistack(m1, 'down');
        uistack(t1, 'bottom');
        legend('boxoff');
         end
    try;subplot(412); hold on
        ylim(YLimits); xlabel(XLabel); ylabel('F (AU)')
        xlim([0,XLim(2)]);title('All midline'); 
        legend('show');
        t2=plot(TimeScale,MaxFSelected,'- .','Color',[0.7,0.7,0.7,0.2],...
            'MarkerEdgeColor',ColorArg(2,:),'HandleVisibility','off');
        m2=plot(TimeScale,nanmean(MaxFSelected,2),'.','MarkerSize',12,...
            'MarkerEdgeColor',ColorArg(6,:),'DisplayName',Legend);
        uistack(m2, 'down');
        uistack(t2, 'bottom');
        legend('boxoff');
        end
    try;subplot(413); hold on
        ylim(YLimits); xlabel(XLabel); ylabel('F (AU)')
        xlim([0,XLim(2)]); title('Short Midline');
        legend('show');
        t3=plot(TimeScale,MaxFSelected,'- .','Color',[0.7,0.7,0.7,0.2],...
            'MarkerEdgeColor',ColorArg(3,:),'HandleVisibility','off');
        m3=plot(TimeScale,nanmean(MaxFSelected,2),'.','MarkerSize',12,...
            'MarkerEdgeColor',ColorArg(7,:),'DisplayName',Legend);
        uistack(m3, 'down');
        uistack(t3, 'bottom');
        legend('boxoff');
        end
    try; subplot(414); hold on
        ylim(YLimits); xlabel(XLabel); ylabel('F (AU)')
        xlim([0,XLim(2)]); title('Long Midline');
        legend('show');
        t4=plot(TimeScale,MaxFSelected,'- .','Color',[0.7,0.7,0.7,0.2],...
            'MarkerEdgeColor',ColorArg(4,:),'HandleVisibility','off');
        m4=plot(TimeScale,nanmean(MaxFSelected,2),'.','MarkerSize',12,...
            'MarkerEdgeColor',ColorArg(8,:),'DisplayName',Legend);
        uistack(m4, 'down');
        uistack(t4, 'bottom');
        legend('boxoff');
    end
          
    end

    if OnlyMeans == 1
        
        subplot(411); hold on
        legend('show');
        Mean = nanmean(MaxFSelected,2);
        SEM = nanstd(MaxFSelected,1,2)./sqrt(sum(~isnan(MaxFSelected),2));
        TimeScalePoly = TimeScale;
        MeanPoly = Mean;
        SEMPoly = SEM;
        TimeScalePoly(isnan(Mean)) = [];
        MeanPoly(isnan(Mean)) = [];
        SEMPoly(isnan(Mean)) = [];
        
        %errorbar(TimeScale,Mean,SEM,'- .','MarkerSize',8,'Color',[ColorArg(6,:),0.5],'DisplayName',[Legend,' n = ',num2str(sum(Midline)),' cells, ',num2str(size(MeanNorm,2)),' embryos'],'CapSize',0,'LineWidth',0.25);
        Poly = polyshape([TimeScalePoly,flip(TimeScalePoly)]',[MeanPoly-SEMPoly;flip(MeanPoly+SEMPoly)]);
        plot(Poly,'LineStyle','none','FaceColor',[ColorArg(6,:),0.5],'HandleVisibility','off');hold on
        plot(TimeScale,Mean,'-','MarkerSize',8,'Color',[ColorArg(6,:),1],'DisplayName',[Legend,' n = ',num2str(size(MaxFSelected,2)),' cells, ',num2str(size(MaxFSelected,2)),' embryos'],'LineWidth',0.5);
        legend boxoff
        legend('Location','NorthWest')
        ax = gca;
        ax.YAxisLocation = 'origin';
        ax.XAxisLocation = 'origin';
        ylim(YLimits); xlabel(XLabel); ylabel('Mean Fluorescence (AU)')
        %xlim(XLim); %yticks([0:500:YLimits(2)]);
        %set(get(gca, 'XLabel'), 'Position', [XLim(2),-0.1*YLimits(2),0]);
        set(get(gca, 'XLabel'), 'HorizontalAlignment', 'right');
        set(get(gca, 'XLabel'), 'VerticalAlignment', 'top');
        set(get(gca, 'YLabel'), 'Position', [-XLim(2)*0.05,YLimits(2)*1.02,0]);
        set(get(gca, 'YLabel'), 'HorizontalAlignment', 'left');
        set(get(gca, 'YLabel'), 'VerticalAlignment', 'bottom');
        set(get(gca,'YLabel'),'Rotation',0)

   
        subplot(412); hold on
        title('Mean + SEM merged movies');
        legend('show');
        %plot(TimeScale,nanmean(MaxFSelected(:,Midline),2),'- .','MarkerSize',12,...
        %    'Color',ColorArg(6,:),'MarkerEdgeColor',ColorArg(6,:),'DisplayName',Legend); 
        Mean = nanmean(MaxFSelected,2);
        SEM = nanstd(MaxFSelected,1,2)./sqrt(sum(~isnan(MaxFSelected),2));
        TimeScalePoly = TimeScale; 
        MeanPoly = Mean; SEMPoly = SEM; 
        TimeScalePoly(isnan(Mean)) = [];
        MeanPoly(isnan(Mean)) = [];
        SEMPoly(isnan(Mean)) = [];
        %errorbar(TimeScale,Mean,SEM,'- .','MarkerSize',8,'Color',[ColorArg(6,:),0.5],'DisplayName',[Legend,' n = ',num2str(sum(Midline)),' cells, ',num2str(size(MeanNorm,2)),' embryos'],'CapSize',0,'LineWidth',0.25);
        Poly = polyshape([TimeScalePoly,flip(TimeScalePoly)]',[MeanPoly-SEMPoly;flip(MeanPoly+SEMPoly)]);
        plot(Poly,'LineStyle','none','FaceColor',[ColorArg(6,:),0.5],'HandleVisibility','off');hold on
        plot(TimeScale,Mean,'-','MarkerSize',8,'Color',[ColorArg(6,:),1],'DisplayName',[Legend,' n = ',num2str(size(MaxFSelected,2)),' cells, ',num2str(size(MeanNorm,2)),' embryos'],'LineWidth',0.5);
        legend boxoff
        legend('Location','NorthWest')
        ax = gca;
        %ax.YAxisLocation = 'origin';
        ax.XAxisLocation = 'origin';
        ylim(YLimits); xlabel(XLabel); ylabel('Mean Fluorescence (AU)')
        xlim(XLim); %yticks([0:500:YLimits(2)]);
        set(get(gca, 'XLabel'), 'Position', [XLim(2),-0.1*YLimits(2),0]);
        set(get(gca, 'XLabel'), 'HorizontalAlignment', 'right');
        set(get(gca, 'XLabel'), 'VerticalAlignment', 'top');
        set(get(gca, 'YLabel'), 'Position', [XLim(1),YLimits(2)*1.02,0]);
        set(get(gca, 'YLabel'), 'HorizontalAlignment', 'left');
        set(get(gca, 'YLabel'), 'VerticalAlignment', 'bottom');
        set(get(gca,'YLabel'),'Rotation',0)

   
        subplot(413); hold on
         title('Mean + SEM each movie');
        legend('show');
        %plot(TimeScale,nanmean(MaxFSelected(:,Midline),2),'- .','MarkerSize',12,...
        %    'Color',ColorArg(6,:),'MarkerEdgeColor',ColorArg(6,:),'DisplayName',Legend); 
        Mean = nanmean(MeanNorm,2);
        SEM = nanstd(MeanNorm,1,2)./sqrt(sum(~isnan(MeanNorm),2));
        %errorbar(TimeScale,Mean,SEM,'- .','MarkerSize',8,'Color',[ColorArg(6,:),0.5],'DisplayName',[Legend,' n = ',num2str(sum(Midline)),' cells, ',num2str(size(MeanNorm,2)),' embryos'],'CapSize',0,'LineWidth',0.25);
         TimeScalePoly = TimeScale; 
        MeanPoly = Mean; SEMPoly = SEM; 
        TimeScalePoly(isnan(Mean)) = [];
        MeanPoly(isnan(Mean)) = [];
        SEMPoly(isnan(Mean)) = [];
        Poly = polyshape([TimeScalePoly,flip(TimeScalePoly)]',[MeanPoly-SEMPoly;flip(MeanPoly+SEMPoly)]);
        plot(Poly,'LineStyle','none','FaceColor',[ColorArg(6,:),0.25],'HandleVisibility','off');hold on
        plot(TimeScale,MeanNorm,':','MarkerSize',8,'Color',[ColorArg(6,:),0.5],'HandleVisibility','off','LineWidth',0.25);
        plot(TimeScale,Mean,'-','MarkerSize',8,'Color',[ColorArg(6,:),1],'DisplayName',[Legend,' n = ',num2str(size(MeanNorm,2)),' embryos'],'LineWidth',1);
        legend boxoff
        legend('Location','NorthWest')
        ax = gca;
        ax.YAxisLocation = 'origin';
        ax.XAxisLocation = 'origin';
        ylim(YLimits); xlabel(XLabel); ylabel('Mean Fluorescence (AU)')
        xlim(XLim); %yticks([0:500:YLimits(2)]);
        set(get(gca, 'XLabel'), 'Position', [XLim(2),-0.1*YLimits(2),0]);
        set(get(gca, 'XLabel'), 'HorizontalAlignment', 'right');
        set(get(gca, 'XLabel'), 'VerticalAlignment', 'top');
        set(get(gca, 'YLabel'), 'Position', [XLim(1),YLimits(2)*1.02,0]);
        set(get(gca, 'YLabel'), 'HorizontalAlignment', 'left');
        set(get(gca, 'YLabel'), 'VerticalAlignment', 'bottom');
        set(get(gca,'YLabel'),'Rotation',0)
       
        
        subplot(414); hold on
        title('Proportion cells ON');
        legend('show');
%         plot(TimeScale,nanmean(MaxFSelected(:,ShortMid),2),'- .','MarkerSize',12,...
%             'Color',ColorArg(7,:),'MarkerEdgeColor',ColorArg(7,:),'DisplayName',Legend); 
%         plot(TimeScale,[nansum(OnOffSelected(:,Midline),2)./nansum(~isnan(OnOffSelected(:,Midline)),2)],...
%             'Color',[ColorArg(6,:),0.5],'DisplayName',Legend,'LineWidth',0.5);
        yyaxis right
        plot(TimeScale,[nansum(OnOffSelected,2)./nansum(~isnan(OnOffSelected),2)],...
            '--','MarkerSize',8,...
             'Color',ColorArg(6,:),'MarkerEdgeColor',ColorArg(6,:),'LineWidth',1,'DisplayName',[Legend,' n = ',num2str(size(MaxFSelected,2))]); 
        legend boxoff
        legend('Location','NorthWest')
        %ax = gca;
        %ax.YAxisLocation = 'origin';
        ylim([0,1]); xlabel(XLabel); ylabel('Fraction of nuclei ON');
        xlim(XLim); 
        set(get(gca, 'YLabel'), 'Position', [XLim(2),1.02,0]);
        set(get(gca, 'YLabel'), 'HorizontalAlignment', 'left');
        set(get(gca, 'YLabel'), 'VerticalAlignment', 'bottom');
        set(get(gca,'YLabel'),'Rotation',0)
        
        yyaxis left
        plot(TimeScale,nansum(OnOffSelected,2),...
            '-','MarkerSize',8,...
             'Color',ColorArg(6,:),'MarkerEdgeColor',ColorArg(6,:),'LineWidth',1,'DisplayName',[Legend,' n = ',num2str(size(MaxFSelected,2))]); 
        legend boxoff
        legend('Location','NorthWest')
        %ax = gca;
        %ax.YAxisLocation = 'origin';
        xlabel(XLabel); ylabel('Number of nuclei ON');
        xlim(XLim); 
        set(get(gca, 'XLabel'), 'Position', [XLim(2),-0.1,0]);
        set(get(gca, 'XLabel'), 'HorizontalAlignment', 'right');
        set(get(gca, 'XLabel'), 'VerticalAlignment', 'top');
        set(get(gca, 'YLabel'), 'Position', [XLim(1),1.02,0]);
        set(get(gca, 'YLabel'), 'HorizontalAlignment', 'left');
        set(get(gca, 'YLabel'), 'VerticalAlignment', 'bottom');
        set(get(gca,'YLabel'),'Rotation',0)
        
 
    end
end