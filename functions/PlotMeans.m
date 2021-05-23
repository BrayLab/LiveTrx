function[Fig] = PlotMeans(MaxF,TimeScale,Selected,Properties,Bits,ColorArg,Fig,OnlyMeans,Legend, XLim, YLimits)
    set(0, 'DefaulttextInterpreter', 'none')
    PropertiesSelected = Properties(Selected,:);
    ShortMid = PropertiesSelected.Type=='ShortMidline';
    LongMid = PropertiesSelected.Type=='LongMidline';
    Midline = PropertiesSelected.Type~='EarlyOnly' & PropertiesSelected.Type~='LateTrack';
    
    MaxFSelected = MaxF(:,Selected);
      
    figure(Fig);
    if OnlyMeans == 0;
    try;subplot(411); hold on
        
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('F (AU)')
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
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('F (AU)')
        xlim([0,XLim(2)]);title('All midline'); 
        legend('show');
        t2=plot(TimeScale,MaxFSelected(:,Midline),'- .','Color',[0.7,0.7,0.7,0.2],...
            'MarkerEdgeColor',ColorArg(2,:),'HandleVisibility','off');
        m2=plot(TimeScale,nanmean(MaxFSelected(:,Midline),2),'.','MarkerSize',12,...
            'MarkerEdgeColor',ColorArg(6,:),'DisplayName',Legend);
        uistack(m2, 'down');
        uistack(t2, 'bottom');
        legend('boxoff');
        end
    try;subplot(413); hold on
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('F (AU)')
        xlim([0,XLim(2)]); title('Short Midline');
        legend('show');
        t3=plot(TimeScale,MaxFSelected(:,ShortMid),'- .','Color',[0.7,0.7,0.7,0.2],...
            'MarkerEdgeColor',ColorArg(3,:),'HandleVisibility','off');
        m3=plot(TimeScale,nanmean(MaxFSelected(:,ShortMid),2),'.','MarkerSize',12,...
            'MarkerEdgeColor',ColorArg(7,:),'DisplayName',Legend);
        uistack(m3, 'down');
        uistack(t3, 'bottom');
        legend('boxoff');
        end
    try; subplot(414); hold on
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('F (AU)')
        xlim([0,XLim(2)]); title('Long Midline');
        legend('show');
        t4=plot(TimeScale,MaxFSelected(:,LongMid),'- .','Color',[0.7,0.7,0.7,0.2],...
            'MarkerEdgeColor',ColorArg(4,:),'HandleVisibility','off');
        m4=plot(TimeScale,nanmean(MaxFSelected(:,LongMid),2),'.','MarkerSize',12,...
            'MarkerEdgeColor',ColorArg(8,:),'DisplayName',Legend);
        uistack(m4, 'down');
        uistack(t4, 'bottom');
        legend('boxoff');
    end
          
    end

    if OnlyMeans == 1
    %try; 
        subplot(411); hold on
        ax = gca;
        ax.YAxisLocation = 'origin';
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('MeanF (AU)')
        xlim(XLim); yticks([0:500:YLimits(2)]);
        set(get(gca, 'XLabel'), 'Position', [XLim(2)*0.9,-0.1*YLimits(2),0]);
        set(get(gca, 'YLabel'), 'Position', [-XLim(2)*0.05,YLimits(2)*1.1,0]);
        %title('All selected');
        legend('show');
        %plot(TimeScale,nanmean(MaxFSelected,2),'- .','MarkerSize',12,...
        %   'Color',ColorArg(5,:),'MarkerEdgeColor',ColorArg(5,:),'DisplayName',Legend);
        errorbar(TimeScale,nanmean(MaxFSelected,2),nanstd(MaxFSelected,1,2)./sqrt(sum(~isnan(MaxFSelected),2)),...
            '- .','MarkerSize',8,'Color',[ColorArg(6,:),0.5],'DisplayName',Legend,'CapSize',0,'LineWidth',0.25);
    
        legend boxoff
        legend('Location','NorthWest')
       
    %     end
    try;
        subplot(412); hold on
        ax = gca;
        ax.YAxisLocation = 'origin';
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('MeanF (AU)')
        xlim([-0.01,XLim(2)]); yticks([0:500:YLimits(2)]);
        set(get(gca, 'XLabel'), 'Position', [XLim(2)*0.9,-0.1*YLimits(2),0]);
        set(get(gca, 'YLabel'), 'Position', [-XLim(2)*0.05,YLimits(2)*1.1,0]);
        %title('All midline');
        legend('show');
        %plot(TimeScale,nanmean(MaxFSelected(:,Midline),2),'- .','MarkerSize',12,...
        %    'Color',ColorArg(6,:),'MarkerEdgeColor',ColorArg(6,:),'DisplayName',Legend); 
        errorbar(TimeScale,nanmean(MaxFSelected(:,Midline),2),nanstd(MaxFSelected(:,Midline),1,2)./sqrt(sum(~isnan(MaxFSelected(:,Midline)),2)),...
            '- .','MarkerSize',8,'Color',[ColorArg(6,:),0.5],'DisplayName',Legend,'CapSize',0,'LineWidth',0.25);
        legend boxoff
        legend('Location','NorthWest')

         end
    try;subplot(413); hold on
        ax = gca;
        ax.YAxisLocation = 'origin';
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('MeanF (AU)')
        xlim([-0.01,XLim(2)]); yticks([0:500:YLimits(2)]);
        set(get(gca, 'XLabel'), 'Position', [XLim(2)*0.9,-0.1*YLimits(2),0]);
        set(get(gca, 'YLabel'), 'Position', [-XLim(2)*0.05,YLimits(2)*1.1,0]);
        %title('Short Midline');
        legend('show');
%         plot(TimeScale,nanmean(MaxFSelected(:,ShortMid),2),'- .','MarkerSize',12,...
%             'Color',ColorArg(7,:),'MarkerEdgeColor',ColorArg(7,:),'DisplayName',Legend); 
        errorbar(TimeScale,nanmean(MaxFSelected(:,Midline),2),nanstd(MaxFSelected(:,Midline),1,2)./sqrt(sum(~isnan(MaxFSelected(:,Midline)),2)),...
            'Color',[ColorArg(6,:),0.5],'DisplayName',Legend,'CapSize',0,'LineWidth',0.5);
        
        legend boxoff
        legend('Location','NorthWest')
         end
    try; subplot(414); hold on
        ax = gca;
        ax.YAxisLocation = 'origin';
        ylim(YLimits); xlabel('time into nc14 (min)'); ylabel('MeanF (AU)')
        xlim([-0.01,XLim(2)]); yticks([0:500:YLimits(2)]);
        set(get(gca, 'XLabel'), 'Position', [XLim(2)*0.9,-0.1*YLimits(2),0]);
        set(get(gca, 'YLabel'), 'Position', [-XLim(2)*0.05,YLimits(2)*1.1,0]);
        %title('Long Midline');
        legend('show');
%         plot(TimeScale,nanmean(MaxFSelected(:,LongMid),2),'- .','MarkerSize',12,...
%             'Color',ColorArg(8,:),'MarkerEdgeColor',ColorArg(8,:),'DisplayName',Legend);
         plot(TimeScale,nanmean(MaxFSelected(:,Midline),2),'- .','MarkerSize',12,...
             'Color',ColorArg(6,:),'MarkerEdgeColor',ColorArg(6,:),'DisplayName',Legend); 
        legend boxoff
        legend('Location','NorthWest')
          end
    end
end