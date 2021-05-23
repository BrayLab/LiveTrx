function[FigH] = PlotHeatmaps(FigH,NormMerged,Selected,PropertiesMerged,TimeScaleMerged,SelectedN,TimeRes,XLim,i,Max,XLabel,varargin)
        try
            SizeEach = varargin{1};
        catch
            SizeEach = [];
        end
        %PropertiesMerged() = [];
        NormMerged(NormMerged < 0) = 1;
        %NormMerged(isnan(NormMerged)) = 0;
        NormMerged = NormMerged(:,Selected);
        PropertiesMerged = PropertiesMerged(Selected,:);
        [~,nc14Merged] = min(abs(round(TimeScaleMerged,4) - XLim(1)));
        [~,End] = min(abs(round(TimeScaleMerged,4) - XLim(2)));
        %nc14Merged = max(1,find());
        %End = find(round(TimeScaleMerged,4) == XLim(2));
        %if isempty(End); End = size(NormMerged,1); end
        %if isempty(DelayMerged); DelayMerged = 0; end
        DelayMerged = round((TimeScaleMerged(nc14Merged) - XLim(1)).*60./TimeRes);
        Bits = 12;
        CMAP = parula(2^Bits+1);
        %CMAP = gray(2^Bits+1);
        CMAP(1,:) = [0,0,0];
        Cases = {'Onset','Total','SelectedTotal'};
        Cases = {'Onset','Total','None'};
        figure(FigH);
        SelectedTotal = 0;
        for n = 1:length(Cases)
            Sort = Cases{n};
            subplot(3,length(SelectedN),(n-1)*length(SelectedN)+i)
            switch Sort
                case 'None'
                    Sorted = NormMerged;
                    YLabel = 'grouped by embryo';
                case 'AP'
                    [~,Order] = sort(PropertiesMerged.NormAP);
                    Sorted = NormMerged(:,Order);
                    YLabel = ['sorted by AP position'];
                case 'Total'
                    [~,Order] = sort(PropertiesMerged.TotalmRNA);
                    Sorted = NormMerged(:,Order);
                    %YLabel = [num2str(size(NormMerged,2)),' cells sorted by total mRNA production'];
                    YLabel = ['sorted by total mRNA production'];
                case 'Onset'
                    [~,Order] = sort(PropertiesMerged.Onset);
                    Sorted = NormMerged(:,Order);
                    YLabel = ['sorted by onset time'];
                case 'DV'
                    [~,Order] = sort(PropertiesMerged.DMSE,'descend');
                    Sorted = NormMerged(:,Order);
                    YLabel = ['sorted by DV position'];
                case 'SelectedTotal'
                    SelectedTotal = 1;
                    [~,Order] = sort(PropertiesMerged.TotalmRNA);
                    Sorted = NormMerged(:,Order);
                    N = size(Sorted,2);
                    ToPlotEach = Sorted(:,[round(N*0.66),round(N*0.95)]);
                    
            end
            if SelectedTotal == 0
                Heatmap = zeros(size(Sorted,2),round((XLim(2)-XLim(1))*60./TimeRes)+1);
                Heatmap(:,(1:(End-nc14Merged+1))+DelayMerged) = Sorted(nc14Merged:End,:)';
                imagesc(XLim,[1,size(Sorted,2)],Heatmap,[0,Max]);
                %xlim([0,(XLim(2)-XLim(1))]*60./TimeRes)
                %XTickLabels = str2double(get(gca,'XTickLabels'));
                %xticks([0:20:(XLim(2)-XLim(1))]*60/TimeRes)
                %xticklabels(cellstr(num2str([XLim(1):20:XLim(2)]')))
                yticklabels('')
                xlabel(XLabel)
                ylabel(YLabel)
                set(get(gca, 'YLabel'), 'Units', 'Normalized','Position', [-0.02,0.5,0],'Rotation',90);
                colormap(CMAP)
                title([SelectedN{i},' (',num2str(size(NormMerged,2)),' cells)'])
                 if strcmp(Sort,'None')
                    hold on
                    BarsEachRep = cumsum(SizeEach);
                    set(gca, 'ColorOrder', [170,170,170;110,110,110]./255);
                    plot(repmat(mean(XLim),2,length(SizeEach)),[0,BarsEachRep(1:end-1);BarsEachRep],'-','LineWidth',1)
                end
            else
                plot(TimeScaleMerged,ToPlotEach); hold on;
                %plot(TimeScale,OnOffSelected(:,(n-1)*30+i).*(2^12-1),'.');
                ylim([-100, 2^12-1]);
                xlim(XLim);
                %title(join(['#',num2str(PropertiesSelected.Label((n-1)*30+i)),' ',PropertiesSelected.Type((n-1)*30+i)]));
                hold off
            end
        end
end