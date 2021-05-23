function[] = BoxplotsActivity(SelectedN,SelectedNames,Which,Info,PathPlots,Exps,Selections, Palette,XLim,varargin)
            ThresholdSplit = 6;
            try
                if ~isempty(varargin{1})
                    ThresholdSplit = varargin{1};
                end
            end
%             SelectedN = {'m5m8peveIII wRi','m5m8peveIII ArmRi','m5m8peveIII aCatRi'}
%             SelectedNames = {'wRi','ArmRi','aCatRi'}
%             Selection = 'MSE'
%             PaletteDefault = [130,130,130;80,80,80; %grey1, grey2. 1,2
%                     142,183,36;105,139,34;62,81,16; %green1,green2,green3. 3:5
%                     61,131,183;54,100,139;33,63,86; %blue1,blue2,blue3. 6:8
%                     250,174,64;226,140,26; %yellow1, yellow2. 9:10
%                     167,51,170;104,15,107; %pink1, pink2. 11:12
%                     78,69,139; 54, 44,96; %purple1, purple2. 13:14
%                     205,55,0]./255; %red. 15
%             Palette = PaletteDefault([1,15,13],:)

            set(0,'defaultAxesColorOrder',Palette)
            set(0,'defaultLegendAutoUpdate','on')

            ToSave = [PathPlots,char(join(SelectedN,'_vs_'))];

            Fig2 = figure('PaperSize',[32 12],'PaperUnits','inches','resize','on', 'visible','on');

            TotalCellsAll = cell(1,length(SelectedN));
            CellsAAll = cell(1,length(SelectedN));
            CellsBAll = cell(1,length(SelectedN));
            TmRNAAll = cell(1,length(SelectedN));
 %           Table = cell2table(cell(0,7), 'VariableNames', {'File','Nickname','Repeat','GastT1', 'GasT2', 'GastT3','TranscT2'}); 

            for i = 1:length(SelectedN)
                        if length(Selections) == length(SelectedN)
                            Selection = Selections{i}
                        else
                            Selection = Selections{1};
                        end
                    %try
                        Index = find(cellfun(@(x) strcmp(x,SelectedN{i}),table2array(Exps(:,Which)))==1)';
                        if strcmp(SelectedN{i},'all'); Index = 1:height(Exps); end
                        for x = 1:length(Index)
                            %try
                            Experiment = [Exps.Nickname{Index(x)},' ',num2str(Exps.Rep(Index(x)))]
                            PathToSave = [Info.Path{Index(x)},Info.File{Index(x)},...
                            Info.Name{Index(x)},Info.File{Index(x)}]; 
                            Parameters = Info(Index(x),:);
                            Table2Vars(Parameters);
                            load([PathToSave,'_Data.mat']);
                            Struct2Vars(Data)
                            Norm = ((MaxF-Baseline').*Baseline(1)./Baseline').*(2.^(12-Bits));
                            Norm(Norm==0) = NaN;
                            SplitEarlyF = max([SplitEarly*60./TimeRes+nc14-Delay,1]);
                            %XLim = [0,50];
                            XLimF = XLim(2)*60./TimeRes+nc14-Delay;
                            ElT = 5000/2000*60;
                            Properties.TotalmRNA = [nansum(Norm(SplitEarlyF:min(size(Norm,1),XLimF),:),1)]'* TimeRes ./ 2 ./ ElT;
                            %Properties.TotalmRNA = [nansum(Norm(SplitEarlyF:size(Norm,1),:),1)]'* TimeRes ./ 2 ./ ElT;;
                            Selected = Properties.Region==Selection & Properties.Onset < XLim(2);
                            TmRNA = Properties(Selected,:).TotalmRNA
                            TmRNA(TmRNA<0) = 0;
                            plot(sort(TmRNA)); hold on
                            TotalCells = sum(Selected);
                            CellsA = sum(TmRNA >= (ThresholdSplit*1000));
                            CellsB = sum(TmRNA < (ThresholdSplit*1000));
                            TotalCellsAll{i} = [TotalCellsAll{i},TotalCells];
                            CellsAAll{i} = [CellsAAll{i},CellsA];
                            CellsBAll{i} = [CellsBAll{i},CellsB];
                            TmRNAAll{i} = [TmRNAAll{i},TmRNA'];
                            %end
                            %NewTable = table(string(File),string(Nickname),Rep,GastT(1),GastT(2),GastT(3),PeakT(2),'VariableNames', {'File','Nickname','Repeat','GastT1', 'GasT2', 'GastT3','TranscT2'});
                            %Table = [Table;NewTable];
                       PlotTotal = 1;     
                       if PlotTotal 
                            try
                                load([PathToSave, '_Stats_tracked.mat']);
                            catch
                                load([PathToSave, '_Stats.mat']);
                                Stats_tracked = Stats_GFP;
                                clear Stats_GFP
                            end
                            F = max(round(35*60/TimeRes)+nc14-Delay,1);
                            if F >= Frames ; F = Frames-1; end
                            %Width = size(Im,1); Height = size(Im,2);
                                try
                                    Merged_meanF_maxGFP = Read3d([PathToSave, '_replacedF_maxGFP_MCP.tiff']);
                                catch
                                    try
                                        Merged_meanF_maxGFP = Read3d([PathToSave, '_replacedF_maxGFP_mean.tiff']);
                                    catch
                                        Merged_meanF_maxGFP = Read3d([PathToSave, '_maxF_maxGFP.tiff']);
                                    end
                                end
                            Im = Merged_meanF_maxGFP(:,:,F)./255;
                            Width = size(Im,1); 
                            Height = size(Im,2);
                            MaxTracked = zeros(Width,Height,Frames);
                             SingleFrame = Stats_tracked{F};
                             Tracked3D = zeros(Width,Height,Slices);
                             for n = 1:size(SingleFrame,1)
                                 Tracked3D(SingleFrame.PixelIdxList{n}) =  SingleFrame.Label(n);
                             end
                             MaxTracked= max(Tracked3D,[],3);
                            SelectedLabels = Properties.Label;
                            TotalProd = zeros(size(MaxTracked));
                            for n = 1:length(SelectedLabels)
                                NewLabels = Properties.NewLabel(n);
                                TotalProd(ismember(MaxTracked,NewLabels{:})) =  Properties(n,:).TotalmRNA;
                            end
                            for n = 1:length(BGLabels)
                                TotalProd(ismember(MaxTracked,BGLabels(n))) =  1; %noBG
                            end
                            %
                            maxTotRNA = ThresholdSplit*2*1000; % 3000 for lateral. 10000 for ventral
                            %TotalProd = cat(3,TotalProdEarly,TotalProdMid);
                            TotalProd (TotalProd < 0) = 1;
                            TotalProd (TotalProd > maxTotRNA) = maxTotRNA;
                            TotalProdRBG = zeros(size(TotalProd,1), size(TotalProd,2),3,size(TotalProd,3));
                            CMAP = colormap(parula(maxTotRNA+1));
                            for f = 1:size(TotalProd,3)
                                TotalProdRBG(:,:,:,f) = label2rgb(round(TotalProd(:,:,f)),CMAP,'k');
                            end
                            WriteRGB(double(TotalProdRBG), PathPlots,[Info.File{Index(x)},'_TotmRNA_',num2str(SplitEarly),'-',num2str(XLim(2)),'_MAX',num2str(ThresholdSplit*2),'.tiff'],'none')
                       end    
                       end 


            end

             figure(Fig2)

            TotalCellsAll = Cell2Mat(TotalCellsAll);
            CellsAAll = Cell2Mat(CellsAAll);
            CellsBAll = Cell2Mat(CellsBAll);
            TmRNAAll = Cell2Mat(TmRNAAll);
            Jitter = 0.5; %Jitter./2 cant be > BarW 
            BarW = 0.35; FaceAlpha = 0.3; DotSize = 15; LineWidth = 2; FontSizeTitle = 16; FontSize = 14;
            set(gca,'FontSize',FontSize)

            subplot(151); hold on     
            plotBoxplot([TotalCellsAll],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,50])
            ylabel('# active MSE cells')

            subplot(153); hold on     
            plotBoxplot([CellsAAll],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,30])
            ylabel(['# high transcribing MSE cells (> ',num2str(ThresholdSplit),' AU)'])

            subplot(154); hold on     
            plotBoxplot([CellsBAll],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,30])
            ylabel(['# low transcribing MSE cells (< ',num2str(ThresholdSplit),' AU)'])

            subplot(1,5,2); hold on
            
            
                plotViolin(TmRNAAll./1000,SelectedNames,SelectedNames,0.25,BarW,'total mRNA (AU)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,0,0)
            
            subplot(1,5,5); hold on
            for i = 1:size(TmRNAAll,2) 
            	Hist = histogram(TmRNAAll(:,i)./1000,'Normalization','Probability','BinWidth',2,'FaceColor',Palette(i,:),'FaceAlpha',0.3,'LineStyle','none')
                plot(Hist.BinEdges(1:end-1)+Hist.BinWidth/2,Hist.Values,'Color',Palette(i,:),'LineWidth',2)
            end
            title(['Accumulated signal ',num2str(SplitEarly),' - ',num2str(XLim(2))])
            ylabel('Proportion of cells')
            xlabel('Total mRNA production (AU)')

            Selection = strjoin(Selections);

            print(Fig2,[ToSave,'_',Selection,'_TotmRNA_',num2str(SplitEarly),'-',num2str(XLim(2)),'_split',num2str(ThresholdSplit),'.pdf'],'-fillpage', '-dpdf');

            close all

end