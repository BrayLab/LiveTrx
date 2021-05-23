function[] = CorrelateGastrulation(SelectedN,SelectedNames,Which,Info,PathPlots,Exps, Selections, Palette,varargin)
            Confirm = 0;
            try
                if ~isempty(varargin{1})
                    Confirm = varargin{1};
                end

            end
            %SelectedN = PairstoSelect{n}
            %SelectedNames = Nicknames{n}
            %Palette = PaletteMain(PaletteEach{n},:);
            %Palette2 = kron(Palette,ones(2,1)); %duplicate each row consecutively
            set(0,'defaultAxesColorOrder',Palette)
            set(0,'defaultLegendAutoUpdate','on')

            ToSave = [PathPlots,char(join(SelectedN,'_vs_'))];

            Fig = figure('PaperSize',[45 30],'PaperUnits','inches','resize','on', 'visible','on');
            Fig2 = figure('PaperSize',[40 12],'PaperUnits','inches','resize','on', 'visible','on');
            %Fig2 = figure('PaperSize',[25 12],'PaperUnits','inches','resize','on', 'visible','on');
            TranscTime = cell(1,length(SelectedN));
            GastTime1 = cell(1,length(SelectedN));
            GastTime2 = cell(1,length(SelectedN));
            GastTime3 = cell(1,length(SelectedN));
            Increase = cell(1,length(SelectedN));
            PeakF1 = cell(1,length(SelectedN));
            Table = cell2table(cell(0,7), 'VariableNames', {'File','Nickname','Repeat','GastT1', 'GasT2', 'GastT3','TranscT2'}); 
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
                            PathData = [PathToSave,'_Data.mat'];
                            load(PathData);
                            clear GastValues
                             if contains(Selection,'|')
                                    SelectionToSave = join(Selections,'');SelectionToSave = SelectionToSave{:};
                                 else
                                    SelectionToSave = Selection;
                             end
                            try
                                GastValues = Data.(['GastValues',SelectionToSave]);
                                if Confirm
                                    Data = DrawGastTimes(Data,PathData,Parameters,Selection,Confirm,SplitEarly);
                                    GastValues = Data.(['GastValues',SelectionToSave]);
                                end
                            catch
                                Data = DrawGastTimes(Data,PathData,Parameters,Selection,Confirm,SplitEarly);
                                GastValues = Data.(['GastValues',SelectionToSave]);
                            end
                            PeakT = GastValues(:,1); PeakT(PeakT<0) = NaN;
                            PeakF = GastValues(:,2);
                            GastT = GastValues(:,3); GastT(GastT<0) = NaN;
                            MeanF = GastValues(:,4);
                            TranscTime{i} = [TranscTime{i},PeakT(2)];
                            GastTime1{i} = [GastTime1{i},GastT(1)];
                            GastTime2{i} = [GastTime2{i},GastT(2)];
                            GastTime3{i} = [GastTime3{i},GastT(3)];
                            Increase{i} = [Increase{i}, PeakF(3)./PeakF(2)];
                            PeakF1{i} = [PeakF1{i},MeanF(1)];
                            %TranscTimeAll{i} = [TranscTimeAll, TranscTime];
                            %GastTime2All{i} = [GastTime2All,GastTime2];
                            %end
                            NewTable = table(string(File),string(Nickname),Rep,GastT(1),GastT(2),GastT(3),PeakT(2),'VariableNames', {'File','Nickname','Repeat','GastT1', 'GasT2', 'GastT3','TranscT2'});
                            Table = [Table;NewTable];
                        end 
                    %end

              figure(Fig)      
              subplot(231); hold on
            PlotScatter(GastTime1{i},TranscTime{i},SelectedNames{i}, 'Start of gastrulation (min into nc14)', 'Transition point (min)',1,0)

             subplot(232); hold on 
             PlotScatter(GastTime2{i},TranscTime{i},SelectedNames{i}, 'Start of ME invagination (min into nc14)', 'Transition point (min)',1,0)

             subplot(233); hold on 
             PlotScatter(GastTime3{i},TranscTime{i},SelectedNames{i}, 'End of gastrulation (min into nc14)', 'Transition point (min)',1,0)

              subplot(234); hold on    
             PlotScatter(GastTime2{i}-GastTime1{i},Increase{i},SelectedNames{i}, 'Length 1st gast half (min)','Increase in levels (AU)',1,0)

             subplot(235); hold on  
             PlotScatter(GastTime3{i}-GastTime2{i},Increase{i},SelectedNames{i}, 'Length 2nd gast half (min)','Increase in levels (AU)',1,0)

             subplot(236); hold on      
             PlotScatter(GastTime3{i}-GastTime1{i},Increase{i},SelectedNames{i}, 'Gastrulation length (min)','Increase in levels (AU)',1,0)


            %  figure(Fig2)
            %  subplot(1,10,[1:4]); hold on
            %   PlotScatter(GastTime2{i},TranscTime{i},SelectedNames{i}, 'Start of ME invagination (min into nc14)', 'Transition in levels (min)',1,0)

            end
             figure(Fig)

              subplot(231); hold on
            PlotScatter([GastTime1{:}],[TranscTime{:}],'All', 'Start of gastrulation (min into nc14)', 'Transition point (min)',0,1)

             subplot(232); hold on 
             PlotScatter([GastTime2{:}],[TranscTime{:}],'All', 'Start of ME invagination (min into nc14)', 'Transition point (min)',0,1)

             subplot(233); hold on 
             PlotScatter([GastTime3{:}],[TranscTime{:}],'All', 'End of gastrulation (min into nc14)', 'Transition point (min)',0,1)

            print(Fig,[ToSave,'_GastCorr',strjoin(Selections),'.pdf'],'-fillpage', '-dpdf');

             figure(Fig2)
             %subplot(1,10,[1:4]); hold on 
             %PlotScatter([GastTime2{:}],[TranscTime{:}],'All', 'Peak of gastrulation (min into nc14)', 'Transition in levels (min)',0,1)


            GastTime3 = Cell2Mat(GastTime3);
            GastTime2 = Cell2Mat(GastTime2);
            GastTime1 = Cell2Mat(GastTime1);
            PeakF1 = Cell2Mat(PeakF1);
            Jitter = 0.5; %Jitter./2 cant be > BarW 
            BarW = 0.35; FaceAlpha = 0.3; DotSize = 15; LineWidth = 2; FontSizeTitle = 16; FontSize = 14;
            set(gca,'FontSize',FontSize)

            subplot(151); hold on     
            plotBoxplot([GastTime2],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[30,70])
            ylabel('Start of ME invagination (min)')

            subplot(152); hold on     
            plotBoxplot([GastTime2-GastTime1],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,30])
            ylabel('Duration of apical constriction (min)')

            subplot(153); hold on     
            plotBoxplot([GastTime3-GastTime2],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,30])
            ylabel('Duration of ME invagination (min)')

            subplot(154); hold on     
            plotBoxplot([GastTime3-GastTime1],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[10,45])
            ylabel('Duration of gastrulation (min)')


            subplot(155); hold on     
            plotBoxplot([PeakF1],SelectedNames,SelectedNames,Jitter,BarW,'',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,700])
            ylabel('Mean Fluorescence 30-50'' (AU)')


            print(Fig2,[ToSave,'_GastCorr',strjoin(Selections),'_min.pdf'],'-fillpage', '-dpdf');
            writetable(Table,[ToSave,'_GastCorr',strjoin(Selections),'_table.txt']);

            close all

end