function[] = CompareGastrulation(SelectedN,Which,Info, PathPlots,Exps, Selections, XLim, Palette, varargin)
    YLimits = [0,0.15];
    WhichTime = 'nc14';
    Confirm = 0;
    SplitEarlyOverride = NaN
    try
        if ~isempty(varargin{1})
            YLimits = varargin{1};
        end
        if ~isempty(varargin{2})
            SplitEarlyOverride = varargin{2}; 
        end
        if ~isempty(varargin{3})
            WhichTime = varargin{3};
        end
        if ~isempty(varargin{4})
            Confirm = varargin{4};
        end

    end
        if strcmp(WhichTime,'Gast')
            Fig1 = figure('PaperSize',[15 50],'PaperUnits','inches','resize','on', 'visible','on');
            Fig1.Renderer='Painters';
        else
            Fig1 = figure('PaperSize',[30 50],'PaperUnits','inches','resize','on', 'visible','on');
            Fig1.Renderer='Painters';
        end
  
    if Which == 3
        Index = find(cellfun(@(x) strcmp(x,SelectedN{1}),table2array(Exps(:,1)))==1)';
        SelectedN = table2cell(Exps(Index,3));
    end
    
     ToSave = [PathPlots,char(join(SelectedN,'_vs_'))];
       
     
    for i = 1:length(SelectedN)
        Index = find(cellfun(@(x) strcmp(x,SelectedN{i}),table2array(Exps(:,Which)))==1)';
        PropertiesMerged = table();
        TimeScaleMerged = [];
        NormMerged = [];
        MeanNormMerged =[];
        OnOffMerged = [];
        GastTimes = [];
        if length(Selections) == length(SelectedN)
            Selection = Selections{i}
        else
            Selection = Selections{1};
        end
        %try
        for x = 1:size(Index,2)
            %try
            Experiment = [Exps.Nickname{Index(x)},' ',num2str(Exps.Rep(Index(x)))]
            PathToSave = [Info.Path{Index(x)},Info.File{Index(x)},...
            Info.Name{Index(x)},Info.File{Index(x)}]; 
            PathData = [PathToSave,'_Data.mat'];
            load(PathData);
            Struct2Vars(Data);
            Parameters = Info(Index(x),:);
            Table2Vars(Parameters);
            try SplitEarly = SplitEarlyOverride; end
            
            SelectionToSave = Selection;
            if strcmp(Selection,'') == 1
                Selected = [Properties.Type ~= 'BG'];
            else
                Selected = Properties.Type ~= 'EarlyOnly' & Properties.Region == Selection;
                if contains(Selection,'|')
                    Selected = [];
                    Selections = strsplit(Selection,'|');
                    SelectionToSave = join(Selections,'');SelectionToSave = SelectionToSave{:};
                    for s = 1:length(Selections)
                    	Selected(s,:) = Properties.Type ~= 'EarlyOnly' & Properties.Region == Selections{s};
                    end
                    Selected = sum(Selected,1) ~= 0;      
                end    
            end
        
            XLabel = 'time into nc14 (min)';
            if strcmp(WhichTime,'Gast')
                %disp('Using gastrulation time as t = 0')
                clear MEinvT
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
                    MEinvT = GastValues(2,3);
                    TimeScale = TimeScale - round(MEinvT);
                    XLim = [-5,20];
                    XLabel = 'time from mesoderm invagination (min)';
            end
            minOn = 5;
            if ~isnan(SplitEarlyOverride)
                Labels = Properties.NewLabel;
                FoldOverBaseline = 1.2;
                [~, ~,PropertiesNew, ~] = DefineExpAll(MaxF,MaxFBG,CentX,Labels,Baseline,FoldOverBaseline,TimeRes,3,60,SplitEarly,nc14,Delay,minOn);
                %close Fig
                Properties.Onset = PropertiesNew.Onset;
                Properties.Type = PropertiesNew.Type;
            end

                MovX = abs(diff(CentX,1,1))*XYRes;
                MovY = abs(diff(CentY,1,1))*XYRes;
                MovZ = abs(diff(CentZ,1,1))*ZRes;
                MovXYZ = sqrt(MovX.^2+MovY.^2+MovZ.^2);
                Mov = MovY./TimeRes;
                Mov(end+1,:) = NaN;
                
            OnOff = CleanOnOff(OnOff,minOn);
            [OnOff] = CleanNaNs(MedFilt,OnOff, minOn*2);

            %Limits = [0, 90];
            
            % uncomment to plot selected for each repeat
%             [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll,BurstSize] = CountBursts(Norm,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
%             PlotSelected(Selected, TimeScale,Baseline,MaxF,MedFilt,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,Limits,[PathPlots,'/',Info.File{Index(x)},'_selected_',Selection,'.ps']);
% %         
            Norm = Mov;
            YLabel = 'Mean DV displacement (um/s)';
            Norm = -(CentY-nanmean(CentY(min(Frames,nc14-Delay+60*60/TimeRes):end,:),1)).*XYRes;
            YLabel = 'Distance from midline (um)';
            Norm = -(CentY-min(CentY,[],1,'omitnan')).*XYRes;
            YLabel = 'Relative DV position (um)';
            %Norm = -(CentY-nanmean(CentY(max(1,nc14-Delay+15*60/TimeRes):max(1,nc14-Delay+30*60/TimeRes),:),1)).*XYRes;
            %YLabel = 'Distance from MSE (um)';
            YLimits = [-inf,inf];

            
            MaxF = MaxF(:,Selected);
            OnOff = OnOff(:,Selected);
            Properties = Properties(Selected,:);
            Properties.NormAP = (Properties.AP_position-min(Properties.AP_position))./max(Properties.AP_position);
            Norm = Norm(:,Selected);
            OnOff(isnan(MaxF)) = NaN;
            MeanNormEach = nanmean(Norm,2);
            
            [NormMerged,~,~] = MergeFMatrix(Norm,NormMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     
            [MeanNormMerged,~,~] = MergeFMatrix(MeanNormEach,MeanNormMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     
            [OnOffMerged,PropertiesMerged,TimeScaleMerged] = MergeFMatrix(OnOff,OnOffMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     
            %GastTimes(:,:,x) = GastT;
        end
        
%        [NormMerged,~,~] = MergeFMatrix([NaN],NormMerged,table(0),PropertiesMerged,[0],TimeScaleMerged,TimeRes);
%        [MeanNormMerged,~,~] = MergeFMatrix([NaN],MeanNormMerged,table(0),PropertiesMerged,[0],TimeScaleMerged,TimeRes);
%        [OnOffMerged,PropertiesMerged,TimeScaleMerged] = MergeFMatrix([NaN],OnOffMerged,table(0),PropertiesMerged,[0],TimeScaleMerged,TimeRes);
%         MeanNormMerged(:,end) = [];
%         
       % ColorArg is 8 value colormap with info for allPoints, allMidPoints,
       % shortPoints, longPoints, allMeans, allMidMeans, shortMeans, longMeans
       % repeat value for each 4 times, it will update with every repeat
       ColorArg = [Palette(i,:);Palette(i,:);Palette(i,:);Palette(i,:);...
       Palette(i,:);Palette(i,:);Palette(i,:);Palette(i,:)];
       Selected = [PropertiesMerged.Type ~= 'BG'];

            Fig1 = PlotMeansFractionShaded(NormMerged,MeanNormMerged,OnOffMerged,TimeScaleMerged,Selected,PropertiesMerged,Bits,ColorArg,Fig1,1,SelectedN{i},XLim, YLimits,XLabel);
        
            % change plot settings from default PlotMeansFractionShaded
            figure(Fig1)
            subplot(411); ylabel(YLabel)
            subplot(412); ylabel(YLabel)
            	%plot([GastT,GastT],[0,400],'--','Color',PalettePoints(i,:))
                %MeanGastT = nanmean(GastTimes,3);
                %SEMGast = nanstd(GastTimes,[],3)./sqrt(sum(~isnan(GastTimes),3));
%                 t=1;
%                 plot([MeanGastT(t),MeanGastT(t)],[0,YLimits(2)],'--','Color',PalettePoints(i,:),'HandleVisibility','off')
%                 Poly = polyshape([MeanGastT(t)-SEMGast(t),MeanGastT(t)-SEMGast(t),MeanGastT(t)+ SEMGast(t),MeanGastT(t)+ SEMGast(t)],[0,YLimits(2),YLimits(2),0]);
%                 plot(Poly,'LineStyle','none','FaceColor',[PalettePoints(i,:),0.3],'HandleVisibility','off');
%                 t=2;
%                 plot([MeanGastT(t),MeanGastT(t)],[0,YLimits(2)],'-','Color',PalettePoints(i,:),'HandleVisibility','off')
%                 Poly = polyshape([MeanGastT(t)-SEMGast(t),MeanGastT(t)-SEMGast(t),MeanGastT(t)+ SEMGast(t),MeanGastT(t)+ SEMGast(t)],[0,YLimits(2),YLimits(2),0]);
%                 plot(Poly,'LineStyle','none','FaceColor',[PalettePoints(i,:),0.3],'HandleVisibility','off');
%                 t=3;
%                 plot([MeanGastT(t),MeanGastT(t)],[0,YLimits(2)],':','Color',PalettePoints(i,:),'HandleVisibility','off')
%                 Poly = polyshape([MeanGastT(t)-SEMGast(t),MeanGastT(t)-SEMGast(t),MeanGastT(t)+ SEMGast(t),MeanGastT(t)+ SEMGast(t)],[0,YLimits(2),YLimits(2),0]);
%                 plot(Poly,'LineStyle','none','FaceColor',[PalettePoints(i,:),0.3],'HandleVisibility','off');
%                     
            subplot(413); ylabel(YLabel)
            subplot(414); ylabel(YLabel)

    
    
    end
    Selection = strjoin(Selections);
    if strcmp(Selection,'') == 0; Selection = ['_',Selection];end
    if strcmp(WhichTime,'Gast'); Selection = [Selection,'_Gast'];end
    print(Fig1,[ToSave,Selection,'_gast.pdf'],'-fillpage', '-dpdf');
   
    
    close all
    %end


    
end
