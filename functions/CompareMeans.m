function[] = CompareMeans(SelectedN,Which,Info, PathPlots,Exps, Selections, MeansOnly, XLim,PalettePoints, PaletteMeans, varargin)
    mkdir(PathPlots)
    Heatmaps = 3000;
    PlotS = 0;
    WhichTime = 'nc14';
    WhichComp = 'MCP';
    Confirm = 0;
    SplitEarlyOverride = NaN;
    try
        if ~isempty(varargin{1})
            YLimits = varargin{1};
        end
        if ~isempty(varargin{2})
            Heatmaps = varargin{2};
        end
        if ~isempty(varargin{3})
            SplitEarlyOverride = varargin{3}; 
        end
        if ~isempty(varargin{4})
            PlotS = varargin{4}; 
        end
        if ~isempty(varargin{5})
            WhichTime = varargin{5};
        end
        if ~isempty(varargin{6})
            Confirm = varargin{6};
        end
        if ~isempty(varargin{7})
            WhichComp = varargin{7};
        end
    end
    %try
    %Selection = '';
    if MeansOnly == 1
        Fig1 = figure('PaperSize',[30 50],'PaperUnits','inches','resize','on', 'visible','on');
        Fig1.Renderer='Painters';
    else
        Fig2 = figure('PaperSize',[30 50],'PaperUnits','inches','resize','on', 'visible','on');
        Fig2.Renderer='Painters';
        
    end
    
    Fig3 = figure('PaperSize',[17 50],'PaperUnits','inches','resize','on', 'visible','on');
    Fig3.Renderer='Painters';

     
    if Which == 3
        Index = find(cellfun(@(x) strcmp(x,SelectedN{1}),table2array(Exps(:,1)))==1)';
        SelectedN = table2cell(Exps(Index,3));
    end
    
     ToSave = [PathPlots,char(join(SelectedN,'_vs_'))];
        FigH = figure('PaperUnits','inches','PaperSize',[4*length(SelectedN) 10],'Units','points','resize','on', 'visible','on','DefaultAxesFontSize', 7);
        FigH.Renderer='Painters';
  
    for i = 1:length(SelectedN)
        Index = find(cellfun(@(x) strcmp(x,SelectedN{i}),table2array(Exps(:,Which)))==1)';
        PropertiesMerged = table();
        TimeScaleMerged = [];
        NormMerged = [];
        MeanNormMerged =[];
        OnOffMerged = [];
        SizeEach = [];

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
            if strcmp(WhichComp,'MCP')
                PathData = [PathToSave,'_Data.mat'];
                load(PathData);
                Struct2Vars(Data);
                MeasuredF = MaxF;
            elseif strcmp(WhichComp,'PCP')
                PathData = [PathToSave,'_PCP_Data.mat']
                load(PathData);
                Struct2Vars(Data);
                MeasuredF = MaxF;
            elseif strcmp(WhichComp,'mean')
                PathData = [PathToSave,'_mean_Data.mat'];
                load(PathData);
                Struct2Vars(Data);
                MeasuredF = MeanF;
            end
            Parameters = Info(Index(x),:);
            Table2Vars(Parameters);
            
            minOn = 5;
            Selected = [];
              
            if ~isnan(SplitEarlyOverride)
                SplitEarly = SplitEarlyOverride;
                Labels = Properties.NewLabel;
                FoldOverBaseline = 1.2;
                [~, ~,PropertiesNew, ~] = DefineExpAll(MeasuredF,MaxFBG,CentX,Labels,Baseline,FoldOverBaseline,TimeRes,3,60,SplitEarly,nc14,Delay,minOn);
                %close Fig
                Properties.Onset = PropertiesNew.Onset;
                Properties.Type = PropertiesNew.Type;%%##########################
            end
            
            
            SelectionToSave = Selection;
            if strcmp(Selection,'') == 1
                Selected = [Properties.Type ~= 'BG'];
            else
                Selected = Properties.Type ~= 'EarlyOnly' & Properties.Type ~= 'LateTrack' & Properties.Region == Selection;
                if contains(Selection,'|')
                    Selected = [];
                    Selections = strsplit(Selection,'|');
                    SelectionToSave = join(Selections,'');SelectionToSave = SelectionToSave{:};
                    for s = 1:length(Selections)
                    	Selected(s,:) = Properties.Type ~= 'EarlyOnly' & Properties.Type ~= 'LateTrack' & Properties.Region == Selections{s};
                    end
                    Selected = [sum(Selected,1) ~= 0]';      
                end    
            end

            
           
            XLabel = 'time into nc14 (min)';
            if strcmp(WhichTime,'Gast') && (sum(Selected) > 0)
                disp('Using gastrulation time as t = 0')
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
                    TimeScale = TimeScale - round(MEinvT.*60./TimeRes).*TimeRes./60;
                    XLim = [-40,15];
                    %XLim = [-55,15];
                    XLabel = 'time from mesoderm invagination (min)';
            end
            
            
            SplitEarlyF = max([SplitEarly*60./TimeRes+nc14-Delay,1]);
            OnOff = CleanOnOff(OnOff,minOn);
            [OnOff] = CleanNaNs(MedFilt,OnOff, minOn*2);
            Norm = ((MeasuredF-Baseline').*Baseline(1)./Baseline').*(2.^(12-Bits));
            %Norm = (MeasuredF./Baseline').*(2.^(12-Bits));
            if strcmp(WhichComp,'mean')
                [~,Time0] = min(abs(TimeScale));
                Norm = (MeasuredF)./nanmean(MeasuredF(max(1,Time0),:)); % #############################
                Norm = (MeasuredF).*(2.^(12-Bits)); % #############################
                %Norm = (MeasuredF-Baseline').*Baseline(1)./Baseline'.*(2.^(12-Bits));
            end
            Norm(Norm==0) = NaN;
            Limits = [0, 90];
            
            %uncomment to plot selected for each repeat
            if PlotS
                mkdir([PathPlots,'/PlotsSelected/'])
                [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll,BurstSize] = CountBursts(Norm,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
                PlotSelected(Selected,TimeScale,Baseline,MeasuredF,MedFilt,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,Limits,[PathPlots,'/PlotsSelected/',File,'_selected_',Selection,'.ps'],{});
            end

             % count only cells with onsets in the time being plotted
            %Selected = [PropertiesMerged.Type ~= 'BG'];%##############
               if strcmp(WhichTime,'Gast'); 
                    Selected = Selected & (Properties.Type ~= 'BG' & Properties.Type ~= 'LateTrack' & Properties.Onset < (XLim(2)+MEinvT));
                    %Selected = [PropertiesMerged.Type ~= 'BG']; %################
               else
                    Selected = Selected & (Properties.Type ~= 'BG' & Properties.Type ~= 'LateTrack' & Properties.Onset < XLim(2));
               end
 
            MeasuredF = MeasuredF(:,Selected);
            MedFilt = MedFilt(:,Selected);
            OnOff = OnOff(:,Selected);
            Properties = Properties(Selected,:);
            Properties.NormAP = (Properties.AP_position-min(Properties.AP_position))./max(Properties.AP_position);
            Norm = Norm(:,Selected);
            OnOff(isnan(MeasuredF)) = NaN;
            MeanNormEach = nanmean(Norm,2);
            SizeEach(x) = size(Norm,2);

            [NormMerged,~,~] = MergeFMatrix(Norm,NormMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     
            [MeanNormMerged,~,~] = MergeFMatrix(MeanNormEach,MeanNormMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     
            [OnOffMerged,PropertiesMerged,TimeScaleMerged] = MergeFMatrix(OnOff,OnOffMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     
        end
        
        TimeResMerged = (TimeScaleMerged(2)-TimeScaleMerged(1)).*60;
        %NormMergedON = (NormMerged .* OnOffMerged);
        SplitEarlyFMerged = find(round(TimeScaleMerged,2) == SplitEarly);
        EndMerged = find(round(TimeScaleMerged,2) == XLim(2));
        if strcmp(WhichTime,'Gast'); 
            SplitEarlyFMerged = find(round(TimeScaleMerged,2) == (SplitEarly-round(MEinvT)));
            EndMerged = find(round(TimeScaleMerged,2) == XLim(2));
        end
        if isempty(EndMerged); EndMerged = length(TimeScaleMerged);end
        PropertiesMerged.TotalmRNA = [nansum(NormMerged(SplitEarlyFMerged:min(size(NormMerged,1),EndMerged),:),1)]';

%        [NormMerged,~,~] = MergeFMatrix(nan(3,1),NormMerged,table(0),PropertiesMerged,[0:2].*TimeResMerged./60,TimeScaleMerged,TimeResMerged);
%        [MeanNormMerged,~,~] = MergeFMatrix(nan(3,1),MeanNormMerged,table(0),PropertiesMerged,[0:2].*TimeResMerged./60,TimeScaleMerged,TimeResMerged);
%        [OnOffMerged,PropertiesMerged,TimeScaleMerged] = MergeFMatrix(nan(3,1),OnOffMerged,table(0),PropertiesMerged,[0:2].*TimeResMerged./60,TimeScaleMerged,TimeResMerged);
%         NormMerged(:,end) = [];
%         MeanNormMerged(:,end) = [];
%         OnOffMerged(:,end) = [];
%         PropertiesMerged(end,:) = [];

       % ColorArg is 8 value colormap with info for allPoints, allMidPoints,
       % shortPoints, longPoints, allMeans, allMidMeans, shortMeans, longMeans
       % repeat value for each 4 times, it will update with every repeat
       ColorArg = [PalettePoints(i,:);PalettePoints(i,:);PalettePoints(i,:);PalettePoints(i,:);...
       PaletteMeans(i,:);PaletteMeans(i,:);PaletteMeans(i,:);PaletteMeans(i,:)];
       
                
       YLimits = [min(0,min(MeasuredF(:))), max(2^Bits-1,max(MeasuredF(:)))];
        try
            YLimits = varargin{1};
        end  
        Selected = [PropertiesMerged.Type ~= 'EarlyOnly'];
       if MeansOnly == 1
            Fig1 = PlotMeansFractionShaded(NormMerged,MeanNormMerged,OnOffMerged,TimeScaleMerged,Selected,PropertiesMerged,Bits,ColorArg,Fig1,1,SelectedN{i},XLim, YLimits,XLabel);
       else
            Fig2 = PlotMeans(NormMerged,OnOffMerged,TimeScaleMerged,Selected,PropertiesMerged,Bits,ColorArg,Fig2,0,SelectedN{i},XLim, YLimits);
       end
     
    %end
    if Heatmaps ~= 0
        MaxTime = 30;
        %Selected = [PropertiesMerged.Type ~= 'EarlyOnly'& Selected];
        Selected = [PropertiesMerged.Type ~= 'EarlyOnly'];
        [NormAligned,TimeScaleAligned] = AlignFMatrixtoOnset(NormMerged,PropertiesMerged,TimeScaleMerged,MaxTime,TimeResMerged);
        [OnOffAligned,~] = AlignFMatrixtoOnset(OnOffMerged,PropertiesMerged,TimeScaleMerged,MaxTime,TimeResMerged);
        Fig3 = PlotMeansFractionShaded(NormAligned,nanmean(NormAligned,2),OnOffAligned,TimeScaleAligned,Selected,PropertiesMerged,Bits,ColorArg,Fig3,1,SelectedN{i},[0,MaxTime], YLimits,'time from onset (min)');
        %Fig3 = PlotMeans(NormAligned,TimeScaleAligned,Selected,PropertiesMerged,Bits,ColorArg,Fig3,1,SelectedN{i},[0,MaxTime], YLimits);
        [FigH] = PlotHeatmaps(FigH,NormMerged,Selected,PropertiesMerged,TimeScaleMerged,SelectedN,TimeResMerged,XLim,i,Heatmaps,XLabel,SizeEach);
    end
    
    
    end
    
      
    %pause(30)
    Selection = strjoin(Selections);
    if strcmp(Selection,'') == 0; Selection = ['_',Selection];end
     if strcmp(WhichTime,'Gast'); Selection = [Selection,'_Gast'];end
    if MeansOnly == 1
        print(Fig1,[ToSave,Selection,'_means.pdf'],'-fillpage', '-dpdf');
    else
        print(Fig2,[ToSave,Selection,'.pdf'],'-fillpage', '-dpdf');
    end
    

    if Heatmaps ~= 0
        print(Fig3,[ToSave,Selection,'_aligned.pdf'],'-fillpage', '-dpdf');
        print(FigH,[ToSave,Selection,'_heatmaps.pdf'],'-fillpage', '-dpdf');
    end
    close all
    %end
end