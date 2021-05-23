function[] = CompareIO(SelectedN,Which,Info, PathPlots,Exps, Selections, XLim,PalettePoints, PaletteMeans, varargin)
    mkdir(PathPlots)
    SplitEarlyOverride = NaN;
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
    %try
    %Selection = '';
        Fig1 = figure('PaperSize',[20 20],'PaperUnits','inches','resize','on', 'visible','on');
        Fig1.Renderer='Painters';


    if Which == 3
        Index = find(cellfun(@(x) strcmp(x,SelectedN{1}),table2array(Exps(:,1)))==1)';
        SelectedN = table2cell(Exps(Index,3));
    end
    
     ToSave = [PathPlots,char(join(SelectedN,'_vs_'))];
  
    for i = 1:length(SelectedN)
        Index = find(cellfun(@(x) strcmp(x,SelectedN{i}),table2array(Exps(:,Which)))==1)';
        
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
            Parameters = Info(Index(x),:);
            Table2Vars(Parameters);
            
            minOn = 5;
            Selected = [];
            
%            if strcmp(WhichComp,'MCP')
               
 %           elseif strcmp(WhichComp,'PCP')
                %PathData = [PathToSave,'_PCP_Data.mat']
                %load(PathData);
                %Struct2Vars(Data);
                %MeasuredF = MaxF;
           % elseif strcmp(WhichComp,'mean')
                PathData = [PathToSave,'_mean_Data.mat'];
                load(PathData);
                Struct2Vars(Data);            
                NormMeanF = (MeanF-Baseline').*Baseline(1)./Baseline'.*(2.^(12-Bits));
                PropertiesMean = Properties;
            %end
             PathData = [PathToSave,'_Data.mat'];
                load(PathData);
                Struct2Vars(Data);
                MCPMaxF = MaxF;
                 SplitEarlyF = max([SplitEarly*60./TimeRes+nc14-Delay,1]);
                OnOff = CleanOnOff(OnOff,minOn);
                [OnOff] = CleanNaNs(MedFilt,OnOff, minOn*2);
                Norm = (MaxF-Baseline').*Baseline(1)./Baseline'.*(2.^(12-Bits));
                PropertiesMCP = Properties;
            

              
            if ~isnan(SplitEarlyOverride)
                SplitEarly = SplitEarlyOverride;
                Labels = PropertiesMCP.NewLabel;
                FoldOverBaseline = 1.2;
                [~, ~,PropertiesNew, ~] = DefineExpAll(MeasuredF,MaxFBG,CentX,Labels,Baseline,FoldOverBaseline,TimeRes,3,60,SplitEarly,nc14,Delay,minOn);
                %close Fig
                PropertiesMCP.Onset = PropertiesNew.Onset;
                PropertiesMCP.Type = PropertiesNew.Type;%%##########################
            end
            
            
            SelectionToSave = Selection;
            if strcmp(Selection,'') == 1
                Selected = [PropertiesMCP.Type ~= 'BG'];
            else
                Selected = PropertiesMCP.Type ~= 'EarlyOnly' & PropertiesMCP.Type ~= 'LateTrack' & PropertiesMCP.Region == Selection;
                if contains(Selection,'|')
                    Selected = [];
                    Selections = strsplit(Selection,'|');
                    SelectionToSave = join(Selections,'');SelectionToSave = SelectionToSave{:};
                    for s = 1:length(Selections)
                    	Selected(s,:) = PropertiesMCP.Type ~= 'EarlyOnly' & PropertiesMCP.Type ~= 'LateTrack' & PropertiesMCP.Region == Selections{s};
                    end
                    Selected = [sum(Selected,1) ~= 0]';      
                end    
            end

            figure(Fig1)
           Onsets =  PropertiesMCP.Onset(Selected);
           MeanMeanFSelected = NormMeanF(:,Selected);
           MeanFwhenOnset = MeanMeanFSelected(sub2ind(size(MeanMeanFSelected),Onsets*60/TimeRes+nc14-Delay,[1:length(Onsets)]'));
           plot(PropertiesMCP.Onset(Selected),MeanFwhenOnset,'*','Color',PalettePoints(i,:)); hold on
           box off
            xlabel('Onset time (min into nc14)')
            ylabel('Normalized NICD levels (AU)')
            xlim(XLim)
            ylim(YLimits)
            
           mkdir([PathPlots,'/PlotsIO/'])
           [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll,BurstSize] = CountBursts(Norm,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
           PlotSelected(Selected,TimeScale,Baseline,MCPMaxF,NormMeanF,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,XLim,[PathPlots,'/PlotsIO/',File,'_IO_',Selection,'.ps'],NormMeanF);

        end
            
    
    end
    
      
    %pause(30)
    Selection = strjoin(Selections);
    if strcmp(Selection,'') == 0; Selection = ['_',Selection];end
     if strcmp(WhichTime,'Gast'); Selection = [Selection,'_Gast'];end

    print(Fig1,[ToSave,Selection,'_scatterIO.pdf'],'-fillpage', '-dpdf');

    close all
    %end
end