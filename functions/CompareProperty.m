function[] = CompareProperty(SelectedN,Which,Info, PathPlots,Exps, Selection, XLim,PalettePoints, PaletteMeans, varargin)
    mkdir(PathPlots)    
    SplitEarly = 15;
    WhichTime = 'nc14';
    Confirm = 0;
    try
        if ~isempty(varargin{1})
            YLimits = varargin{1};
        end
        if ~isempty(varargin{2})
            Var = varargin{2};
        end
        if ~isempty(varargin{3})
            SplitEarly = varargin{3};
        end
        if ~isempty(varargin{4})
            WhichTime = varargin{4};
        end
        if ~isempty(varargin{5})
            Confirm = varargin{5};
        end
    end
    %try
    %Selection = '';
        Fig1 = figure('PaperSize',[20 50],'PaperUnits','inches','resize','on', 'visible','on');
        Fig1.Renderer='Painters';
     
    if Which == 3
        Index = find(cellfun(@(x) strcmp(x,SelectedN{1}),table2array(Exps(:,1)))==1)';
        SelectedN = table2cell(Exps(Index,3));
    end
    
         ToSave = [PathPlots,char(join(SelectedN,'_vs_'))];
 
  
    for i = 1:length(SelectedN)
        Index = find(cellfun(@(x) strcmp(x,SelectedN{i}),table2array(Exps(:,Which)))==1)';
        PropertiesMerged = table();
        TimeScaleMerged = [];
        %NormMerged = [];
        PropertyMerged = [];
        MeanEachMerged = [];
        %MeanNormMerged =[];
        %OnOffMerged = [];
        if length(Selections) == length(SelectedN)
            Selection = Selections{i}
        else
            Selection = Selections{1};
        end
        %try
        for x = 1:size(Index,2)
            %try
            Experiment = [Exps.Nickname{Index(x)},' ',num2str(Exps.Rep(Index(x)))];
            PathToSave = [Info.Path{Index(x)},Info.File{Index(x)},...
            Info.Name{Index(x)},Info.File{Index(x)}]; 
            load([PathToSave,'_Data3D.mat']);
            %disp(['loaded ',PathToSave,'_Data3D.mat'])
            Property = Data3D.(Var);
            %size(Property)
            Parameters = Info(Index(x),:);
            Table2Vars(Parameters);
            try
                load([PathToSave,'_Data.mat']);
                Struct2Vars(Data);
            catch
                TimeScale = ([1:Frames]-nc14+Delay+From-1)*TimeRes/60;
                %TimeScale = TimeScale(From:end);
                %size(TimeScale)
            end
            
            if strcmp(Var,'Volume')
                Factor = XYRes.^3;
                Units = 'um^3';
            elseif  strcmp(Var,'SurfaceArea') | ...
                   strcmp(Var,'Area1') | strcmp(Var,'Area2') |strcmp(Var,'Area3')
                Factor = XYRes.^2;
                Units = 'um^2';
            elseif strcmp(Var,'PrincipalAxisLength_1') | strcmp(Var,'PrincipalAxisLength_2') | strcmp(Var,'PrincipalAxisLength_3') |... 
                strcmp(Var,'Perimeter1') | strcmp(Var,'Perimeter2') | strcmp(Var,'Perimeter3') | ...
                strcmp(Var,'EquivDiameter1') | strcmp(Var,'EquivDiameter2') | strcmp(Var,'EquivDiameter3') 
                Factor = XYRes;
                Units = 'um';
            else
                Factor = 1;
                Units = 'dimensionless';
            end
            
            Property = Property.*Factor;
%             if SplitEarly ~= 15
%                 %SplitEarly =30; 
%                 minOn = 5;
%                 Labels = Properties.NewLabel;
%                 [~, ~,PropertiesNew, ~] = DefineExpAll(MaxF,MaxFBG,CentX,Labels,Baseline,TimeRes,3,60,SplitEarly,nc14,Delay,minOn);
%                 %close Fig
%                 Properties.Onset = PropertiesNew.Onset;
%                 Properties.Type = PropertiesNew.Type;
%             end
           
           Properties = table();
           Properties.Type = string(repmat('All',size(Property,2),1));
           Properties.Label = [1:size(Property,2)]';

            
            
            
            F = max(round(35*60/TimeRes)+nc14-Delay,1);
        if F >= Frames ; F = Frames-1; end
        
        
        try
            ImLab = CentroidsF;
        catch
            load([PathToSave, '_Stats.mat']);
            try
                Merged_meanF_maxGFP = Read3d([PathToSave, '_maxF_maxGFP.tiff']);
                Im = Merged_meanF_maxGFP(:,:,F)./255;
            catch
                Merged_meanF_maxGFP = Read3dRGB([PathToSave, '_segmented_tracked_boundaries_RGB_clean.tiff']);
                Im = Merged_meanF_maxGFP(:,:,1,F)./255;
            end
            Width = size(Im,1);
            Height = size(Im,2);
             % rescue boundariesL from StatsGFP
            cmap = jet(1000000);
            cmap_shuffled = cmap(randperm(size(cmap,1)),:);
            boundariesL = zeros(Height,Width,Frames);
            for f = 1:length(Stats_GFP)
                 disp(num2str(f))
                 SingleFrame = Stats_GFP{f};
                 Tracked3D = zeros(Height,Width,Slices);
                 for n = 1:size(SingleFrame,1)
                     Tracked3D(SingleFrame.PixelIdxList{n}) =  SingleFrame.Label(n);
                 end
                 %toTrack(:,:,:,f) = Tracked3D;
                 TrackedBoundL = MAX_proj_3D(Tracked3D);
                 [BoundariesBW, ~, ~] = BoundariesTracked_3D(Tracked3D,cmap_shuffled,'off');
                 boundariesL(:,:,f) = BoundariesBW .* TrackedBoundL;
            end   
            
            ImLab = boundariesL(:,:,max(1,F-5):min(F+5,Frames-1));
            for f = 1:size(ImLab,3)
                Stats{f} = regionprops('table',ImLab(:,:,f),'Centroid');
                Stats{f}.Label = [1:height(Stats{f})]';
            end
            [AllF] = MergeAll(Stats, TimeRes);
            AllF = splitvars(AllF);
            Labels = unique(AllF.Label); LabelsOld = Labels;
            [CentXF] = Reshape(AllF,Frames,Labels,'Centroid_1','Label');
            [CentYF] = Reshape(AllF,Frames,Labels,'Centroid_2','Label');
            PosX = nanmean(CentXF,1);
            PosY = nanmean(CentYF,1);
             PosX(isnan(PosX)) = [];
             PosY(isnan(PosY)) = [];
%             StatsF = Stats_GFP{F};
%             CentXF = StatsF.Centroid(:,1);
%             CentYF = StatsF.Centroid(:,2);
            Indices = sub2ind([size(Im,1),size(Im,2)],floor(PosY), floor(PosX));
            ImLab = zeros(size(Im,1),size(Im,2));
            ImLab(Indices) = 1;
            Data.CentroidsF = ImLab;
            save([PathToSave,'_Data.mat'],'Data');
        end
        
        
        PathToSaveR = [Path,File,Name,'regions/',File]; 
        Regions = imread([PathToSaveR,'_regions.tiff']);
               
        CentX = Data3D.Centroid_1;
        CentY = Data3D.Centroid_2;
        
        Selected = [1:size(Property,2)];
        %LabelsSelected = Labels(Selected);
        PosX = floor(nanmean(CentX(max(1,F-5): min(F+5,Frames-1),Selected),1));
        PosY = floor(nanmean(CentY(max(1,F-5): min(F+5,Frames-1),Selected),1));
        %Indices = sub2ind(size(Regions),PosY,PosX);
        % added to prevent too many NaNs when cells are not tracked in the
        % F-5:F+5 window
        PosX(isnan(PosX)) = floor(nanmean(CentX(max(1,nc14): min(Frames-1,nc14-Delay+40*60./TimeRes),Selected(isnan(PosX))),1));
        PosY(isnan(PosY)) = floor(nanmean(CentY(max(1,nc14): min(Frames-1,nc14-Delay+40*60./TimeRes),Selected(isnan(PosY))),1));
        Indices = sub2ind(size(Regions),PosY,PosX);
        Selected = Selected(~isnan(Indices));
        %LabelsSelected = LabelsSelected(~isnan(Indices));
        PosX = PosX(~isnan(Indices));
        PosY = PosY(~isnan(Indices));
        Indices = Indices(~isnan(Indices));
        RegionsInd = Regions(Indices);
        ME = Selected(RegionsInd==1);
        MSE = Selected(RegionsInd==2);
        NE = Selected(RegionsInd==3);
        DE = Selected(RegionsInd==4);
        Properties.Region = string(repmat('NaN',length(Properties.Type),1));
        Properties.Region(ME) = string(repmat('ME',length(ME),1));
        Properties.Region(MSE) = string(repmat('MSE',length(MSE),1));
        Properties.Region(NE) = string(repmat('NE',length(NE),1));
        Properties.Region(DE) = string(repmat('DE',length(DE),1));
        
            
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

            TimeScale = ([1:size(Property,1)]-nc14+Delay+From-1)*TimeRes/60;
            SplitEarlyF = max([SplitEarly*60./TimeRes+nc14-Delay,1]);
            minOn = 5;
            Limits = [0, 90];
            
            XLabel = 'time into nc14 (min)';
            if strcmp(WhichTime,'Gast')
                disp('Using gastrulation time as t = 0')
                clear MEinvT
                    try
                        GastValues = Data.(['GastValues',SelectionToSave]);
                        if Confirm
                        	Data = DrawGastTimes(Data,Parameters,Selection,Confirm,SplitEarly);
                            GastValues = Data.(['GastValues',SelectionToSave]);
                        end
                    catch
                        Data = DrawGastTimes(Data,Parameters,Selection,Confirm,SplitEarly);
                        GastValues = Data.(['GastValues',SelectionToSave]);
                    end
                    MEinvT = GastValues(2,3);
                    TimeScale = TimeScale - round(MEinvT);
                    XLim = [-40,20];
                    XLabel = 'time from mesoderm invagination (min)';
            end
            
            Property = Property(:,Selected);
            MeanEach = nanmean(Property,2);
            [MeanEachMerged,~,~] = MergeFMatrix(MeanEach,MeanEachMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     
            [PropertyMerged,PropertiesMerged,TimeScaleMerged] = MergeFMatrix(Property,PropertyMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);     


        end
        
        TimeResMerged = (TimeScaleMerged(2)-TimeScaleMerged(1)).*60;
        SplitEarlyFMerged = find(round(TimeScaleMerged,2) == SplitEarly);
        EndMerged = find(round(TimeScaleMerged,2) == XLim(2));
        if strcmp(WhichTime,'Gast'); 
            SplitEarlyFMerged = find(round(TimeScaleMerged,2) == (SplitEarly-round(MEinvT)));
            EndMerged = find(round(TimeScaleMerged,2) == XLim(2));
        end
        if isempty(EndMerged); EndMerged = length(TimeScaleMerged);end
%         
%         TimeResMerged = (TimeScaleMerged(2)-TimeScaleMerged(1)).*60;
%        [MeanEachMerged,~,~] = MergeFMatrix([NaN],MeanEachMerged,table(0),PropertiesMerged,[0],TimeScaleMerged,TimeResMerged);
%        [PropertyMerged,PropertiesMerged,TimeScaleMerged] = MergeFMatrix([NaN],PropertyMerged,table(0),PropertiesMerged,[0],TimeScaleMerged,TimeResMerged);
%         MeanEachMerged(:,end) = [];
%         PropertyMerged(:,end) = [];
%         PropertiesMerged(end,:) = [];

        ColorArg = [PalettePoints(i,:);PalettePoints(i,:);PalettePoints(i,:);PalettePoints(i,:);...
       PaletteMeans(i,:);PaletteMeans(i,:);PaletteMeans(i,:);PaletteMeans(i,:)];
       Selected = logical(ones(1,size(PropertyMerged,2)));
        try
            YLimits = varargin{1};
        end  
        
        PropertyMerged(PropertyMerged==0) = NaN;
        MeanEachMerged(MeanEachMerged==0) = NaN;
        XLabel = 'time into nc14 (min)';
       Fig1 = PlotMeansFractionShaded(PropertyMerged,MeanEachMerged,zeros(size(PropertyMerged)),TimeScaleMerged,Selected,PropertiesMerged,Bits,ColorArg,Fig1,1,SelectedN{i},XLim, YLimits,XLabel);

    
    end
    
    
     figure(Fig1); hold on
       subplot(411); ylabel([Var,' (',Units,')']);
       subplot(412); ylabel([Var,' (',Units,')']);
       subplot(413); ylabel([Var,' (',Units,')']);
       subplot(414); ylabel([Var,' (',Units,')']);

    Selection = strjoin(Selections);
    if strcmp(Selection,'') == 0; Selection = ['_',Selection];end
     if strcmp(WhichTime,'Gast'); Selection = [Selection,'_Gast'];end
 
        print(Fig1,[ToSave,Selection,'_',Var,'_means.pdf'],'-fillpage', '-dpdf');



    close all
    %end
end