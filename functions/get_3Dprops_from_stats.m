function[]=get_3Dprops_from_stats(Parameters)

%Parameters = info(x,:);
Table2Vars(Parameters);
Flip = str2double(strsplit(Flip,','));
PathToSave = [Path,File,Name,File]; 
reader = bfGetReader([Path,File]);
[Bits,Width,Height, Channels, Slices, Frames, Frames0, XRes, YRes, ZRes,Zoom,TimeRes,Settings] = readMetadataBFOME(reader);
try
    XYRes=mean(round(mean([XRes,YRes]),2));
    ZRes = mean(round(ZRes,2));
end
if strcmp(To,'NA')==1 | To == inf;  To=Frames-1;
elseif isnumeric(To) == 0; To = str2num(To); end;
Frames = To-From+1; 
if strcmp(Channel1,'nucmemb') | strcmp(Channel2,'nucmemb')
    NucMemb = 1;
else
    NucMemb = 0;
end

%% rescue toTrack from Stats3D and delete discarded objects
% and GET 3D PROPERTIES OF NUCLEI AND 2D PROPERTIES IN SLICES

%toTrack = zeros(Height,Width, round(Slices*ZRes/XYRes), Frames);
load([Path,File,Name,File,'_Stats.mat']);

% clean objects tracked for less than minimum after tracking and deleting 
% them from toTrack (CHANGE SO THAT IT DOES IT AFTER FIX TRACKING)

minNumb = 10;

[AllF] = MergeAll(Stats_GFP, TimeRes);
AllF = splitvars(AllF);
        %
Labels = unique(AllF.Label); LabelsOld = Labels;
[Area] = Reshape(AllF,Frames,Labels,'Area','Label');
       
noNaNnum = sum(~isnan(Area),1);
ToDelete = (noNaNnum <= minNumb);
% delete tracks less than minNumb
LabelsToDelete = Labels(ToDelete);
Labels(ToDelete) = [];
Area(:,ToDelete) = [];


Stats_tracked_3D = cell(Frames,1);
Stats_tracked_2D_Q1 = cell(Frames,1);
Stats_tracked_2D_Q2 = cell(Frames,1);
Stats_tracked_2D_Q3 = cell(Frames,1);
PropObjects = zeros(Frames,round(Slices*ZRes/XYRes));

Fig = figure;
set(Fig,'Visible', 'on','PaperOrientation','landscape');
cmapgrad = parula(Frames*4);
set(gcf,'defaultAxesColorOrder',cmapgrad)

cmap = jet(1000000);
cmap_shuffled = cmap(randperm(size(cmap,1)),:);
Boundaries_RGB_clean = zeros(Height,Width,3,Frames);


for f = 1:length(Stats_GFP)
    disp(num2str(f))
    SingleFrame = Stats_GFP{f};
    SingleFrame(isnan(SingleFrame.Label),:) = [];
    if NucMemb
        Tracked3D = zeros(Height,Width,round(Slices*ZRes/XYRes));
    else
        Tracked3D = zeros(Height,Width,Slices);
    end
    for n = 1:size(SingleFrame,1)
        Tracked3D(SingleFrame.PixelIdxList{n}) =  SingleFrame.Label(n);
    end
    PixelsToRemove = ismember(Tracked3D,LabelsToDelete);
    Tracked3D(PixelsToRemove) = 0;
    if ~ NucMemb
        Tracked3D = imresize3(Tracked3D,[Width,Height,round(Slices*ZRes/XYRes)],'nearest');
    end
%     if NucMemb
%         toTrack(:,:,:,f) = Tracked3D;
%     else
%         toTrack(:,:,:,f) = imresize3(Tracked3D,[Width,Height,round(Slices*ZRes/XYRes)],'nearest');
%     end
%end        
      %save second clean
    [~, Boundaries_RGB_clean(:,:,:,f), ~] = BoundariesTracked_3D(Tracked3D,cmap_shuffled,'off');
%imshow(Boundaries_RGB_clean(:,:,:,f))


%for f = 1:Frames
        disp(['running regionprops f',num2str(f),'...']);
        %TL = toTrack(:,:,:,f);
        TL = Tracked3D;
        Stats3 = regionprops3(TL,'BoundingBox','Centroid','ConvexVolume',...
           'Orientation','PrincipalAxisLength','Solidity','SurfaceArea',...
           'Volume','VoxelIdxList');
        
        Stats3.Label = [1:length(Stats3.Volume)]';
        ToRemove = find(Stats3.Volume == 0);
        Stats3(ToRemove,:) = [];
        Stats_tracked_3D{f,1} = Stats3;
        
        %P = permute(TL==0,[3,2,1]);
        %P = sum(P,3); P = 1-sum(P,2)./max(sum(P,2));
        %plot(P); hold on
        %drawnow
        %[~,CenterSlice] = max(P); % center slice is where biggest proportion of objects in the frame
        for s = 1:size(TL,3)
            P(s) = sum(unique(TL(:,:,s))~=0);
        end
        
        if sum(P~=0) > 0
            P = P./max(P);
            plot(P); hold on
            Indices = find(P>0.5);
            PropObjects(f,:) = P;
            %CenterSlice = Indices(round(length(Indices)/2)); % center slice is the slice in the center of the top 50% of unique labels
            Slice25(f) = Indices(round(length(Indices)*0.25)); % center slice is the slice in the center of the top 50% of unique labels
            Slice50(f) = Indices(round(length(Indices)*0.5)); % center slice is the slice in the center of the top 50% of unique labels
            Slice75(f) = Indices(round(length(Indices)*0.75)); % center slice is the slice in the center of the top 50% of unique labels

            plot(Slice25(f),0.25,'.')
            plot(Slice50(f),0.5,'.')
            plot(Slice75(f),0.75,'.')

            drawnow

            Stats25 = regionprops('table',TL(:,:,Slice25(f)),'Area','BoundingBox',...
                'Centroid','ConvexArea','Eccentricity','EquivDiameter',...
                'Extent','Perimeter','Solidity');
            Stats25.Label = [1:length(Stats25.Area)]';
            ToRemove = find(Stats25.Area == 0);
            Stats25(ToRemove,:) = [];
            Stats_tracked_2D_Q1{f,1} = Stats25;

            Stats50 = regionprops('table',TL(:,:,Slice50(f)),'Area','BoundingBox',...
                'Centroid','ConvexArea','Eccentricity','EquivDiameter',...
                'Extent','Perimeter','Solidity');
            Stats50.Label = [1:length(Stats50.Area)]';
            ToRemove = find(Stats50.Area == 0);
            Stats50(ToRemove,:) = [];
            Stats_tracked_2D_Q2{f,1} = Stats50;

            Stats75 = regionprops('table',TL(:,:,Slice75(f)),'Area','BoundingBox',...
                'Centroid','ConvexArea','Eccentricity','EquivDiameter',...
                'Extent','Perimeter','Solidity');
            Stats75.Label = [1:length(Stats75.Area)]';
            ToRemove = find(Stats75.Area == 0);
            Stats75(ToRemove,:) = [];
            Stats_tracked_2D_Q3{f,1} = Stats75;
        else
            Slice25(f) = NaN;
            Slice50(f) = NaN;
            Slice75(f) = NaN;
            PropObjects(f,:) = NaN;
        end
end

save([Path,File,Name,File,'_Stats3D.mat'],'Stats_tracked_3D','-v7.3');
save([Path,File,Name,File,'_Stats2D_Q1.mat'],'Stats_tracked_2D_Q1','-v7.3');
save([Path,File,Name,File,'_Stats2D_Q2.mat'],'Stats_tracked_2D_Q2','-v7.3');
save([Path,File,Name,File,'_Stats2D_Q3.mat'],'Stats_tracked_2D_Q3','-v7.3');

print(Fig,[PathToSave,'_Slices.pdf'],'-fillpage', '-dpdf');

TableQ = array2table([(1:length(Stats_GFP))',Slice25', Slice50', Slice75'],'VariableNames',{'Frame','SliceQ1','SliceQ2','SliceQ3'});
writetable(TableQ,[PathToSave,'_Slices.txt'],'Delimiter','\t')
%writetable(PropObjects,[PathToSave,'_PropObjects.txt'],'Delimiter','\t')
Data3D.PropObjects = PropObjects;
Data3D.QuartileSlices = TableQ;

Data3D = NuclearMembraneProperties(0, Path, File, Name, Stats_tracked_3D, TimeRes, Frames, nc14, Delay, Data3D,minNumb);
save([Path,File,Name,File,'_Data3D.mat'],'Data3D','-v7.3');
clear Stats_GFP
clear Data3D
%clear toTrack

% SAVE SECOND CLEAN??
        WriteRGB(Boundaries_RGB_clean, PathToSave, '_segmented_tracked_boundaries_RGB_clean2.tiff','none')

end

