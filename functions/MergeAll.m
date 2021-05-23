function[AllF] = MergeAll(Stats_GFP,TimeRes)
AllF = Stats_GFP{1,1};
try
    AllF.PixelIdxList = [];
    AllF.MaxIntensity = [];
    AllF.SpotPositionX = cellfun(@(x) x(1),Stats_GFP{1}.SpotPosition);
     AllF.SpotPositionY = cellfun(@(x) x(2),Stats_GFP{1}.SpotPosition);
     AllF.SpotPositionZ = cellfun(@(x) x(3),Stats_GFP{1}.SpotPosition);
end
AllF.Frame = repmat(1,size(AllF,1),1);
AllF.Time = repmat(1*TimeRes/60,size(AllF,1),1);
for f = 2:size(Stats_GFP)
    Table = Stats_GFP{f,1};
    try
        Table.PixelIdxList = [];
        Table.MaxIntensity = [];
        Table.SpotPositionX = cellfun(@(x) x(1),Stats_GFP{f}.SpotPosition);
        Table.SpotPositionY = cellfun(@(x) x(2),Stats_GFP{f}.SpotPosition);
        Table.SpotPositionZ = cellfun(@(x) x(3),Stats_GFP{f}.SpotPosition);
    end
    Table.Frame = repmat(f,size(Table,1),1);
    Table.Time = repmat(f*TimeRes/60,size(Table,1),1);
    try
        AllF = vertcat(AllF, Table);
    end
end
% AllF = table();
% try
% AllF.CentroidX = Stats_GFP{1}.Centroid(:,1);
% AllF.CentroidY = Stats_GFP{1}.Centroid(:,2);
% AllF.Label = Stats_GFP{1}.Label;
% AllF.Frame = repmat(1,size(AllF,1),1);
% AllF.Time = repmat(1*TimeRes/60,size(AllF,1),1);
%     try
%     AllF.Max = Stats_GFP{1}.Max;
%     AllF.CentroidZ = Stats_GFP{1}.Centroid(:,3);
%     AllF.SpotX = cellfun(@(x) x(1),Stats_GFP{1}.SpotPosition);
%     AllF.SpotY = cellfun(@(x) x(2),Stats_GFP{1}.SpotPosition);
%     AllF.SpotZ = cellfun(@(x) x(3),Stats_GFP{1}.SpotPosition);
%     AllF.NewLabel = Stats_GFP{1}.NewLabel;
%     end
% 
% end
% % same with spot position
% for f = 2:size(Stats_GFP)
%     try
%     Table = table();
%     Table.CentroidX = Stats_GFP{f}.Centroid(:,1);
%     Table.CentroidY = Stats_GFP{f}.Centroid(:,2);
%     Table.Label = Stats_GFP{f}.Label;
%     Table.Frame = repmat(f,size(Table,1),1);
%     Table.Time = repmat(f*TimeRes/60,size(Table,1),1);
%         try
%         Table.Max = Stats_GFP{f}.Max;
%         Table.CentroidZ = Stats_GFP{f}.Centroid(:,3);
%         Table.SpotX = cellfun(@(x) x(1),Stats_GFP{f}.SpotPosition);
%         Table.SpotY = cellfun(@(x) x(2),Stats_GFP{f}.SpotPosition);
%         Table.SpotZ = cellfun(@(x) x(3),Stats_GFP{f}.SpotPosition);
%         Table.NewLabel = Stats_GFP{f}.NewLabel;
%         end
%     AllF = vertcat(AllF, Table);
%     end
% end
end