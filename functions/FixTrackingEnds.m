function[LabelsE,LabelsS] = FixTrackingEnds(MaxF,Labels,TimeRes,CentX,CentY,CentZ,minPDis,nc14)
%% fix tracking
NaNString = join(string(double(~isnan(MaxF(nc14:end,:)))),'',1);
StartF = nan(1,length(Labels));
Beginings = strfind(NaNString,'01');
GoodOnes = cellfun(@(x) ~isempty(x), Beginings);
StartL = Labels(GoodOnes);
Dummy = cellfun(@(x) x(1)+1, Beginings(GoodOnes));
StartF(GoodOnes) = Dummy;
EndF = nan(1,length(Labels));
Ends = strfind(NaNString,'10');
GoodOnes = cellfun(@(x) ~isempty(x), Ends);
EndL = Labels(GoodOnes);
Dummy = cellfun(@(x) x(end), Ends(GoodOnes));
EndF(GoodOnes) = Dummy;
%%
StartT = StartF(~isnan(StartF)).* TimeRes./60; % in min
IndicesS = sub2ind(size(MaxF),StartF(~isnan(StartF)),find(~isnan(StartF)));
StartX = CentX(IndicesS); %um
StartY = CentY(IndicesS); %um
StartZ = CentZ(IndicesS); %um
EndT = EndF(~isnan(EndF)).* TimeRes./60;
IndicesE = sub2ind(size(MaxF),EndF(~isnan(EndF)),find(~isnan(EndF)));
EndX = CentX(IndicesE);
EndY = CentY(IndicesE);
EndZ = CentZ(IndicesE);
% calculate eucleudian distance with X, Y, Z and T, min distances in 4D
% all starts with all ends
%%
Eucledian = sqrt((StartT-EndT').^2 + (StartX-EndX').^2 +(StartY-EndY').^2 +(StartZ-EndZ').^2);
OEuc = Eucledian;
%imagesc(Euc); colormap gray
Rows = []; Cols = []; Mins = [];
Min = 0;
% get row and col labels from min eucledian distance and remove row and
% column succesively until all are gone
while size(Eucledian) >1 & Min < minPDis
    [Min, Index] = min(Eucledian(:));
    [R,C] = find(OEuc==Min,1);
    Rows = [Rows, R]; Cols = [Cols, C]; Mins = [Mins, Min];
    [R,C] = ind2sub(size(Eucledian),Index);
    Eucledian(R,:) = []; 
    Eucledian(:,C) = [];
end
Indices = sub2ind(size(OEuc),Rows, Cols);
%hist(OEuc(Indices))
%OEuc(Indices) = 10000; 
%imagesc(OEuc)
LabelsS = StartL(Cols); 
LabelsE = EndL(Rows);
TimesS = StartT(Cols); 
TimesE = EndT(Rows);
Eq = LabelsS == LabelsE;
LabelsS(Eq) = [];
LabelsE(Eq) = [];
TimesS(Eq) = [];
TimesE(Eq) = [];
%[LabelsE,LabelsS,TimesE',TimesS']
end