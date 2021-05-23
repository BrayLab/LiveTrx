function[Labels, Mat] = ReplaceLabels(Labels,LabelsE,LabelsS,Mat)
    toRemove = [];
    %Properties.NewLabel = Properties.Label;
    for x = 1:length(LabelsE);
        IndE = Labels == LabelsE(x);
        IndS = Labels == LabelsS(x);
        if Labels(IndE)==LabelsS(1:x)
            IndE = LabelsE(Labels(IndE) == LabelsS(1:x));
            disp('already replaced!')
        end
        Mat(:,IndE) = nanmean([Mat(:,IndE),Mat(:,IndS)],2);
        toRemove = [toRemove, find(IndS)];
        %Properties.NewLabel(IndS) = Properties.Label(IndE);
    end
    Labels(toRemove) = [];
    Mat(:,toRemove) = [];
    %PropRemoved = Properties(:,toRemove);
    %Properties(:,toRemove) = [];
end