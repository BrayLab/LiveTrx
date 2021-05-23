function[Properties] = ReplaceLabelsProps(Labels,LabelsE,LabelsS,Properties)
    %[LabelsE,LabelsS]
    %toRemove = [];
    Properties.NewLabel = num2cell(Properties.Label);
    for x = 1:length(LabelsE);
        disp(x)
        EndL = LabelsE(x);
        IndE = Labels == EndL;
        counter = 0;
        while isempty(Labels(IndE))
            EndL = LabelsE(EndL == LabelsS);
            IndE = Labels == EndL;
            counter = counter+1;
            if counter>100
                break
            end
        end  
        if counter < 100
            Properties.NewLabel{IndE} = [Properties.NewLabel{IndE},LabelsS(x)];
        end
    end
end