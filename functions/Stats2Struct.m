function[Data3D] = Stats2Struct(Data3D,Stats_tracked_3D, TimeRes,Frames,minNumb, Path, File, Name,VarSuffix,FileSuffix)

        [AllF] = MergeAll(Stats_tracked_3D, TimeRes);
        AllF = splitvars(AllF);
        %
        Labels = unique(AllF.Label); LabelsOld = Labels;
        
       
        
        VariableNames = AllF.Properties.VariableNames;
        for i = [1:length(VariableNames)-3]
            try
            	Var = VariableNames{i};
                [Variable] = Reshape(AllF,Frames,Labels,Var,'Label');
                if i == 1
                    noNaNnum = sum(~isnan(Variable),1);
                    ToDelete = (noNaNnum <= minNumb);
                    % delete tracks less than minNumb
                end
                Variable(:,ToDelete) = [];
                writetable(array2table(Variable),[Path,File,Name,'tables/',File,'_',FileSuffix,'_',Var,'.txt'],'Delimiter','\t','WriteVariableNames',0)
                Data3D.([Var,VarSuffix]) = Variable;
            end
        end
        LabelsToDelete = Labels(ToDelete);
        Labels(ToDelete) = [];
        writetable(array2table(Labels),[Path,File,Name,'tables/',File,'_',FileSuffix,'_labels.txt'],'Delimiter','\t','WriteVariableNames',0)
          
end
