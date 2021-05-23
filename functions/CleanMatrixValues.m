function[CentX] = CleanMatrixValues(AllF,Frames,LabelsOld,String,ToDelete,BG)
            [CentX] = Reshape(AllF,Frames,LabelsOld,String,'Label');
            CentX(:,ToDelete) = [];
            CentX = CentX(:,~BG);
    end