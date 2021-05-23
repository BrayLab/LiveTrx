function[Reshaped] = Reshape(AllF,Frames,Labels,String,By)

    Reshaped = nan(Frames,size(Labels,1));
    for L = 1:size(Labels,1)
        try
            Times = AllF{AllF{:,By} == Labels(L),'Frame'};
            Values = AllF{AllF{:,By} == Labels(L),String};
        catch
            Times = AllF{strcmp(AllF{:,By},Labels{L}),'Frame'};
            Values = AllF{strcmp(AllF{:,By},Labels{L}),String};
        end
        Reshaped(Times,L) = Values;
    end

end
