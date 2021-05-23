function[Out] = Cell2Mat(BNumAll)
    MaxLength = max(cellfun(@(x) length(x),BNumAll));
    Out = nan(MaxLength, size(BNumAll,2));
    for i = 1:size(BNumAll,2)
        Out(1:length(BNumAll{i}),i) = [BNumAll{i}];
    end


end
