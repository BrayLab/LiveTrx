function[Out] = FindInMetadata(Info,String)
    Out = NaN;
    try
        Index = find(cellfun(@(x) isempty(x), strfind(Info,String)) == 0);
        D = Info{Index}; D = strsplit(D,'='); Out = str2num(D{2});
    end
end