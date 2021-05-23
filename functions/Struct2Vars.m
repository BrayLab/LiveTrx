function[] = Struct2Vars(Data)
    %Struct = Data.Data;
    for x = 1:length(fieldnames(Data))
        Vars = fieldnames(Data);
        %NewField = Struct(Vars{x});
        %command = [Vars{x}, ' = Data.Data.',Vars{x},''];
        command = [Vars{x}, ' = Data.',Vars{x},';'];
        %evalc(['global ',Vars{x}]);
        evalin('caller',command);
    end
end