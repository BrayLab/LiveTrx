function[] = Table2Vars(Table)
        for x = [1:length(Table.Properties.VariableNames)]
            Var = Table.Properties.VariableNames{x};
            Par = Table{1,x};
            if iscell(Par) == 1
                command = [Var, ' = ''', Par{:},''''];
            end
            if isnumeric(Par) == 1
                command = [Var, ' = ', num2str(Par)];
            end
             if islogical(Par) == 1
                command = [Var, ' = ', num2str(Par)];
            end
            %eval(['global ',Var]);
            evalin('caller',[command,';']);
        end
end