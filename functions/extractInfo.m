function[D] = extractInfo(Info, String)
Index = find(cellfun(@(x) isempty(x), regexp(Info,String,'match')) == 0);
D = Info(Index); D = cellfun(@(x) strsplit(x,'='),D,'UniformOutput',false);
D = cell2mat(cellfun(@(x) str2double(strrep((x{2}), '"', '')),D,'UniformOutput',false));
end