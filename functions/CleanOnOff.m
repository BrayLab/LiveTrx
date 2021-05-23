function[OnOff] = CleanOnOff(OnOff, MinOn)
OnOffString = join(string(double(OnOff)),'',1);
for m = 1:MinOn
    Beginings = strfind(OnOffString,['0',repmat('1',1,m),'0']);
    for l = 1:size(Beginings,2)
        for n = 0:m+1
         OnOff([Beginings{l}+n],l) = 0;
        end
    end
end

