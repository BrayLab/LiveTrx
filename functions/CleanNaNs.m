function[OnOff] = CleanNaNs(MedFilt,OnOff, minOn)
OnOff = double(OnOff);
OnOff(isnan(MedFilt)) = NaN;
OnOffString = string(OnOff);
OnOffString = fillmissing(OnOffString,'constant',"N");
OnOffString = join(OnOffString,'',1);

for m = 1:minOn
    Beginings = strfind(OnOffString,['0',repmat('N',1,m),'0']);
    for l = 1:size(Beginings,2)
        for n = 0:m+1
         OnOff([Beginings{l}+n],l) = 0;
        end
    end
end

for m = 1:minOn
    Beginings = strfind(OnOffString,['1',repmat('N',1,m),'1']);
    for l = 1:size(Beginings,2)
        for n = 0:m+1
         OnOff([Beginings{l}+n],l) = 1;
        end
    end
end

OnOff(isnan(OnOff)) = 0;
end