function[FTL_tracked,FTL_tracked_divisions_RGB,Phylogeny, Families] = PhylogenyDivisions(varargin)
Stats_tracked = varargin{1};
FTL_tracked = varargin{2};
From = varargin{3};
To = varargin{4};
cmap_shuffled = varargin{5};
show = varargin{6};

AllParents = [];
Families = table(NaN,[NaN NaN],NaN, 'VariableNames', {'Parent', 'Daughters','Frame'});
for f=From:To
    Parent = Stats_tracked{f,1}.Parent;
    Daughters = Stats_tracked{f,1}.Daughters;
    Frame = repmat(f,size(Parent,1),1);
    Family = table (Parent, Daughters, Frame);
    Families=vertcat(Families,Family);
end
AllParents = Families.Parent;
AllParents(isnan(AllParents)|AllParents==0)=[];
AllParents = unique(AllParents');
Families(1,:) = [];
%%
Phylogeny = cell(0,1);
for P=AllParents
    %disp(P)
    Offspring = Families(Families.Parent==P,:).Daughters;
    Offspring(Offspring==0 | isnan(Offspring))=[];
    Offspring=unique(Offspring);
    %if parent is already a daughter, add to previous phylogeny
    if isempty(find([Phylogeny{:}]==P))
            Phylogeny{P} = Offspring;
    else
        %find its parent
        Table =table(Phylogeny(:));
        dummy = cellfun(@(X) sum(ismember(X,P))==1,Phylogeny, 'UniformOutput',false);
        Grandmother = find([dummy{:}]==1);
          Grandmother = Grandmother(1);
%         % add offpring to previous one
        Phylogeny{Grandmother} = unique([Phylogeny{Grandmother},Offspring]);
    end            
end
%%
Daughters = Families{:,2};
Families(isnan(Daughters(:,2)),:) = [];
Daughters = Families{:,2};
Families(Daughters(:,2)==0,:) = [];
Families(Families.Parent==0,:) = [];
%%
OriginalParents = table([find(cellfun(@isempty,Phylogeny)==0)]','VariableNames',{'OriginalParent'});
OriginalParents.NewIndex = [1:size(OriginalParents,1)]';
Phylogeny(cellfun(@isempty,Phylogeny)==1)=[];
%% replace labels in each phylogeny
%FTL_tracked_divisions = zeros(size(FTL_tracked));
for P=OriginalParents.NewIndex'
    %disp(P)
    L=OriginalParents.OriginalParent(P,1);
    toReplace = Phylogeny{1,P};
    FTL_tracked(ismember(FTL_tracked,toReplace))=L;
end
Dim = length(size(FTL_tracked));
for f=1:size(FTL_tracked,Dim)
    if Dim == 3
        FTL_tracked_divisions_RGB(:,:,:,f) = label2rgb(FTL_tracked(:,:,f), cmap_shuffled, 'k', 'noshuffle');
    end
    if Dim == 4
        FTL_tracked_divisions_RGB(:,:,:,f) = label2rgb(MAX_proj_3D(FTL_tracked(:,:,:,f)), cmap_shuffled, 'k', 'noshuffle');
    end
end
if strcmp(show,'on')==1;
    mov = immovie(FTL_tracked_divisions_RGB);
            implay(mov)
end
end