function[Stats_GFP] = Retrack2Spots(Stats_GFP)

    % retrack spots within each frame, or find a better way to do this
    for f = 2:size(Stats_GFP,1)
        Stats_0 = Stats_GFP{f-1,1};
        Stats_1 = Stats_GFP{f,1};
        % find in f the ones that have 2 spots and correct if assigned
        % wrong. if before only 1, call 1 the closest. if none then just
        % leave the rendom assignment and it will be adjusted int he ones
        % after
        Labels2Spots = Stats_1.Label(~isnan(Stats_1.MaxGauss2));
        for Label = Labels2Spots'
            Index1 = find(Stats_1.Label==Label);
            Index0 = find(Stats_0.Label==Label);
            if isnan(Index0) == 0
                Positions1 = [Stats_1.PositionGaussX(Index1),Stats_1.Position2GaussX(Index1);...
                    Stats_1.PositionGaussY(Index1),Stats_1.Position2GaussY(Index1);...
                    Stats_1.PositionGaussZ(Index1),Stats_1.Position2GaussZ(Index1)];
                Positions0 = [Stats_0.PositionGaussX(Index0),Stats_0.Position2GaussX(Index0);...
                    Stats_0.PositionGaussY(Index0),Stats_0.Position2GaussY(Index0);...
                    Stats_0.PositionGaussZ(Index0),Stats_0.Position2GaussZ(Index0)];
                Distances = [sqrt(sum((Positions0 - Positions1).^2,1)),sqrt(sum((Positions0 - flip(Positions1,2)).^2,1))];
                % order distances: 11, 22, 12, 21. first is in f0, second
                % in f2
                if find(min(Distances)==Distances) > 2
                    Temp = Stats_1.MaxGauss(Index1);
                    Stats_1.MaxGauss(Index1) = Stats_1.MaxGauss2(Index1);
                    Stats_1.MaxGauss2(Index1) = Temp;
                    Stats_1.PositionGaussX(Index1) = Positions1(1,2);
                    Stats_1.Position2GaussX(Index1) = Positions1(1,1);
                    Stats_1.PositionGaussY(Index1) = Positions1(2,2);
                    Stats_1.Position2GaussY(Index1) = Positions1(2,1);
                    Stats_1.PositionGaussZ(Index1) = Positions1(3,2);
                    Stats_1.Position2GaussZ(Index1) = Positions1(3,1);
                end
                
            end
        end
        Stats_GFP{f-1,1} = Stats_0;
        Stats_GFP{f,1} = Stats_1;
    end
end