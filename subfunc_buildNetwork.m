function [ cellsPosition ] = subfunc_buildNetwork(L,rc)
if nargin < 1
    L = 48;
end
if nargin < 2
    rc = 800;
end
% Build netwrok ----------------------------------------------------------
        numOfLevels = 8;
        numOfCellsInEachLevel = [4 5 6 7 8 7 6 5];
        cellsIndex            = [ 1  2  5  6, ....
                                  4  3  8  7 24,...
                                  9 10 13 14 17 18,...
                                 12 11 16 15 20 19 40,...
                                 21 22 25 26 29 30 33 34,...
                                 23 28 27 32 31 36 35,...
                                 37 38 41 42 45 46,...
                                 39 44 43 48 47];
        % build y position of each level
        cellsPosition = zeros(L,2);
        Ind    = 0;
        xStart = sqrt(3/2)*rc;
        yStart = 0;
        for l = 1 : numOfLevels
            range = Ind + (1:numOfCellsInEachLevel(l));
            Ind   = Ind + numOfCellsInEachLevel(l);
            % build x positions
            if l <= 5
                xStart = xStart - sqrt(3)*rc/2;
            else
                xStart = xStart + sqrt(3)*rc/2;
            end
            cellsPosition(cellsIndex(range),1) = xStart + (0:(numOfCellsInEachLevel(l)-1))*sqrt(3)*rc;
            % build y positions
            cellsPosition(cellsIndex(range),2) = yStart;
            yStart = yStart - (3/2)*rc;
        end 

    if L == 48
        % do nothing
    elseif L <= 7
        cellsPosition = cellsPosition([3 4 1 2 8 13 10],:);        
        cellsPosition = cellsPosition(1:L,:);
    else
        error('L = %d is not supported.',L)
    end
%-------------------------------------------------------------------------- 
end