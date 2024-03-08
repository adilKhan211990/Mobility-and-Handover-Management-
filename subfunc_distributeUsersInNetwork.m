function [usersPositions,usersDist2Cells] = subfunc_distributeUsersInNetwork(L,rc,rh,ru,cellsPosition,usersCellInd)
if nargin < 1
    L = 48;
end
if nargin < 2
    rc = 800;
end
if nargin < 3
    rh = 100;
end
if nargin < 4
    cellsPosition = zeros(L,2);
end

    usersNum       = size(usersCellInd(:),1);
    usersPositions = zeros(usersNum,2);
    K0             = 0; % user counter
    for l = 1 : L
        % generate user positions relative to their cells position
        % sqrt(3)/2 can be multiplied to keep all the nodes inside the cell
        % uniform positions
        K       = nnz(usersCellInd == l);
        randPos = zeros(K,2);
        for kk = 1 : K
            Continue = true; % it is a flag
            while Continue
                Pos = rc*(2*rand(1,2)-1);
                % check if the user position is with the radius rh
                % and rh
                if norm(Pos) > rh && norm(Pos) < rc*ru
                    Continue      = false;   % stp search
                    randPos(kk,:) = Pos; % choose this position
                end
            end
        end
        % shift relative users position to their absolute value by adding
        % the cell positions
        usersPositions(K0+(1:K),:) = [randPos(:,1)+cellsPosition(l,1),...
                                      randPos(:,2)+cellsPosition(l,2)];
        % update user counter
        K0 = K0 + K;
    end
%-------------------------------------------------------------------------- 
    % Memory allocation
    usersDist2Cells = zeros(usersNum,L);
    for l = 1 : L
        usersDist2Cells(:,l) = sqrt(...
        (usersPositions(:,1) - cellsPosition(l,1)).^2 + ...
        (usersPositions(:,2) - cellsPosition(l,2)).^2);    
    end
    usersDist2Users = zeros(usersNum,usersNum);
    for kk = 1 : usersNum
        usersDist2Users(:,kk) = sqrt(...
        (usersPositions(:,1) - usersPositions(kk,1)).^2 + ...
        (usersPositions(:,2) - usersPositions(kk,2)).^2);    
    end
    cellsDist2Cells = zeros(L);
    for l = 1 : L
        cellsDist2Cells(:,l) = sqrt(...
        (cellsPosition(:,1) - cellsPosition(l,1)).^2 + ...
        (cellsPosition(:,2) - cellsPosition(l,2)).^2); 
    end
end