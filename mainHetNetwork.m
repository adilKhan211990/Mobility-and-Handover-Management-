%% Clear ------------------------------------------------------------------
%if nargin == 0
    clear                     % clear all workspace variables
   %close all                 % close all open figures
    rng(1)                    % reset random genetator seed
    clc                       % clear command window    
%end
%-------------------------------------------------------------------------- 
%% Parameters -------------------------------------------------------------
    % control parameters
    Ctrl.ifShow_Figure1_NetSetup = true;
    Ctrl.ifShow_Figure1_Update   = true;

    % network parameters
    M   = 2;             % number of macro-cells (up to 7)
    rcM = 15e2;          % macro radius
    rhM = 1e1;           % Minimum distance of a user from a macro (Meter)
    ruM = sqrt(3)/2;     % ratio maximum user distance from the center of macro
    fM = 900e6;          % operating frequency of the macro
    eM = 2.1;            % path-loss exponent of macro in dB
    pM = 36;             % reference attenuation
    PM = 0;              % transmit power in dB
    
    L   = 7;             % number of BS (always 7)
    rcL = 5e2;           % BS radius
    rhL = 1e1;           % Minimum distance of a user from a BS (Meter)
    ruL = sqrt(3)/2;     % ratio maximum user distance from the center of BS
    fL = 2.4e6;          % operating frequency of the BS
    eL = 3.0;            % path-loss exponent of BS in dB
    pL = 22;             % reference attenuation
    PL = 0;              % transmit power in dB
    
    F  = 5;              % number of femto cells in each BS
    rcF = 1e2;           % femto radius
    rhF = 1e1;           % Minimum distance of a user from a femto (Meter)
    ruF = sqrt(3)/2;     % ratio maximum user distance from the center of femto
    fF = 18e9;           % operating frequency of the femto
    eF = 4.0;            % path-loss exponent of femto in dB
    pF = 20;             % reference attenuation
    PF = 0;              % transmit power in dB

    
    % related to users
    usersMaxNumber       = 20*L*M;   % total number of users in the network
    disp (usersMaxNumber);
    % related to handover
    usersCandidate             = 1;      % candidate user which moves in a network
    usersCandidateVelocity     =  [20 5]*1;  % m/s x and y 
    usersCandidateAcceleration =  [-1  0];  % m/s x and y 
    ksi                        = 0.5;
    sigma_v                    = 1;

    macroMaxUser  = 30;  % maximum number of user that can be handled by each macro
    bsMaxUser     = 20;  % maximum number of user that can be handled by each bs
    femtoMaxUser  =  5;  % maximum number of user that can be handled by each femto
    % simulation parameters
    simEndTime    = 100; % one second
    simStepTime   = 5;
    
%--------------------------------------------------------------------------
%% dependent parameters ---------------------------------------------------
    % create macro-cell positions
    [ macroCellsPosition ] = subfunc_buildNetwork(M,rcM);
    
    % create BS positions
    bsPositions = zeros(L,2,M);
    for m = 1 : M
        [ bsPositions((1:L),:,m) ] = subfunc_buildNetwork(L,rcL);
          bsPositions((1:L),1,m)   = bsPositions((1:L),1,m) - bsPositions(1,1,m) + macroCellsPosition(m,1);
          bsPositions((1:L),2,m)   = bsPositions((1:L),2,m) - bsPositions(1,2,m) + macroCellsPosition(m,2);          
    end
    
    % create femto-cell positions at random
    femtocellPositions = zeros(F,2,L,M);
    for m = 1 : M
        for l = 1 : L
            for f = 1 : F
                % raise a flag
                ifDistanceIsOk = false;
                while ifDistanceIsOk == false
                    % generate random position in cell l of macro-cell m
                    femtocellPositions(f,:,l,m) = (rcL-2*rcF)*(2*rand(1,2)-1);
                    if f == 1 % check if the distance from other femto-cell is ok
                        ifDistanceIsOk = true;
                    else
                        % calulate the distance from previous femto-cells
                        dist = sqrt( ...
                            (femtocellPositions(f,1,l,m)-femtocellPositions(1:f-1,1,l,m)).^2 + ...
                            (femtocellPositions(f,2,l,m)-femtocellPositions(1:f-1,2,l,m)).^2 );
                        if all(dist > 2*rcF)
                            ifDistanceIsOk = true;
                        end
                    end
                end
            end
            femtocellPositions(:,1,l,m) = femtocellPositions(:,1,l,m) + bsPositions(l,1,m);
            femtocellPositions(:,2,l,m) = femtocellPositions(:,2,l,m) + bsPositions(l,2,m);
        end
    end
   
    % allocate users in cells
    usersCellInd  = randi(L*M,usersMaxNumber,1);    % random allocation of users in bs-cells

    % find users positions
    bs2DPositions = zeros(L*M,2);
    for m = 1 : M
        bs2DPositions((m-1)*L + (1:L),:) = bsPositions(:,:,m);
    end
    usersPositions = ...
        subfunc_distributeUsersInNetwork(L*M,rcL,rhL,ruL,bs2DPositions,usersCellInd);
    
    % candidate user position
    usersCandidatePosition = usersPositions(usersCandidate,:);
    
    % total number of cells in the network for global cell allocation 
    nCells = M + M*L + M*L*F;
    % collect all cell positions
    cellsPosition = zeros(nCells,2);
    cellsId       = zeros(nCells,1);    % 1 for macro, 2 for BS, 3 for femto
    kCell = 0;      % cell counter
    for m = 1 : M                               % bs
        kCell = kCell + 1;
        cellsPosition(kCell,:) = macroCellsPosition(m,:);  % macro cells
        cellsId(kCell) = 1;
    end
    for m = 1 : M                               % bs
        for l = 1 : L
            kCell = kCell + 1;
            cellsPosition(kCell,:) = bsPositions(l,:,m);
            cellsId(kCell) = 2;
        end
    end
    for m = 1 : M                               % femto cells
        for l = 1 : L
            for f = 1 : F
                kCell = kCell + 1;
                cellsPosition(kCell,:) = femtocellPositions(f,:,l,m);
                cellsId(kCell) = 3;
            end
        end
    end
    
    % find users distance from each cell (macro,bs,femto)
    usersDist2Cells = zeros(usersMaxNumber,nCells);
    for k = 1 : usersMaxNumber
        for c = 1 : nCells
            usersDist2Cells(k,c) = norm(usersPositions(k,:)-cellsPosition(c,:));
        end
    end
    % find user power from each cell
    usersRxPow = zeros(usersMaxNumber,1);
    e  = [eM,eL,eF];                         % path-loss exponent
    rd = [rhM,rhL,rhF];                      % reference distance
    p  = [pM,pL,pF];                         % reference power
    P  = [PM,PL,PF];                         % transmit power
    for k = 1 : usersMaxNumber
        for c = 1 : nCells
            d = max(rd(cellsId(c)),usersDist2Cells(k,c));
            usersRxPow(k,c) = P(cellsId(c))-p(cellsId(c))-10*e(cellsId(c))*log10(d/rd(cellsId(c)));
        end
    end
    
    % allocate users to femtos, BS or macros based on power
    [~,usersSelectedCells] = max(usersRxPow,[],2);
    % find the load of each femto,BS, or macro
    cellsMaxUser = [macroMaxUser,bsMaxUser,femtoMaxUser];
    cellsHist = hist(usersSelectedCells,1:nCells);
    cellsLoad = cellsHist./cellsMaxUser(cellsId);
%--------------------------------------------------------------------------
%% plot network setup -----------------------------------------------------
if Ctrl.ifShow_Figure1_NetSetup
    f1 = figure(1);
    f1.CloseRequestFcn = '';
        clf

        % plot macrocells
        boundaries = rcM*[sin(0:pi/100:2*pi).',...
                          cos(0:pi/100:2*pi).'];
        kCell = 0;
        for m = 1 : M
                plot(macroCellsPosition(m,1)/1e3 + boundaries(:,1)/1e3,...
                     macroCellsPosition(m,2)/1e3 + boundaries(:,2)/1e3,'m:');    
                hold on 
                kCell = kCell + 1;
                text(      macroCellsPosition(m,1)/1e3 ,...
                       max(macroCellsPosition(m,2)/1e3 + boundaries(:,2)/1e3),sprintf('%d',kCell))
        end
        
        % plot BSs
        boundaries = rcL*[sin(0:pi/3:2*pi).',...
                          cos(0:pi/3:2*pi).'];
        Colors = [0 1 0];
        for m = 1 : M
            for l = 1 : L
                fill(bsPositions(l,1,m)/1e3 + boundaries(:,1)/1e3,...
                     bsPositions(l,2,m)/1e3 + boundaries(:,2)/1e3,...
                        Colors,'FaceAlpha',0.1);                
                plot(bsPositions(l,1,m)/1e3 + boundaries(:,1)/1e3,...
                     bsPositions(l,2,m)/1e3 + boundaries(:,2)/1e3,'k');   
                 
                kCell = kCell + 1;
                text(      bsPositions(l,1,m)/1e3 ,...
                       max(bsPositions(l,2,m)/1e3 + boundaries(:,2)/1e3),sprintf('%d',kCell))
                   
                hold on
            end
        end
        % plot femto-cells
        boundaries = rcF*[sin(0:pi/100:2*pi).',...
                          cos(0:pi/100:2*pi).'];
        Colors = [0 0 1];
        for m = 1 : M
            for l = 1 : L
                for f = 1 : F
                    fill(femtocellPositions(f,1,l,m)/1e3 + boundaries(:,1)/1e3,...
                         femtocellPositions(f,2,l,m)/1e3 + boundaries(:,2)/1e3,...
                            Colors,'FaceAlpha',0.1);                
                    plot(femtocellPositions(f,1,l,m)/1e3 + boundaries(:,1)/1e3,...
                         femtocellPositions(f,2,l,m)/1e3 + boundaries(:,2)/1e3,'k');    
                    kCell = kCell + 1;
                    text(      femtocellPositions(f,1,l,m)/1e3 ,...
                           max(femtocellPositions(f,2,l,m)/1e3 + boundaries(:,2)/1e3),sprintf('%d',kCell))
                end
            end
        end
        
        % plot users
        hold on
        plot(usersPositions(:,1)/1e3,usersPositions(:,2)/1e3,'r.')
        
        % plot candidate users
        hCandidate = plot(usersCandidatePosition(1)/1e3,...
                          usersCandidatePosition(2)/1e3,'om');
        
        
        axis equal
        drawnow
        xlabel('x-axis [Km]')
        ylabel('y-axis [Km]')      
        title('network setup (Macro,BS,femto and users')
        
    figure(2)
        clf
        colors = 'brg'; % color for macro, bs, and femto
        for n = 0 : 3
            subplot(2,1,1)
                if n == 0
                    bar(cellsHist)
                else
                    Ind = find(cellsId == n); % find the indices of the cell
                	h1(n) = bar(Ind,cellsHist(Ind),colors(n));%#ok
                end
                hold on
                xlabel('cell ID')
                ylabel('Num of user in cell')
            subplot(2,1,2)            
                if n == 0
                    bar(cellsLoad)
                else
                    h2(n) = bar(Ind,cellsLoad(Ind),colors(n));%#ok
                end
                xlabel('cell ID')
                ylabel('cells Load')
                hold on
        end
        legend(h1,'Macro BS','Micro','Femto')
        legend(h2,'Macro BS','Micro','Femto')
end

%--------------------------------------------------------------------------
%% Simulation -------------------------------------------------------------
    simNumOfSteps = ceil(simEndTime./simStepTime);
    
    % extract information
    simT     = zeros(simNumOfSteps,1);
    simPos   = zeros(simNumOfSteps,2);
    simAvgV  = zeros(simNumOfSteps,2);
    simDir   = zeros(simNumOfSteps,2);
    simPow   = zeros(simNumOfSteps,nCells);
    simDist  = zeros(simNumOfSteps,nCells);
    simAngle = zeros(simNumOfSteps,nCells);
    simEstP  = zeros(simNumOfSteps,nCells)*NaN;
    simEx   = zeros(simNumOfSteps,nCells)*NaN;
    
    % throuput for the proposed and LET system
    simThrPro  = zeros(simNumOfSteps,1)*NaN;    
    simThrEx  = zeros(simNumOfSteps,1)*NaN;
    % load of the selected cell
    simLoadPro  = zeros(simNumOfSteps,1)*NaN;    
    simLoadEx  = zeros(simNumOfSteps,1)*NaN;
    
    % decision factor 
    usersCandidateAvgVelocity = norm(usersCandidateVelocity);
    
    % update load of the cell which candidate user belongs to 
    candidateCellIndex = usersSelectedCells(usersCandidate);
    candidateCellType  = cellsId(candidateCellIndex);
    % remove candidate user from this cell and update its load, and let the
    % user select also this cell in case it interested.
    cellsLoad(candidateCellIndex) = ...
        (cellsLoad(candidateCellIndex)*cellsMaxUser(candidateCellType)-1)/...
        cellsMaxUser(candidateCellType);
    cellsHist(candidateCellIndex) = cellsHist(candidateCellIndex) - 1;
    
    simTime = 0;
    if Ctrl.ifShow_Figure1_Update  == true
        figure(1)
    end
    for simCnt = 1 : simNumOfSteps
        
        % update user position
        usersCandidatePosition = usersCandidatePosition + ...
            usersCandidateVelocity.*simStepTime + ... 
            usersCandidateAcceleration.*simStepTime.^2;
        % update user velocity 
        usersCandidateVelocity = usersCandidateVelocity + ...
            usersCandidateAcceleration.*simStepTime;
        % update moving direction 
        usersCandidateDirection = atan2(usersCandidateVelocity(2),...
                                        usersCandidateVelocity(1));
        
        % update average velocity 
        usersCandidateAvgVelocity = .9*usersCandidateAvgVelocity + ...
            0.1*norm(usersCandidateVelocity);
                                    
        % update distance to cells
        usersCandidateDistanceToCells = sqrt(...
            (usersCandidatePosition(:,1)-cellsPosition(:,1)).^2 + ...
            (usersCandidatePosition(:,2)-cellsPosition(:,2)).^2 );
        
        % update rx power
        usersCandidatePowFromCells = zeros(nCells,1);
        for c = 1 : nCells
            d = max(rd(cellsId(c)),usersCandidateDistanceToCells(c));
            usersCandidateDistanceToCells(c) = d;
            usersCandidatePowFromCells(c) = P(cellsId(c))-p(cellsId(c))-...
                10*e(cellsId(c))*log10(d/rd(cellsId(c)));
        end
        
        % update angle to cells     
        usersCandidateAngleFromCells = ...
            atan2((usersCandidatePosition(:,2)-cellsPosition(:,2)) , ...
                  (usersCandidatePosition(:,1)-cellsPosition(:,1)) );
        
        % handover algorithm         
        DataP = estDataP(...
                      usersCandidatePowFromCells(1:nCells),...
                      0*usersCandidatePowFromCells(1),...
                      usersCandidateDirection*180/pi,...
                      usersCandidateAngleFromCells(1:nCells)*180/pi,...
                      cellsLoad(1:nCells),nCells,ksi,sigma_v,...
                      usersCandidateAvgVelocity,...
                      usersCandidateDistanceToCells(1:nCells));
        
        Ex = usersCandidatePowFromCells(1:nCells) - ...
              0*usersCandidatePowFromCells(1);  % only power
       
        [~,cellSelectedProposed] = max(DataP);
        [~,cellSelectedEx]      = max(Ex);        
        
        % calculate SIR for the candidate user in the selected cell
        % it is assumed that all macro cells are working with the same frequency
        % it is assumed that all base stations are working with the same frequency
        % it is assumed that all femto cells are working with the same frequency
        cellIdProposed = cellsId(cellSelectedProposed);
        cellIdEx      = cellsId(cellSelectedEx);
        % find power of the interferers power
        IntProposed = usersCandidatePowFromCells(cellsId == cellIdProposed & (1:nCells).' ~= cellSelectedProposed);
        IntEx      = usersCandidatePowFromCells(cellsId == cellIdEx      & (1:nCells).' ~= cellSelectedEx);
        
        SIR_Proposed = 10.^(0.1*usersCandidatePowFromCells(cellSelectedProposed)) ./ sum(...
            10.^(0.1*IntProposed));
        SIR_Ex      = 10.^(0.1*usersCandidatePowFromCells(cellSelectedEx)) ./ sum(...
            10.^(0.1*IntEx));
        
        % store information 
        simT(simCnt,:)           = simTime;
        simPos(simCnt,:)         = usersCandidatePosition;
        simAvgV(simCnt,:)        = usersCandidateAvgVelocity;
        simDir(simCnt,:)         = usersCandidateDirection;
        simPow(simCnt,:)         = usersCandidatePowFromCells;
        simDist(simCnt,:)        = usersCandidateDistanceToCells;
        simAngle(simCnt,:)       = usersCandidateAngleFromCells;
        simEstP(simCnt,1:nCells) = DataP;
        simEx(simCnt,1:nCells)  = Ex;
        simThrPro(simCnt,1)      = log2(1+SIR_Proposed);
        simThrEx(simCnt,1)      = log2(1+SIR_Ex);
        simLoadPro(simCnt,1)     = (cellsHist(cellSelectedProposed)+1)/cellsMaxUser(cellIdProposed);
        simLoadEx(simCnt,1)     = (cellsHist(cellSelectedEx)+1)/cellsMaxUser(cellIdEx);
        % update time
        simTime = simTime + simStepTime;
        
        % update network figure
        if Ctrl.ifShow_Figure1_Update  == true
            hCandidate.XData = usersCandidatePosition(1)/1e3;
            hCandidate.YData = usersCandidatePosition(2)/1e3;
            drawnow
        end
    end
    if Ctrl.ifShow_Figure1_NetSetup
        figure(1)
            hold on
            plot(simPos(:,1)/1e3,simPos(:,2)/1e3,'k--')
            f1.CloseRequestFcn = 'closereq';
    end
%--------------------------------------------------------------------------
%% plot results -----------------------------------------------------------
    color = 'yrgkmcb';
    line1 = '-------';
    line2 = ' -. -. ';
    figure(3)        
    clf
        subplot(2,1,1)
        [~,cellSelected] = max(simEstP.');
        stairs(simT,cellSelected,'linewidth',1.5)
        xlabel('simulation time [s]')
        ylabel('selected cell')
        title('proposed paper')
        grid on
        
        %ylim([0 L+1])
        
        subplot(2,1,2)
        [~,cellSelected] = max(simEx.');
        stairs(simT,cellSelected,'linewidth',1.5)
        xlabel('simulation time [s]')
        ylabel('selected cell')
        title('Existing')
        grid on
        %ylim([0 L+1])
  
    figure(4)
    clf

        plot(simT,simLoadPro,'b')
        hold on 
        plot(simT,simLoadEx,'r--')
        grid on 
        grid minor
        xlabel('simulation time [s]')
        ylabel('cell load')
        legend('proposed','Existing')
        ylim([0 1])
%-------------------------------------------------------------------------- 
   