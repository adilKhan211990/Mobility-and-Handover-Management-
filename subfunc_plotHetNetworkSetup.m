function [h,hcu] = subfunc_plotHetNetworkSetup(L,rc,cellsPosition,usersPositions,usersCandidate)
figure(1)
        clf
    % set a color
        Colors = [1 1 1];
        hCnt   = 1;
    % plot cells' central positions
        h(hCnt) = plot(cellsPosition(:,1)/1e3*NaN,cellsPosition(:,2)/1e3*NaN,'xr','MarkerSize',4);
        hCnt = hCnt + 1;
        hold on
    for cnt = 1 : L
        % plot cell boundaries
        boundaries = rc(cnt)*[sin([0,.01:.01:2*pi,2*pi]).',...
                              cos([0,.01:.01:2*pi,2*pi]).'];
        % Note that 1e3 is used to covert from m to Km
        h(hCnt) =  fill(cellsPosition(cnt,1)/1e3 + boundaries(:,1)/1e3,...
                        cellsPosition(cnt,2)/1e3 + boundaries(:,2)/1e3,...
                        Colors,'FaceAlpha',0.8);
        hCnt = hCnt + 1;
    end                     
    for cnt = 1 : L
        boundaries = rc(cnt)*[sin([0,.01:.01:2*pi,2*pi]).',...
                              cos([0,.01:.01:2*pi,2*pi]).'];
        % Note that 1e3 is used to covert from m to Km
        h(hCnt) = plot(cellsPosition(cnt,1)/1e3 + boundaries(:,1)/1e3,...
                       cellsPosition(cnt,2)/1e3 + boundaries(:,2)/1e3,'k:');
        hCnt = hCnt + 1;           
    end
    for cnt = 1 : L 
        text(cellsPosition(cnt,1)/1e3-rc(cnt)*sqrt(3)/20/1e3,...
             cellsPosition(cnt,2)/1e3-0*rc(cnt)*sqrt(3)/20/1e3,...
             num2str(cnt),'FontSize',8,'FontWeight','Bold','Color','k')            
    end

    % Plot users positions
    h(hCnt) = plot(usersPositions(usersCandidate,1)/1e3,...
                   usersPositions(usersCandidate,2)/1e3,'ro','MarkerFaceColor','r');
    hcu  = hCnt;
    hCnt = hCnt + 1;           
    % Plot users positions
    h(hCnt) = plot(usersPositions(:,1)/1e3,usersPositions(:,2)/1e3,'r.');
    axis equal
    xlabel('x-axis [Km]')
    ylabel('y-axis [Km]')
    drawnow 
end