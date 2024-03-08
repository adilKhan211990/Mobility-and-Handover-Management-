function P = estDataP(RSRP,RSRP_Src,Theta_UE,Theta_UE_C,CNr,nCells,ksi,sigma_v,V_avg,Dc)
    Ctrl.ifRandomJc = false;
    % find RSRP_g and dc
    RSRP_g      = RSRP(:) - RSRP_Src;
    dc          = cos((Theta_UE(:)-Theta_UE_C(:))*pi/180);    % value of UE's direction
    [~,Ind]     = max(RSRP_g);
    dc(Ind)     = -1/2;    % the cell with maximum power in the cell that user is located in
    
    % find random sojurn time in cell
    Jc    = zeros(nCells-0,1);
    dt    = 1e-3;  % time resolution
    tEnd  = 200;   % final time
    sigma = 1/2;   % normlization factor in (1)
    for n = 1 : nCells-0
        % find distribution explained in 1 for cell n
        f = @(t) ksi*Dc(n)*sigma^-1./(t.^2*sqrt(2*pi)).* ...
            exp(-((Dc(n)./t-V_avg)./(sqrt(2)*sigma_v)).^2);
        
        t = dt:dt:tEnd;
        if Ctrl.ifRandomJc 
            % find the cdf of f
            F = cumsum(f(t).*dt);
            % generate a uniform random number 
            u = rand;
            % find corresponding value from inversfunction 
            % see https://uk.mathworks.com/help/stats/generate-random-numbers-using-the-uniform-distribution-inversion-method.html
            [~,Ind] = min(abs(F-u)); % find the index of the t that minimizes ||F(t) - u||            
        else
            [~,Ind] = max(abs(f(t))); % max of the distribution
        end
        Jc(n)   = 1/1*t(Ind);
        %figure(1)
        %plot(t,F)
        %hold on         
    end
    
    % calculate P
    P = 10.^(0.1*RSRP_g) .* (-dc) .* Jc .* (1-CNr(:));
end