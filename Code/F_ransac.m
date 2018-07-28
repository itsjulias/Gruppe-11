function [Korrespondenzen_robust] = F_ransac(Korrespondenzen, varargin)
    % Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
    % robusten Korrespondenzpunktpaaren
    %% Input parser
    p = inputParser;
    valid_epsilon = @(x) isnumeric(x) & (x<1) & (x>0);
    valid_p = @(x) isnumeric(x) & (x>0) & (x<1);
    valid_tolerance = @(x) isnumeric(x);
    addOptional(p,'epsilon',0.50,valid_epsilon);
    addOptional(p,'p',0.5,valid_p);
    addOptional(p,'tolerance',0.01,valid_tolerance);
    parse(p, varargin{:});
    
    x1_pixel = [Korrespondenzen(1:2,:); ones(1,length(Korrespondenzen(1,:)))];
    x2_pixel = [Korrespondenzen(3:4,:); ones(1,length(Korrespondenzen(1,:)))];
    epsilon = p.Results.epsilon;
    tolerance = p.Results.tolerance;
    p = p.Results.p;
    
    
    %% RANSAC Algorithmus Vorbereitung
    k = 8;
    s = log(1-p)/log(1-(1-epsilon)^k);
    largest_set_size = 0;
    largest_set_dist = inf;
    largest_set_F = zeros(3,3);
    
   for i=1:floor(s)
        idx = ceil(rand(1,k)*length(Korrespondenzen(1,:)));
        F = achtpunktalgorithmus(Korrespondenzen(:,idx));
        sd = sampson_dist(F, x1_pixel, x2_pixel);
        
        consensus_set = Korrespondenzen(:,sd<tolerance);
        num_consensus_set = length(consensus_set(1,:));
        dist_consensus_set = sum(sd(sd<tolerance));
        if (num_consensus_set == largest_set_size)
            if(dist_consensus_set < largest_set_dist)
                largest_set_size = num_consensus_set;
                largest_set_dist = dist_consensus_set;
                Korrespondenzen_robust = consensus_set;
                largest_set_F = F;
            end
        elseif(num_consensus_set > largest_set_size)
            largest_set_size = num_consensus_set;
            largest_set_dist = dist_consensus_set;
            Korrespondenzen_robust = consensus_set;
            largest_set_F = F;
        end

    end
end