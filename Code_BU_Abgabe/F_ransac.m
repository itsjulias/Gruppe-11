function [Korrespondenzen_robust] = F_ransac(Korrespondenzen, varargin)
    % Diese Funktion implementiert den RANSAC-Algorithmus zur Bestimmung von
    % robusten Korrespondenzpunktpaaren
%     %% Input parser
%     p = inputParser;
%     valid_epsilon = @(x) isnumeric(x) & (x<1) & (x>0);
%     valid_p = @(x) isnumeric(x) & (x>0) & (x<1);
%     valid_tolerance = @(x) isnumeric(x);
%     addOptional(p,'epsilon',0.50,valid_epsilon);
%     addOptional(p,'p',0.5,valid_p);
%     addOptional(p,'tolerance',0.01,valid_tolerance);
%     parse(p, varargin{:});
%     
%     x1_pixel = [Korrespondenzen(1:2,:); ones(1,length(Korrespondenzen(1,:)))];
%     x2_pixel = [Korrespondenzen(3:4,:); ones(1,length(Korrespondenzen(1,:)))];
%     epsilon = p.Results.epsilon;
%     tolerance = p.Results.tolerance;
%     p = p.Results.p;
%     
%     
%     %% RANSAC Algorithmus Vorbereitung
%     k = 8;
%     s = log(1-p)/log(1-(1-epsilon)^k);
%     largest_set_size = 0;
%     largest_set_dist = inf;
%     largest_set_F = zeros(3,3);
%     
%    for i=1:floor(s)
%         idx = ceil(rand(1,k)*length(Korrespondenzen(1,:)));
%         F = achtpunktalgorithmus(Korrespondenzen(:,idx));
%         sd = sampson_dist(F, x1_pixel, x2_pixel);
%         
%         consensus_set = Korrespondenzen(:,sd<tolerance);
%         num_consensus_set = length(consensus_set(1,:));
%         dist_consensus_set = sum(sd(sd<tolerance));
%         if (num_consensus_set == largest_set_size)
%             if(dist_consensus_set < largest_set_dist)
%                 largest_set_size = num_consensus_set;
%                 largest_set_dist = dist_consensus_set;
%                 Korrespondenzen_robust = consensus_set;
%                 largest_set_F = F;
%             end
%         elseif(num_consensus_set > largest_set_size)
%             largest_set_size = num_consensus_set;
%             largest_set_dist = dist_consensus_set;
%             Korrespondenzen_robust = consensus_set;
%             largest_set_F = F;
%         end
% 
%     end

    % Input parser
    par = inputParser;
    % Standardwerte definieren
    defaultepsilon = 0.5;
    defaultp = 0.5;
    defaulttolerance = 0.01;
    
    % Bedingungen für gültige Inputs definieren:
    validepsilon = @(x) isnumeric(x) && x<1 && x>0;
    validp = @(x) isnumeric(x) && x<1 && x>0;
    validtolerance = @(x) isnumeric(x);
   
    % Optionale Parameterübergabe definieren:
    addOptional(par,'epsilon',defaultepsilon, validepsilon);
    addOptional(par,'p',defaultp, validp);
    addOptional(par,'tolerance',defaulttolerance, validtolerance);
    
    % Parsen durchführen.
    parse(par,varargin{:});
    x1_pixel = [Korrespondenzen(1:2,:); ones(1,length(Korrespondenzen))];
    x2_pixel = [Korrespondenzen(3:4,:); ones(1,length(Korrespondenzen))];
    epsilon = par.Results.epsilon;p = par.Results.p;
    tolerance = par.Results.tolerance;


%% RANSAC Algorithmus Vorbereitung
    k = 8;
    s = log(1-p)/log(1-(1-epsilon)^k);
    largest_set_size = 0;
    largest_set_dist = inf;
    largest_set_F = zeros(3);

    
        %% RANSAC Algorithmus
    % Indizes des besten Consesus-Sets. (Spaltenindizes von Korrespondenzen)
    best_ind_cons_set = 1;
    i = 1;
    while i <= s
        % 8 zufällige Korrespondenzpunktpaare auswählen.
        kor_rand = Korrespondenzen(:,randi(size(Korrespondenzen,2),[8,1]));
        % Fundamentale Matrix schätzen.
        F = achtpunktalgorithmus(kor_rand);
        % Sampson-Distanz für alle Korrespondenzpunktpaare berechnen.
        sd = sampson_dist(F, x1_pixel, x2_pixel);
        % Indizes des aktuellen Consensus-Sets bestimmen.
        ind_cons_set = find(sd < tolerance);
        %  Aktuelle Anzahl der im Consensus-Set enthaltenen Paare und akutelle Set-Distanz bestimmen.
        anz_cons = length(ind_cons_set);
        set_dist = sum(sd(ind_cons_set));
        % Vergleichen Sie diese mit der bisher groessten Set-Groesse gespeichert in largest_set_size.
        % Ist das aktuelle Set groesser wird dieses als des neue groesste uebernommen.
        % Sind beide gleich gross vergleichen Sie die Set-Distanz mit der des bisher groessten Sets gespeichter in largest_set_size.
        % Das Set mit der kleineren absolut Sampson-Distanz soll das neue groesste Set werden.
        if anz_cons > largest_set_size || (anz_cons == largest_set_size && set_dist < largest_set_dist)
            largest_set_size = anz_cons;
            largest_set_dist = set_dist;
            % Indizes des besten Consesus-Sets abspeichern. (Spaltenindizes von Korrespondenzen)
            best_ind_cons_set = ind_cons_set;
            % Speichern Sie auch die Fundamentalmatrix F mit der das neue groesste Set generiert wurde in largest_set_F.
            largest_set_F = F;
        end
        i = i+1;
    end
    Korrespondenzen_robust = Korrespondenzen(:,best_ind_cons_set);
    
function sd = sampson_dist(F, x1_pixel, x2_pixel)
    % Diese Funktion berechnet die Sampson Distanz basierend auf der
    % Fundamentalmatrix F
    % Vorher programmierten Dach-Operator verwenden.
    e3_dach = dach([0 0 1]);
    f1 = e3_dach*F*x1_pixel;
    f2 = (x2_pixel'*F*e3_dach)';
    sd = (sum( (x2_pixel'*F)'.*x1_pixel ).^2) ./ (sum(f1.^2)+sum(f2.^2));
end

function W = dach(w)
    % Diese Funktion implementiert den ^-Operator.
    % Sie wandelt einen 3-Komponenten Vektor in eine
    % schiefsymetrische Matrix um.
    W = [0 -w(3) w(2);
        w(3) 0 -w(1);
        -w(2) w(1) 0];
end




end
