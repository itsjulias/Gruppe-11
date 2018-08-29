function [min_disparity,max_disparity] = get_min_max_disparity(Korrespondenzen,T_rec_1,T_rec_2)
%GET_LINE_SEARCH_INTERVALL Berechnet grob minimale und maximale Disparity
%   um den Liniensuchbereich zurück zu geben.

% Umwandlung in homogene Pixelkoordinaten
anz_kor = size(Korrespondenzen,2);
Kor_L = [Korrespondenzen(1:2,:);ones(1,anz_kor)];
Kor_R = [Korrespondenzen(3:4,:);ones(1,anz_kor)];

% Wende die Rectification Transformation auf die Korrespondenzpunktpaare an
NEW_Kor_L = T_rec_1*Kor_L;
NEW_Kor_R = T_rec_2*Kor_R;

% Bestimme disparities: d = x_rec_L - x_rec_R
% Korrespondenzpunkt des linken Bildes liegt jeweils weiter rechts als der
% des rechten Bildes
disparities = NEW_Kor_L(1,:) - NEW_Kor_R(1,:);
min_disparity = min(disparities);
max_disparity = max(disparities);

end

