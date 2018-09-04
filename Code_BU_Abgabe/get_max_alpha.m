function alpha = get_max_alpha(R)
% GET_ALPHA Berechnet aus der Rotationsmatrix den Winkel der Drehung um die
% y-Achse

% Laut Vorgabe wird nur eine Rotation um die y-Achse durchgeführt, das
% heißt die Rotationsmatrix hat die Form:
% R_v = [cos(alpha),0,sin(alpha);
%       0,1,0;
%       -sin(alpha),0,cos(alpha)];

% Berechne aus jedem Eintrag, der alpha enthält, den Winkel und gebe den
% Mittelwert zurück
alpha = (acosd(R(1,1))+asind(R(1,3))+asind(-R(3,1))+acosd(R(3,3)))/4;

end

