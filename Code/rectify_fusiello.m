function [TL,TR] = rectify_fusiello(PoL,PoR,K,dL,dR)
%RECTIFY_FUSIELLO Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    dL = [0;0];
    dR = [0;0];
end

% Translation ist bereits berechnet worden: Wieso nicht v1 = T bzw.
% im unkalibrierten Fall T durch T' = K*T ersetzen.
% c1 berechnen unnötig, da gleich 0?
c1 = -inv(PoL(:,1:3))*PoL(:,4); % is equal to (0,0,0)
c2 = -inv(PoR(:,1:3))*PoR(:,4); % is equal to -R'*T

v1 = (c2-c1); % new x-axis along translation vector
v2 = cross([0 0 1]',v1); % new y-axis is orthogonal to x- and z-axis
v3 = cross(v1,v2); % z-axis is crossproduct of x-/

% new extrinsic (translation unchanged)
R_rect = [v1'/norm(v1)
          v2'/norm(v2)
          v3'/norm(v3)];

% translate image centers
KL = K;
KR = K;
% Justieren des Ursprungs in x-Richtung
KL(1,3)=KL(1,3)+dL(1);
% Justieren des Ursprungs in y-Richtung
KL(2,3)=KL(2,3)+dL(2);
% Justieren des Ursprungs in x-Richtung
KR(1,3)=KR(1,3)+dR(1);
% Justieren des Ursprungs in y-Richtung
KR(2,3)=KR(2,3)+dR(2);      
      
% new projection matrices
% Überflüssig? (siehe unten)
% Woher kommen diese Gleichungen?
PnL = KL * [R_rect -R_rect*c1 ];
PnR = KR * [R_rect -R_rect*c2 ];

% rectifying image transformation
% Home Position / Umwandlung Weltkoordinaten X in Pixelkoordinaten x: x = K*X
% Rotation by matrix R for left image and by matrix R_rect*R' for right
% image.
% Umwandlung Weltkoordinaten X in rotierte Pixelkoordinaten x' für linkes
% Bild: x' = KL*R_rect*X
% Umwandlung Weltkoordinaten X in rotierte Pixelkoordinaten x' für rechtes
% Bild: x' = KR*R_rect*R'*X
% Einsetzen der Gleichung für die Home Position in die Gleichungen für die
% rotierten Pixelkoordinaten, führt zu:
% Linkes Bild: x' = KL*R_rect*inv(K)*X
% Rechtes Bild: x' = KR*R*_rect*R'*inv(K)*X
% Daraus ergeben sich die folgenden Homographien:
% TL = KL*R_rect*inv(K)
% TR = KR*R_rect*inv(K*R) = KR*R_rect*inv(R)*inv(K) => In rectification.pdf
% steht Formel: mit R_r = R*R_rect ~= R_rect*R' => Fehler in pdf?
TL = PnL(1:3,1:3)* inv(PoL(1:3,1:3));
TR = PnR(1:3,1:3)* inv(PoR(1:3,1:3));


end

