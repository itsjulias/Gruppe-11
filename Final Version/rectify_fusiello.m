function [TL,TR,R_rect] = rectify_fusiello(PoL,PoR,K,R_rect,dL,dR)
%RECTIFY_FUSIELLO Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    dL = [0;0];
    dR = [0;0];
    
    c1 = -inv(PoL(:,1:3))*PoL(:,4); % is equal to (0,0,0)
    % Epipol e1
    c2 = -inv(PoR(:,1:3))*PoR(:,4); % is equal to -R'*T

    v1 = (c2-c1); % new x-axis along translation vector
    v2 = cross([0 0 1]',v1); % new y-axis is orthogonal to x- and z-axis
    v3 = cross(v1,v2); % z-axis is crossproduct of x-/

    % new extrinsic (translation unchanged)
    R_rect = [v1'/norm(v1)
              v2'/norm(v2)
              v3'/norm(v3)];
end


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
% Linkes Bild: x' = KL*R_rect*inv(K)*x
% Rechtes Bild: x' = KR*R_rect*R'*inv(K)*x
% Daraus ergeben sich die folgenden Homographien:
% TL = KL*R_rect*inv(K)
% TR = KR*R_rect*inv(K*R) = KR*R_rect*inv(R)*inv(K)
TL = KL * R_rect * inv(PoL(1:3,1:3));
TR = KR * R_rect * inv(PoR(1:3,1:3));


end

