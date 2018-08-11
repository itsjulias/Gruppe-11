function [TL,TR] = rectify_fusiello(PoL,PoR,K,dL,dR)
%RECTIFY_FUSIELLO Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    dL = [0;0];
    dR = [0;0];
end

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
KL(1,3)=KL(1,3)+dL(1);
KL(2,3)=KL(2,3)+dL(2);
KR(1,3)=KR(1,3)+dR(1);
KR(2,3)=KR(2,3)+dR(2);      
      
% new projection matrices
PnL = KL * [R_rect -R_rect*c1 ];
PnR = KR * [R_rect -R_rect*c2 ];

% rectifying image transformation
TL = PnL(1:3,1:3)* inv(PoL(1:3,1:3));
TR = PnR(1:3,1:3)* inv(PoR(1:3,1:3));


end

