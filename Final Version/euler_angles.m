function [theta, psi, phi] = euler_angles(R)
% Quelle: https://www.learnopencv.com/rotation-matrix-to-euler-angles/
sy = sqrt(R(1,1)*R(1,1)+R(2,1)*R(2,1));
if(sy>1e-6)
    psi = atan2d(R(3,2),R(3,3));
    theta = atan2d(-R(3,1),sy);
    phi = atan2d(R(2,1),R(1,1));
else
    psi = atan2d(-R(2,3),R(2,2));
    theta = atan2d(-R(3,1),sy);
    phi = 0;
end
end

% http://www.gregslabaugh.net/publications/euler.pdf
% if(R(3,1) ~= 1 || R(3,1)~=-1)
%     theta1 = -asin(R(3,1));
%     theta2 = pi-theta1;
%
%     psi1 = atan2(R(3,2)/cos(theta1),R(3,3)/theta1);
%     psi2 = atan2(R(3,2)/cos(theta2),R(3,3)/theta2);
%     phi1 = atan2(R(2,1)/cos(theta1),R(1,1)/theta1);
%     phi2 = atan2(R(2,1)/cos(theta2),R(1,1)/theta2);
%
%     theta = theta1;
%     psi = psi1;
%     phi = phi1;
% else
%     phi = 0;
%     if(R(3,1) == -1)
%         theta = pi/2;
%         psi = theta + atan2(R(1,2),R(1,3));
%     else
%         theta = -pi/2;
%         psi = -theta+atan2(-R(1,2),-R(1,3));
%     end
% end

% theta = 180/pi*theta;
% psi = 180/pi*psi;
% phi = 180/pi*phi;