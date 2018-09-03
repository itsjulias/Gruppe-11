function rotation = euler_rotation(phi_deg, theta_deg, psi_deg)
% http://www.gregslabaugh.net/publications/euler.pdf    
% rads = deg2rad([yaw_deg, pitch_deg, roll_deg]);
    cos_theta = cosd(theta_deg);
    sin_theta = sind(theta_deg);
    cos_psi = cosd(psi_deg);
    sin_psi = sind(psi_deg);
    cos_phi = cosd(phi_deg);
    sin_phi = sind(phi_deg);

    r_x = eye(3);
    r_x(2,2) = cos_psi;
    r_x(3,3) = cos_psi;
    r_x(2,3) = -sin_psi;
    r_x(3,2) = sin_psi;

    r_y = eye(3);
    r_y(1,1) = cos_theta;
    r_y(3,3) = cos_theta;
    r_y(1,3) = sin_theta;
    r_y(3,1) = -sin_theta;

    r_z = eye(3);
    r_z(1,1) = cos_phi;
    r_z(2,2) = cos_phi;
    r_z(1,2) = -sin_phi;
    r_z(2,1) = sin_phi;

    rotation = r_z * r_y * r_x;    
end