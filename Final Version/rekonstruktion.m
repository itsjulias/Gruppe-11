function [T_cell, R_cell,T,R, d_cell, x1, x2] = rekonstruktion(T1, T2, R1, R2, Korrespondenzen, K, varargin)
%% Input parser
p = inputParser;
valid_do_plot = @(x) islogical(x);
addOptional(p,'do_plot',false,valid_do_plot);

parse(p, varargin{:});
do_plot = p.Results.do_plot;

%% Preparation
T_cell = {T1, T2, T1, T2};
R_cell = {R1, R1, R2, R2};
x1_ = [Korrespondenzen(1:2,:); ones(1,length(Korrespondenzen(1,:)))];
x2_ = [Korrespondenzen(3:4,:); ones(1,length(Korrespondenzen(1,:)))];
x1 = K\x1_;
x2 = K\x2_;
d = zeros(length(Korrespondenzen(1,:)),2);
d_cell = {d,d,d,d};

for i = 1:length(T_cell)
    m1 = mat2cell(cross(x2,R_cell{i}*x1,1),3,ones(1,length(x2)));
    M1 = blkdiag(m1{:});
    m1e = cross(x2,T_cell{i}.*ones(size(x2)));
    M1(:,end+1) = m1e(:);
    
    m2 = mat2cell(cross(x1,R_cell{i}'*x2,1),3,ones(1,length(x1)));
    M2 = blkdiag(m2{:});
    m2e = -cross(x1,R_cell{i}'*T_cell{i}.*ones(size(x1)));
    M2(:,end+1) = m2e(:);
    
    
    [U1,S1,V1] = svd(M1);
    [U2,S2,V2] = svd(M2);
    
    d1 = V1(:,end)/V1(end,end);
    d2 = V2(:,end)/V2(end,end);
    d_cell{i} = [d1(1:end-1),d2(1:end-1)];
end
num_pos_val = sum([d_cell{:}]>0,1);
idx_sum = 1:2:(length(num_pos_val)-1);
[val, idx]=max(num_pos_val(idx_sum)+num_pos_val(idx_sum+1));
T = T_cell{idx};
R = R_cell{idx};
lambda = d_cell{idx};






if(do_plot)
    % Plotten der Raumpunkte
    P1 = lambda(:,1)'.*x1;
    figure
    scatter3(P1(1,:),P1(2,:),P1(3,:),'.');
    hold all;
    text(P1(1,:),P1(2,:),P1(3,:),num2str((1:length(P1(1,:)))'))
    %Berechnung der Kamera-frames/Plotting
    camC1 = [[-0.2;0.2;1],[0.2;0.2;1],[0.2;-0.2;1],[-0.2;-0.2;1]];
    camC2 = R\(camC1-T);
    plot3([camC1(1,:) camC1(1,1)],[camC1(2,:) camC1(2,1)],[camC1(3,:) camC1(3,1)],'b');
    text(camC1(1,end),camC1(2,end),camC1(3,end),'Cam1');
    plot3([camC2(1,:) camC2(1,1)],[camC2(2,:) camC2(2,1)],[camC2(3,:) camC2(3,1)],'r');
    text(camC2(1,end),camC2(2,end),camC2(3,end),'Cam2');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    % xlim([-0.6 1]); zlim([0 5]);
    % Plotten der Kameraframes
    figure
    plot3([camC1(1,:) camC1(1,1)],[camC1(2,:) camC1(2,1)],[camC1(3,:) camC1(3,1)],'b');
    text(camC1(1,end),camC1(2,end),camC1(3,end),'Cam1');
    hold all
    plot3([camC2(1,:) camC2(1,1)],[camC2(2,:) camC2(2,1)],[camC2(3,:) camC2(3,1)],'r');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal
    grid
    
    camup([0 -1 0])
    campos([43 -22 -87])
    grid on
end
end