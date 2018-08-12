function depth_map = depth_estimation(img1_rectified,img2_rectified,K,T,varargin)
%DEPTH_MAP Summary of this function goes here
%   Detailed explanation goes here

% % Toleranz Verschiebung der Liniensuche in y-Richtung
% tol = 0;
% Depth Map Matrix: Zeile entspricht x-Koordinate, Spalte entspricht
% y-Koordinate, Wert enspricht der Tiefe.
% Vorbelegung mit Nullen
depth_map = zeros(size(img1_rectified,2),size(img1_rectified,1));
for i = 1:size(img1_rectified,1)
    leng = size(img1_rectified,2);
    zeilenkoord = [1:leng;i*ones(1,leng)];
    kor = punkt_korrespondenzen(img1_rectified,img2_rectified,...
        zeilenkoord,zeilenkoord,...
        'window_length',25,'min_corr',0.5);
    depth = (K*T)/(kor(1:2,:)-kor(3:4,:));
    for j = 1:size(kor,2)
        p = kor(1:2,i);
        depth_map(p(1),p(2)) = depth(j);
    end
end
figure;
[X,Y] = meshgrid(1:size(depth_map,1),1:size(depth_map,2));
surf(X,Y,depth_map)
end

