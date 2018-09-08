function [repro_error, x2_repro] = rueckprojektion(Korrespondenzen, P1, Image2, T, R, K)
    % Diese Funktion berechnet den mittleren Rueckprojektionsfehler der 
    % Weltkooridnaten P1 aus Bild 1 im Cameraframe 2 und stellt die 
    % korrekten Merkmalskoordinaten sowie die rueckprojezierten grafisch dar.
   
    P2 = K*(R*P1+T);
    x2_repro = (P2)./P2(end,:);
    N =length(Korrespondenzen(1,:));
    N_ = 1:N';
    repro_error = 1/N*...
        sum(sqrt((Korrespondenzen(3,:)-x2_repro(1,:)).^2+(Korrespondenzen(4,:)-x2_repro(2,:)).^2));

    %Plotting
    figure
    imshow(Image2)
    hold all
    plot(x2_repro(1,:),x2_repro(2,:),'*b')
    text(x2_repro(1,:),x2_repro(2,:),num2str(N_),'color','b')
    plot(Korrespondenzen(3,:),Korrespondenzen(4,:),'*r')
    text(Korrespondenzen(3,:),Korrespondenzen(4,:),num2str(N_),'color','r')
end