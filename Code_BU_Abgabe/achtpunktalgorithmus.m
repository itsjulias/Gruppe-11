function [EF] = achtpunktalgorithmus(Korrespondenzen, K)
    % Diese Funktion berechnet die Essentielle Matrix oder Fundamentalmatrix
    % mittels 8-Punkt-Algorithmus, je nachdem, ob die Kalibrierungsmatrix 'K'
    % vorliegt oder nicht

    pixel_pic1 = [Korrespondenzen(1:2,:); ones(1,length(Korrespondenzen(1,:)))];
    pixel_pic2 = [Korrespondenzen(3:4,:); ones(1,length(Korrespondenzen(1,:)))];
    
    if nargin == 2
        pixel_pic1 = K\pixel_pic1;
        pixel_pic2 = K\pixel_pic2;
    end
    
    A = [];
    
    for i = 1:length(Korrespondenzen(1,:))
       A = [A;kron(pixel_pic1(:,i),pixel_pic2(:,i))']; 
    end
    
    [~,~,V] = svd(A);
    %% Schaetzung der Matrizen
    G = reshape(V(:,end),3,3);
    [Ug,SIGMAg,Vg] =svd(G);
    if(nargin == 2)
        SIGMAg = [1 0 0; 0 1 0; 0 0 0];
    else
        SIGMAg(3,3) = 0;
    end
    EF = Ug*SIGMAg*Vg';
end