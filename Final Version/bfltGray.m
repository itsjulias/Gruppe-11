% Implements bilateral filtering for grayscale images.
function B = bfltGray(A,w,sigma_d,sigma_r)
% Quelle: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/12191/versions/1/previews/Bilateral%20Filtering/bfilter2.m/index.html
% A: grayscale image, double precision matrix of size NxMx1 with normalized values in the closed interval [0,1];
% w: half-size of the Gaussian bilateral filter window,
% sigma: standard diviation of the bilateral filter; 
%   sigma_d: spatial-domain standard deviation,
%   sigma_r: intensity-domain standard deviation

% Pre-compute Gaussian distance weights.
[X,Y] = meshgrid(-w:w,-w:w);
G = exp(-(X.^2+Y.^2)/(2*sigma_d^2));

% Create waitbar.
h = waitbar(0,'Applying bilateral filter...');
set(h,'Name','Bilateral Filter Progress');

% Apply bilateral filter.
dim = size(A);
B = zeros(dim);
for i = 1:dim(1)
    for j = 1:dim(2)
        
        % Extract local region.
        iMin = max(i-w,1);
        iMax = min(i+w,dim(1));
        jMin = max(j-w,1);
        jMax = min(j+w,dim(2));
        I = A(iMin:iMax,jMin:jMax);
        
        % Compute Gaussian intensity weights.
        H = exp(-(I-A(i,j)).^2/(2*sigma_r^2));
        
        % Calculate bilateral filter response.
        F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1);
        sum_F(i,j) = 1/sum(F(:));
        sum_F_I(i,j) = sum(F(:).*I(:));
%         B(i,j) = sum(F(:).*I(:))/sum(F(:));
        
    end
    waitbar(i/dim(1));
end
B = sum_F_I.*sum_F;
close(h)
end
