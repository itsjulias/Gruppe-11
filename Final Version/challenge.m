%% Computer Vision Challenge
close all;
clear all;

% Groupnumber:
group_number = 11;

% Groupmembers:
members = {'Julia Stroebel','Oliver Doege','Andre Thommessen','Tim Janssen',...
           'Sebastian Huegler'};

% Email-Adress (from Moodle!):
mail = {'sebastian.huegler@tum.de',...
        'ga53tex@mytum.de',...
        'julia.stroebel@tum.de',...
        'oliver.doege@tum.de',...
        'tim.janssen@tum.de'};



%% Load images
img1 = imread([pwd '/img/L1.JPG']);
img2 = imread([pwd '/img/R1.JPG']);

%% Free Viewpoint Rendering
% start execution timer -> tic;
p = 0.2;
tic
img_viewpoint = free_viewpoint(img1,img2,p);

% stop execution timer -> toc;
elapsed_time = toc;

%% Display Output
% Display Virtual View
figure
imshow(img_viewpoint)

