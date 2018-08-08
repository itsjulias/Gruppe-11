%% Computer Vision Challenge
close all;
clear all;
% Groupnumber:
group_number = 11;

% Groupmembers:
% members = {'Max Mustermann', 'Johannes Daten'};
members = {'Julia Stroebel','Oliver Koege','Andre Thommessen','Tim Janßen','Sebastian Huegler'};

% Email-Adress (from Moodle!):
% mail = {'ga99abc@tum.de', 'daten.hannes@tum.de'};
mail = {'sebastian.huegler@tum.de','andrethom@hotmail.de','julia.stroebel@tum.de'};

%% Load images
img1 = imread([pwd '/img/L2.jpg']);
img2 = imread([pwd '/img/R2.jpg']);

%% Free Viewpoint Rendering
% start execution timer -> tic;
p = 0.5;
tic
img_viewpoint = free_viewpoint(img1,img2,p);

% stop execution timer -> toc;
elapsed_time = toc;

%% Display Output
% Display Virtual View
figure
imshow(img_viewpoint)
