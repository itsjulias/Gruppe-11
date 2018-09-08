%% Load images
name_img1 = 'L2';
name_img2 = 'R2';

img1 = imread([pwd '/img/' name_img1 '.JPG']);
img2 = imread([pwd '/img/' name_img2 '.JPG']);

p = [0.2, 0.45, 0.7, 1];
tic
parfor i=1:length(p)
    p(i)
    img_viewpoint = free_viewpoint(img1,img2,p(i));
    figure
    imshow(img_viewpoint);
    imwrite(img_viewpoint,[pwd '/img_created_new/FV_' name_img1 '_' name_img2 '_' num2str(p(i)) '.jpg'])
end
elapsed_time = toc;
