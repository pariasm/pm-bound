% this script test the patch algorithm in an inpainting-like 
% context: a single image (u1 = u2)
%          a region msk1
%          msk2 is the complement of msk1
% 

%images
u1 = double(imread('118.png')) ;
u2 = double(imread('118.png'));

u1 = mean(u1,3);
u2 = mean(u2,3);

randn('state',0);
u1 = u1 + 00*randn(size(u1)); u1(u1 < 0) = 0;
u2 = u2 + 00*randn(size(u2)); u2(u2 < 0) = 0;


% sizes
sz1 = size(u1);
sz2 = size(u2);


% patch size <------------------------------
psz = [5,5];
hpsz = floor(psz(1)/2);

% masks
msk1 = double(imread('118_msk.png'));
hl = sum(msk1(:));

msk2 = 1 - imdilate(msk1,ones(psz));
msk2(1:hpsz,:) = 0; msk2(end-hpsz+1:end,:) = 0;
msk2(:,1:hpsz) = 0; msk2(:,end-hpsz+1:end) = 0;



% call nnf
prms.mxit = 5;    % nnf iterations
prms.list = 3;    % number of candidates
prms.psz = psz;
window = gauss_kernel(prms.psz,3);
window = ones(prms.psz);


addpath ../bin
[nnfy,nnfx] = nnfield(u1,msk1,u2,msk2,window,prms);
[nnfy,nnfx] = nnfield(u1,msk1,u2,msk2,window,prms, nnfy, nnfx);


% display results
figure, imagesc([nnfx(:,:,1) nnfy(:,:,1)]), axis equal
title('1st nearest neighbor field (x and y components)')

figure, imagesc([nnfx(:,:,2) nnfy(:,:,2)]), axis equal
title('2nd nearest neighbor field (x and y components)')

figure, imagesc([nnfx(:,:,3) nnfy(:,:,3)]), axis equal
title('3rd nearest neighbor field (x and y components)')

figure, imagesc([255*msk1 u1 apply_nnf(u1,msk1,nnfx,nnfy,1)]), axis equal, colormap gray, axis off
title('mask, original image and warping done through the nnf')


