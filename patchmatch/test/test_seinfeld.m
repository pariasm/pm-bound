% images
u1 = double(imread('seinfeld-129.png')) ;
u2 = double(imread('seinfeld-130.png'));

% u1 = u1(10:90,30:170,:);
% u2 = u2(10:90,30:170,:);

% sizes
sz1 = [size(u1,1),size(u1,2)];
sz2 = [size(u2,1),size(u2,2)];


% patch size <------------------------------
psz = [7,7];

hpsz = floor(psz(1)/2);

% masks
msk1 = ones(sz1); 
msk1(1:hpsz,:) = 0; msk1(end-hpsz+1:end,:) = 0;
msk1(:,1:hpsz) = 0; msk1(:,end-hpsz+1:end) = 0;
hl = sum(msk1(:));

msk2 = ones(sz2);
msk2(1:hpsz,:) = 0; msk2(end-hpsz+1:end,:) = 0;
msk2(:,1:hpsz) = 0; msk2(:,end-hpsz+1:end) = 0;


% call nnf
prms.mxit = 50;    % nnf iterations
prms.list = 10;    % number of candidates
prms.psz = psz;
window = gauss_kernel(prms.psz,30);
window = ones(prms.psz);

addpath ../bin
tic
[nnfy0,nnfx0] = nnfield(u1,msk1,u2,msk2,window,prms);
toc
%[nnfy0,nnfx0] = nnfield_extpatches(u1,msk1,u2,msk2,window,prms);


% display results
figure, imagesc([nnfx0(:,:,1) nnfy0(:,:,1)]), axis equal
title('1st nearest neighbor field (x and y components)')

% figure, imagesc([nnfx0(:,:,2) nnfy0(:,:,2)]), axis equal
% title('2nd nearest neighbor field (x and y components)')
% 
% figure, imagesc([nnfx0(:,:,3) nnfy0(:,:,3)]), axis equal
% title('3rd nearest neighbor field (x and y components)')



