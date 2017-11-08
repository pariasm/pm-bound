
im1 = double(imread('~/Work/denoising/data/derf/mobile_mono/001.png'));
im2 = double(imread('~/Work/denoising/data/derf/mobile_mono/031.png'));

% compute multiple energies
params.iters = 5;
params.psz = 4;
params.k = 3;
params.verbose = 0;
params.reconstruct = 0;

for i = 1:20,
	disp(i)
	[ns,ds(:,:,:,i)] = patchmatch(im1(51:150, 51:170, :), im2, params);
end

