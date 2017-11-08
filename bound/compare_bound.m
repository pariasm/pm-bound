% load a pair of images --------------------------------------------------------
u1 = mean(double(imread('~/Work/denoising/data/derf/bus/001.png')),3);
u2 = mean(double(imread('~/Work/denoising/data/derf/bus/010.png')),3);

u12 = randn(20,20);
u1 = 128 + 20*[randn(10, 50); randn(20,10) u12 randn(20,20); randn(20, 50)];
u2 = 128 + 20*[randn(50,100); randn(20,50) u12 randn(20,30); randn(30,100)];
%u1 = 128 + 20*randn(50,50);
%u2 = 128 + 20*randn(100,100);

%u1 = imfilter(u1, fspecial('gaussian',5,2));
%u2 = imfilter(u2, fspecial('gaussian',5,2));

%crop = [61 200 21 140];
crop = [1 50 1 50];
u1 = u1(crop(1):crop(2),crop(3):crop(4));

sz1 = [size(u1,1),size(u1,2)];
sz2 = [size(u2,1),size(u2,2)];

% set parameters ---------------------------------------------------------------
prms.psz   = [5,5];  % patch size
prms.list  = 1;      % number of nearest neighbors
prms.iters = 200;     % patchmatch iterations

% select a set of points in image 1, and some energy levels --------------------
xrange = [5:5:25]; % x are rows
yrange = [5:5:25]; % y are cols
[xx,yy] = meshgrid(yrange,xrange);
points = [yy(:) xx(:)]; clear xx yy
values = [1, 2, 5, 10, 15];
npoints = length(points);
nvalues = length(values);

% compute theoretical convergence bound at selected points ---------------------
C0 = zeros(npoints, nvalues);
for i = 1:length(points),
	prms.rsz = 10;
	tmp = bound(u1, u2, points(i,:), values, prms);
	C0(i,:) = prod(tmp);
end

% compute empirical convergence bound at selected points -----------------------
trials = 4000;
P = [];
%exdata = load('../test/bus_nnf_exhaustive');
%sc0 = exdata.scex(crop(1):crop(2),crop(3):crop(4))/prod(prms.psz); clear exdata;
sc0 = [];
[P,sc0] = empirical_bound(u1, u2, prms, points, values, trials, sc0, P, trials);

% visualize results comparing bounds -------------------------------------------

% generate theoretical energy decay
Pt = zeros(size(P));
Pt(:,1,:) = P(:,1,:);
for i = 2:prms.iters,
	Pt(:,i,:) = C0.*squeeze(Pt(:,i-1,:));
end
