% load a pair of images --------------------------------------------------------
u1 = mean(double(imread('~/Work/denoising/data/derf/bus/001.png')),3);
u2 = mean(double(imread('~/Work/denoising/data/derf/bus/010.png')),3);

u12 = randn(24,24);
u1 = 128 + 20*u12;%[randn(10, 50); randn(20,10) u12 randn(20,20); randn(20, 50)];
u2 = 128 + 20*[randn(48,100); randn(24,48) u12 randn(24,28); randn(28,100)];
%u12 = randn(5,5);
%u1 = 128 + 20*[randn(10, 50); randn( 5,10) u12 randn( 5,35); randn(35, 50)];
%u2 = 128 + 20*[randn(50,100); randn( 5,50) u12 randn( 5,45); randn(45,100)];
%u1 = 128 + 20*randn(50,50);
%u2 = 128 + 20*randn(100,100);

%u1 = imfilter(u1, fspecial('gaussian',5,2));
%u2 = imfilter(u2, fspecial('gaussian',5,2));

%crop = [61 200 21 140];
crop = [1 size(u1,1) 1 size(u1,2)];
u1 = u1(crop(1):crop(2),crop(3):crop(4));

sz1 = [size(u1,1),size(u1,2)];
sz2 = [size(u2,1),size(u2,2)];

% set parameters ---------------------------------------------------------------
prms.psz   = [5,5];  % patch size
prms.list  = 1;      % number of nearest neighbors
prms.iters = 250;     % patchmatch iterations
prms.rsz = 20;

% select a set of points in image 1, and some energy levels --------------------
%xrange = [5,12,14,13]; % x are rows
%yrange = [5,12,14,13]; % y are cols
%xrange = [3,22]; % x are rows
%yrange = [3,22]; % y are cols
%[xx,yy] = meshgrid(yrange,xrange);
%points = [yy(:) xx(:)]; clear xx yy
points = [3,3; 22, 22];
npoints = length(points);
%values = [1, 2, 5, 10, 15];
values = [.5];
nvalues = length(values);
		prms.rsz = 20;

% compute theoretical convergence bound at selected points ---------------------
Pi = zeros(npoints, nvalues);
for i = 1:npoints,
	Pi(i,:) = initial_probas(u1, u2, points(i,:), values, prms);
end

if 0,
	C0f = zeros(npoints, nvalues);
	C1f = zeros(npoints, nvalues);
	for i = 1:npoints,
		C0x = bound(u1, u2, points(i,:), values, prms);
		C0f(i,:) = prod(C0x);
		C1x = tighter_bound(u1, u2, points(i,:), values, prms);
		C1f(i,:) = prod(prod(C1x));
	end

	u1 = fliplr(flipud(u1));
	u2 = fliplr(flipud(u2));
	points = repmat(sz1+1,[npoints,1]) - points;

	C0b = zeros(npoints, nvalues);
	C1b = zeros(npoints, nvalues);
	for i = 1:npoints,
		C0x = bound(u1, u2, points(i,:), values, prms);
		C0b(i,:) = prod(C0x);
		C1x = tighter_bound(u1, u2, points(i,:), values, prms);
		C1b(i,:) = prod(prod(C1x));
	end

	u1 = fliplr(flipud(u1));
	u2 = fliplr(flipud(u2));
	points = repmat(sz1+1,[npoints,1]) - points;
end


% compute empirical convergence bound at selected points -----------------------
if 0,
	trials = 500;
	P = [];
	varP = [];
	%exdata = load('../test/bus_nnf_exhaustive');
	%sc0 = exdata.scex(crop(1):crop(2),crop(3):crop(4))/prod(prms.psz); clear exdata;
	sc0 = [];
	[P, varP, sc0] = empirical_bound(u1, u2, prms, points, values, trials, ...
	                                 sc0, P, varP, trials);
end

% visualize results comparing bounds -------------------------------------------

% generate theoretical energy decay
Pt0 = zeros(npoints, prms.iters + 1, nvalues, prms.list);
Pt1 = zeros(npoints, prms.iters + 1, nvalues, prms.list);
Pt0(:,1,:) = Pi;
Pt1(:,1,:) = Pi;
for i = 2:prms.iters+1,
	if mod(i,2)
		Pt0(:,i,:) = C0b.*squeeze(Pt0(:,i-1,:));
		Pt1(:,i,:) = C1b.*squeeze(Pt1(:,i-1,:));
	else
		Pt0(:,i,:) = C0f.*squeeze(Pt0(:,i-1,:));
		Pt1(:,i,:) = C1f.*squeeze(Pt1(:,i-1,:));
	end
end
