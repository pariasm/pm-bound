% load a pair of images --------------------------------------------------------
u1 = mean(double(imread('001.png')),3);
u2 = mean(double(imread('002.png')),3);

u1 = u1(1:100,1:100);
u2 = u2(1:100,1:100);

% u12 = randn(20,20);
% u1 = 128 + 20*[randn(10, 50); randn(20,10) u12 randn(20,20); randn(20, 50)];
% u2 = 128 + 20*[randn(50,100); randn(20,50) u12 randn(20,30); randn(30,100)];
% %u12 = randn(5,5);
% %u1 = 128 + 20*[randn(10, 50); randn( 5,10) u12 randn( 5,35); randn(35, 50)];
% %u2 = 128 + 20*[randn(50,100); randn( 5,50) u12 randn( 5,45); randn(45,100)];
% %u1 = 128 + 20*randn(50,50);
% %u2 = 128 + 20*randn(100,100);
% 
% %u1 = imfilter(u1, fspecial('gaussian',5,2));
% %u2 = imfilter(u2, fspecial('gaussian',5,2));
% 
% %crop = [61 200 21 140];
% crop = [1 50 1 50];
% u1 = u1(crop(1):crop(2),crop(3):crop(4));
% 
% sz1 = [size(u1,1),size(u1,2)];
% sz2 = [size(u2,1),size(u2,2)];

% set parameters ---------------------------------------------------------------
prms.psz   = [5,5];  % patch size
prms.list  = 1;      % number of nearest neighbors
prms.iters = 250;     % patchmatch iterations
prms.rsz = 20;
prms.transition_kernel = 'acceptance';
%prms.transition_kernel = 'regular';
prms.nradii = 1; % 1 means uniform search ~ 0 means patchmatch search

% select a set of points in image 1, and some energy levels --------------------
%xrange = [5,12,14,13]; % x are rows
%yrange = [5,12,14,13]; % y are cols
xrange = [10,50]; % x are rows
yrange = [10,50]; % y are cols
[xx,yy] = meshgrid(yrange,xrange);
points = [yy(:) xx(:)]; clear xx yy
values = [1, 2, 5, 10, 15];
npoints = length(points);
nvalues = length(values);

% compute theoretical convergence bound at selected points ---------------------
Pi = zeros(npoints, nvalues);
for i = 1:npoints,
	Pi(i,:) = initial_probas(u1, u2, points(i,:), values, prms);
end

if 1,
	C0 = zeros(npoints, nvalues);
	C1 = zeros(npoints, nvalues);
    C1_2 = zeros(npoints, nvalues);
    C2 = zeros(npoints, nvalues);
    C3 = zeros(npoints, nvalues);
	Pi = zeros(npoints, nvalues);
	for i = 1:npoints,
		prms.rsz = 10;
		[C0x, ll] = bound(u1, u2, points(i,:), values, prms);
		C0(i,:) = prod(C0x);
        C2x = expected_bound(u1, u2, points(i,:), values, prms);
		C2(i,:) = prod(C2x);
		[C1x, C1_2x, Pix, leps] = tighter_bound(u1, u2, points(i,:), values, prms);
		C1(i,:) = prod(prod(C1x));
        C1_2(i,:) = prod(prod(C1_2x));
        %C3x = tight_bound(u1, u2, points(i,:), values, prms);
		%C3(i,:) = prod(prod(C3x));
		Pi(i,:) = Pix';
	end
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
figure

% generate theoretical energy decay
Pt0 = zeros(size(P));
Pt1 = zeros(size(P));
Pt1_2 = zeros(size(P));
Pt2 = zeros(size(P));
Pt3 = zeros(size(P));
Pt0(:,1,:) = Pi;
Pt1(:,1,:) = Pi;
Pt1_2(:,1,:) = Pi;
Pt2(:,1,:) = Pi;
Pt3(:,1,:) = Pi;
for i = 2:prms.iters+1,
	Pt0(:,i,:) = C0.*squeeze(Pt0(:,i-1,:));
	Pt1(:,i,:) = C1.*squeeze(Pt1(:,i-1,:));
    Pt1_2(:,i,:) = C1_2.*squeeze(Pt1_2(:,i-1,:));
    Pt2(:,i,:) = C2.*squeeze(Pt2(:,i-1,:));
    Pt3(:,i,:) = C3.*squeeze(Pt3(:,i-1,:));
end
hold on;
plot([0:prms.iters],squeeze(Pt0(end,:,:,1)),'.-')
plot([0:prms.iters],squeeze(Pt1(end,:,:,1)),':')
plot([0:prms.iters],squeeze(Pt2(end,:,:,1)),'.')
plot([0:prms.iters],squeeze(Pt3(end,:,:,1)),'.')
% 
% figure
% hold on
% plot([0:prms.iters],squeeze(P(end,:,1:3,1)))
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot([0:prms.iters],squeeze(Pt1_2(end,:,1:3,1)))%,'MarkerIndices',1:10:length(Pt1(end,:,1)))
% ax.ColorOrderIndex = 1;
% plot([0:5:prms.iters],squeeze(Pt1_2(end,1:5:end,1:3,1)),'*')%,'MarkerIndices',1:10:length(Pt1(end,:,1)))
% ax.ColorOrderIndex = 1;
% plot([0:prms.iters],squeeze(Pt0(end,:,1:3,1)),'')%,'MarkerIndices',1:10:length(Pt0(end,:,1)))
% ax.ColorOrderIndex = 1;
% plot([0:5:prms.iters],squeeze(Pt0(end,1:5:end,1:3,1)),'o')%,'MarkerIndices',1:10:length(Pt1(end,:,1)))

blue=[0    0.4470    0.7410];
red=[0.8500    0.3250    0.0980];
yellow=[0.9290    0.6940    0.1250];

%% Figure with true bound
figure
hold on
plot([0:prms.iters],squeeze(P(end,:,1:3,1)))
line_fewer_markers([0:prms.iters],squeeze(Pt1(end,:,1,1)), 20, '*-', 'spacing', 'curve', 'color', blue)
line_fewer_markers([0:prms.iters],squeeze(Pt1(end,:,2,1)), 20, '*-', 'spacing', 'curve', 'color', red)
line_fewer_markers([0:prms.iters],squeeze(Pt1(end,:,3,1)), 20, '*-', 'spacing', 'curve', 'color', yellow)


line_fewer_markers([0:prms.iters],squeeze(Pt0(end,:,1,1)), 20, 'o-', 'spacing', 'curve', 'color', blue)
line_fewer_markers([0:prms.iters],squeeze(Pt0(end,:,2,1)), 20, 'o-', 'spacing', 'curve', 'color', red)
line_fewer_markers([0:prms.iters],squeeze(Pt0(end,:,3,1)), 20, 'o-', 'spacing', 'curve', 'colo', yellow)


%% Figure with "expected" bound
figure
hold on
plot([0:prms.iters],squeeze(P(end,:,1:3,1)))
line_fewer_markers([0:prms.iters],squeeze(Pt1_2(end,:,1,1)), 20, '*-', 'spacing', 'curve', 'color', blue)
line_fewer_markers([0:prms.iters],squeeze(Pt1_2(end,:,2,1)), 20, '*-', 'spacing', 'curve', 'color', red)
line_fewer_markers([0:prms.iters],squeeze(Pt1_2(end,:,3,1)), 20, '*-', 'spacing', 'curve', 'color', yellow)


line_fewer_markers([0:prms.iters],squeeze(Pt2(end,:,1,1)), 20, 'o-', 'spacing', 'curve', 'color', blue)
line_fewer_markers([0:prms.iters],squeeze(Pt2(end,:,2,1)), 20, 'o-', 'spacing', 'curve', 'color', red)
line_fewer_markers([0:prms.iters],squeeze(Pt2(end,:,3,1)), 20, 'o-', 'spacing', 'curve', 'colo', yellow)