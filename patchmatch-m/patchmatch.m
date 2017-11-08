% A version PatchMatch algorithm
function [nnf, ds] = patchmatch(im1, im2, params)

% default values
iters = 2; % iterations (each iter comprises forward and backward passes)
psz   = 4; % patch size
k     = 1; % number of nearest neighbors searched for
verbo = 0; % verbose
recon = 0; % if k>0, show im1 reconstructed using the k nearest neighbor

% override defaults with user given values
if isfield(params, 'iters'  ), iters = params.iters;   end
if isfield(params, 'psz'    ), psz   = params.psz;     end
if isfield(params, 'k'      ), k     = params.k;       end
if isfield(params, 'verbose'), verbo = params.verbose; end
if isfield(params, 'recon'  ), recon = params.recon;   end


[h1,w1,c1] = size(im1);
[h2,w2,c2] = size(im2);

h1p = h1 - psz + 1; w1p = w1 - psz + 1; wh1p = h1p * w1p; sz1p = [h1p, w1p, k];
h2p = h2 - psz + 1; w2p = w2 - psz + 1; wh2p = h2p * w2p; sz2p = [h2p, w2p, k];

% patch domain
pdom = [0:psz-1];

% reconstructed im1
if verbo && recon, im12 = zeros(h1, w1, c1); ag12 = zeros(h1, w1,  1); end
	
% random initialization
if verbo, disp('initializing ...'), end
nnf.indx = randi(wh2p, h1p, w1p, k) - 1;
nnf.dist = zeros(sz1p);

% initialize distances
for x1 = 1:w1p,
for y1 = 1:h1p,

	% extract patch at x1 y1
	P = im1(y1+pdom, x1+pdom, :);

	dbest = inf; % for viz
	ibest = inf;
	for i = 1:k,

		% position of candidate
		x2 = floor(nnf.indx(y1, x1, i) / h2p) + 1;
		y2 =   mod(nnf.indx(y1, x1, i) , h2p) + 1;

		% retrieve candidate patch
		Q = im2(y2+pdom, x2+pdom, :);

		% L2 distance
		nnf.dist(y1,x1,i) = mean((P(:) - Q(:)).^2);

		if nnf.dist(y1,x1,i) < dbest, % for viz
			dbest = nnf.dist(y1,x1,i);
			ibest = i;
		end
	end

	if verbo && recon,
		x2 = floor(nnf.indx(y1, x1, recon) / h2p) + 1;
		y2 =   mod(nnf.indx(y1, x1, recon) , h2p) + 1;
		im12(y1+pdom, x1+pdom, :) = im12(y1+pdom, x1+pdom, :) + ...
		                            im2 (y2+pdom, x2+pdom, :);
		ag12(y1+pdom, x1+pdom) = ag12(y1+pdom, x1+pdom) + 1;
	end

end
end

% % sort initial nnf
% [xx1, yy1] = meshgrid(1:w1p, 1:h1p);
% xx1 = repmat(xx1, [1 1 k]);
% yy1 = repmat(yy1, [1 1 k]);
% 
% [nnf.dist, ii] = sort(nnf.dist, 3);
% ii = sub2ind(sz1p, yy1(:), xx1(:), ii(:));
% nnf.indx = reshape(nnf.indx(ii), sz1p); 

ds(:,:,1) = max(nnf.dist,[],3);


if verbo,
	if recon,
		imagesc(uint8(im12./repmat(ag12,[1 1 c1])), [0 255]), colormap gray
	else
		imagesc(sqrt(mean(nnf.dist,3)), [0 255]), colormap default
	end
	axis equal, axis off, 
	pause(.1)
	drawnow
end

% random search radii
radii = max(sz2p);
while radii(end) > 1, radii = [radii floor(radii(end)/2)]; end
nradii = length(radii);

% iterate
p_distindx = zeros(3*k, 2);
s_distindx = zeros(k + length(radii), 2);
for it = 1:2*iters,

	% set direction of propagation
	if mod(it,2) == 1, xrange = [1,w1p]; yrange = [1,h1p]; step =  1;
	else               xrange = [w1p,1]; yrange = [h1p,1]; step = -1;
	end

	if verbo && step > 0, disp(sprintf('% 2d  forward pass',floor((it-1)/2)+1)), end
	if verbo && step < 0, disp(sprintf('% 2d backward pass',floor((it-1)/2)+1)), end

	% patch averaging image
	if verbo && recon, 
		im12_0 = im12./repmat(ag12,[1 1 c1]); 
		im12 = zeros(h1, w1, c1); 
		ag12 = zeros(h1, w1,  1); 
	end

	for x1 = xrange(1):step:xrange(2),
	for y1 = yrange(1):step:yrange(2),

		P = im1(y1+pdom, x1+pdom, :);

		% propagation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		p_distindx(1:k, 1) = squeeze(nnf.dist(y1, x1, :));
		p_distindx(1:k, 2) = squeeze(nnf.indx(y1, x1, :));

		% propagate from x1 - step
		if (step > 0 && x1 > 1) || (step < 0 && x1 < w1p),
			indx = squeeze(nnf.indx(y1, x1-step, :));
			x2 = floor(indx / h2p)+1 + step;
			y2 =   mod(indx , h2p)+1;
			p_distindx(k+[1:k],2) = (x2-1)*h2p + y2-1;

			for i = 1:k,
				if x2(i) < 1 || x2(i) > w2p, p_distindx(k+i,1) = inf;
				else
					Q = im2(y2(i)+pdom, x2(i)+pdom, :);
					p_distindx(k+i,1) = mean((P(:) - Q(:)).^2);
				end
			end
		else p_distindx(k+[1:k],:) = inf;
		end

		% propagate from y2 - step
		if (step > 0 && y1 > 1) || (step < 0 && y1 < h1p),
			indx = squeeze(nnf.indx(y1-step, x1, :));
			x2 = floor(indx / h2p)+1;
			y2 =   mod(indx , h2p)+1 + step;
			p_distindx(2*k+[1:k],2) = (x2-1)*h2p + y2-1;

			for i = 1:k,
				if y2(i) < 1 || y2(i) > h2p, p_distindx(2*k+i,1) = inf;
				else
					Q = im2(y2(i)+pdom, x2(i)+pdom, :);
					p_distindx(2*k+i,1) = mean((P(:) - Q(:)).^2);
				end
			end
		else p_distindx(2*k+[1:k],:) = inf;
		end

		% keep k best propagated candidates
		p_distindx = unique(p_distindx,'rows');
		nnf.dist(y1, x1, :) = p_distindx(1:k, 1);
		nnf.indx(y1, x1, :) = p_distindx(1:k, 2);

		% random search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		s_distindx(1:k, 1) = squeeze(nnf.dist(y1, x1, :));
		s_distindx(1:k, 2) = squeeze(nnf.indx(y1, x1, :));
		x2c = floor(s_distindx(1, 2) / h2p)+1;
		y2c =   mod(s_distindx(1, 2) , h2p)+1;
		dc  = s_distindx(1, 1);
		for i = 1:nradii, r = radii(i);

			% draw uniform sample
			x2 = randi([max(1, x2c-r), min(w2p, x2c+r)]);
			y2 = randi([max(1, y2c-r), min(h2p, y2c+r)]);

			Q = im2(y2+pdom, x2+pdom, :);
			d = mean((P(:) - Q(:)).^2);
			s_distindx(k+i,1) = d;
			s_distindx(k+i,2) = (x2-1)*h2p + y2-1;

			% update search center
			if d < dc, dc  = d; x2c = x2; y2c = y2; end
		end

		% keep k best propagated candidates
		p_distindx = unique(s_distindx,'rows');
		nnf.dist(y1, x1, :) = s_distindx(1:k, 1);
		nnf.indx(y1, x1, :) = s_distindx(1:k, 2);

		if verbo && recon,
			x2 = floor(nnf.indx(y1, x1, recon) / h2p)+1;
			y2 =   mod(nnf.indx(y1, x1, recon) , h2p)+1;
			im12(y1+pdom, x1+pdom, :) = im12(y1+pdom, x1+pdom, :) + ...
			                            im2 (y2+pdom, x2+pdom, :);
			ag12(y1+pdom, x1+pdom) = ag12(y1+pdom, x1+pdom) + 1;
%			im12(y1, x1, :) = im12(y1, x1, :) + im2(y2, x2, :);
%			ag12(y1, x1) = ag12(y1, x1) + 1;
		end


	end
	% visual verbosity to an insane level
%	v = im12./repmat(ag12,[1 1 c1]);
%	idx = find(repmat(ag12, [1 1 c1]) == 0);
%	v(idx) = im12_0(idx);
%	imagesc(uint8(v), [0 255]), axis equal, axis off, colormap gray
%%	imagesc(sqrt(mean(nnf.dist,3)), [0 255]), axis equal, axis off, colormap default
%	pause(.01)
%	drawnow
	end

	ds(:,:,1+it) = nnf.dist(:,:,k);

	if verbo,
		if recon,
			imagesc(uint8(im12./repmat(ag12, [1 1 c1])), [0 255]), colormap gray
		else
			imagesc(sqrt(mean(nnf.dist,3)), [0 255]), colormap default
		end
		axis equal, axis off 
		pause(.1)
		drawnow
	end
end

