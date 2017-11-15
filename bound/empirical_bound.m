function [P, V, sc0] = empirical_bound(u1, u2, prms, ...
                                  points, values, trials, sc0, P, V, prev_trials)

if ~exist('trials') || isempty(trials), trials = 100; end
if ~exist('P'), P = []; end

sz1 = size(u1);
sz2 = size(u2);

iters = prms.iters;
window = ones(prms.psz)/prod(prms.psz);

npoints = size(points,1);
nvalues = length(values);

if size(points,2) == 2,
	% points are given as [row/y, col/x] coordinates: compute indices
	points = (points(:,2)-1)*sz1(1) + points(:,1);
end

% masks
hpsz = floor(prms.psz(1)/2);
msk1 = ones(sz1); 
msk1(1:hpsz,:) = 0; msk1(end-hpsz+1:end,:) = 0;
msk1(:,1:hpsz) = 0; msk1(:,end-hpsz+1:end) = 0;

msk2 = ones(sz2);
msk2(1:hpsz,:) = 0; msk2(end-hpsz+1:end,:) = 0;
msk2(:,1:hpsz) = 0; msk2(:,end-hpsz+1:end) = 0;

% compute energy minima
if ~exist('sc0') || isempty(sc0),
	[nnfy0,nnfx0,sc0] = nnfield_exh(u1,msk1,u2,msk2,window,prms);
end

% do many trials and compute statistics
py = round(50/4); px = round(110/4);
scs = single(zeros(npoints, iters+1, prms.list));

if ~prev_trials || ~exist('P') || isempty(P),
	P = single(zeros(npoints, iters+1, nvalues, prms.list));
	V = single(zeros(npoints, iters+1, nvalues, prms.list));
	prev_trials = 0;
end

for i = prev_trials + [1:trials],

	% random initialization
	prms.mxit = 0;
	[nnfy,nnfx,sc] = nnfield(u1,msk1,u2,msk2,window,prms);

	for l = 1:prms.list,
		tmp = single(sc(:,:,l));
		scs(:,1,l) = tmp(points);
	end

	% fw-bw propagation passes
	prms.mxit = +1;
	for j = 1:iters
	 	% negative iters means to start with backward prop
		prms.mxit = -1*prms.mxit;
		[nnfy,nnfx,sc] = nnfield(u1,msk1,u2,msk2,window,prms,nnfy,nnfx);
		sc = sqrt(sc) - sqrt(sc0(:,:,1:prms.list));

		for l = 1:prms.list,
			tmp = single(sc(:,:,l));
			scs(:,j+1,l) = tmp(points);
		end
	end

	% update probabilities and variances
	for v = 1:length(values), 
		Pold = P(:,:,v,:);
		Pi = reshape(squeeze(scs > values(v)), size(P(:,:,v,:)));
		P(:,:,v,:) = 1/i*((i-1)*P(:,:,v,:) + Pi);
		V(:,:,v,:) = 1/i*((i-1)*V(:,:,v,:) + (Pi - Pold).*(Pi - P(:,:,v,:)));
	end

	plot([0:iters],squeeze(P(end,:,:,2)),'.-'), hold on
	plot([0:iters],squeeze(P(end,:,:,2) + 3*sqrt(V(end,:,:,2)/i)),':'),
	plot([0:iters],squeeze(P(end,:,:,2) - 3*sqrt(V(end,:,:,2)/i)),':'), hold off
	ylim([-.1, 1.1]), xlim([0,iters]), grid on
	title(sprintf('%d trials', i)), drawnow
	
end

