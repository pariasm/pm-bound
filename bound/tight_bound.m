function [C,ll,Uim] = compute_patchmatch_bound(a0, b, p, epsilons, prms)

shift = 1;

rsz = prms.rsz; % size of region
psz = prms.psz(1);  % size of patch

hpsz = floor(psz/2);

x = p(1)-hpsz;
y = p(2)-hpsz; % coordinates of a point

% crop a smaller image 
a = a0(max(1,x  - rsz + 1) : min(size(a0, 1), x + psz -1),...
       max(1,y  - rsz + 1) : min(size(a0, 2), y + psz -1));

sza = size(a);
szb = size(b);

A = im2col(a, [psz,psz], 'sliding'); % patches of A
B = im2col(b, [psz,psz], 'sliding'); % patches of B

% L2 distances between all patches in A and B
U = real(sqrt(bsxfun(@plus, sum(B.^2,1), bsxfun(@plus, sum(A.^2,1)', -2*A'*B))/psz/psz));
Umin = min(U,[],2);
U = bsxfun(@minus, U, Umin);

% energy as a szb x sza 4D image: Uim(:,:,i,j) is the energy map for ij in a
Uim = permute(reshape(U,[sza,szb]-psz+1),[3,4,1,2]);
% U = permute(reshape(Uim,[prod(sza-psz+1),prod(szb-psz+1)],[1,2,3,4]); % inverse


ll = nan*ones(size(Uim,3), size(Uim,4), length(epsilons));
for ieps = 1:length(epsilons), eps = epsilons(ieps);

	ll(end,end,ieps) = 0;

	y = size(Uim,4);
	for x = size(Uim,3)-1:-1:1,
		if shift,
			idx = find([(Uim(2:end,:,x+1,y) > eps - ll(x+1,y,ieps)); ones (1,size(Uim,2))]);
%			idx = find([(Uim(2:end,:,x+1,y) > eps - ll(x+1,y,ieps)); zeros(1,size(Uim,2))]);
		else
			idx = find(Uim(:,:,x+1,y) > eps - ll(x+1,y,ieps));
		end
		Uxy = Uim(:,:,x,y);
		ll(x,y,ieps) = eps - min(Uxy(idx));
	end

	x = size(Uim,3);
	for y = size(Uim,4)-1:-1:1,
		if shift,
			idx = find([(Uim(:,2:end,x,y+1) > eps - ll(x,y+1,ieps)); ones (size(Uim,1),1)]);
%			idx = find([(Uim(:,2:end,x,y+1) > eps - ll(x,y+1,ieps)), zeros(size(Uim,1),1)]);
		else
			idx = find(Uim(:,:,x,y+1) > eps - ll(x,y+1,ieps));
		end
		Uxy = Uim(:,:,x,y);
		ll(x,y,ieps) = eps - min(Uxy(idx));
	end

	for x = size(Uim,3)-1:-1:1,
	for y = size(Uim,4)-1:-1:1,

		if shift,
			idx = find([(Uim(2:end,:,x+1,y) > eps - ll(x+1,y,ieps)); ones(1,size(Uim,2))]);
		else
			idx = find(Uim(:,:,x+1,y) > eps - ll(x+1,y,ieps));
		end
		Uxy = Uim(:,:,x,y);
		if ~isempty(idx), ll_x = min(Uxy(idx));
		else              ll_x = inf;
		end

		if shift,
			idx = find([(Uim(:,2:end,x,y+1) > eps - ll(x,y+1,ieps)), ones(size(Uim,1),1)]);
		else
			idx = find(Uim(:,:,x,y+1) > eps - ll(x,y+1,ieps));
		end
		Uxy = Uim(:,:,x,y);
		if ~isempty(idx), ll_y = min(Uxy(idx));
		else              ll_y = inf;
		end

		ll(x,y,ieps) = eps - min(ll_x, ll_y);

	end
	end

end

C = ones(size(ll,1)*size(ll,2),length(epsilons));
%Cprod = ones(1,length(epsilons));

for ieps = 1:length(epsilons), eps = epsilons(ieps);

	disp(eps)

	% propagation neighborhood
	propidx = find(ll(:,:,ieps) < eps);
	[propx, propy] = ind2sub(size(ll(:,:,ieps)),propidx);
	nprop = length(propx);

	% computes C for points in the neighborhood
	for ip = 1:nprop,
		px = propx(ip);
		py = propy(ip);
		C(propidx(ip),ieps) = worst_case_trans(Uim(:,:,px,py) + ll(px,py,ieps), eps);
		if (C(propidx(ip),ieps) == 0),
			% no need to compute anymore
			break
		end
	end

	for iip = ip+1:nprop,
		C(propidx(iip),ieps) = 0;
	end
end

