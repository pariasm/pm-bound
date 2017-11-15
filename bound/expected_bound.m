function [C,ll,Uim] = compute_patchmatch_bound(a0, b, p, epsilons, prms)

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
U = sqrt(bsxfun(@plus, sum(B.^2,1), bsxfun(@plus, sum(A.^2,1)', -2*A'*B))/psz/psz);
Umin = min(U,[],2);
U = bsxfun(@minus, U, Umin);

% energy as a szb x sza 4D image: Uim(:,:,i,j) is the energy map for ij in a
Uim = permute(reshape(U,[sza,szb]-psz+1),[3,4,1,2]);
% U = permute(reshape(Uim,[prod(sza-psz+1),prod(szb-psz+1)],[1,2,3,4]); % inverse

% compute the energy similarity weights for adjacent pixel in a
wx = zeros(size(Uim,3), size(Uim,4));
wy = zeros(size(Uim,3), size(Uim,4));
for x = 1:size(Uim,3)-1,
for y = 1:size(Uim,4)-1,
	wx(x,y) = max(max(abs(Uim(2:end,:,x+1,y) - Uim(1:end-1,:,x,y))));
	wy(x,y) = max(max(abs(Uim(:,2:end,x,y+1) - Uim(:,1:end-1,x,y))));
end
end

y = size(Uim,4);
for k = 1:size(Uim,3)-1,
	wx(k,y) = max(max(abs(Uim(2:end,:,k+1,y) - Uim(1:end-1,:,k,y))));
end

x = size(Uim,3); 
for k = 1:size(Uim,4)-1,
	wy(x,k) = max(max(abs(Uim(:,2:end,x,k+1) - Uim(:,1:end-1,x,k))));
end

% run Dijstra's algorithm to computing geodesics
% ll contains the geodesic distance from every point in a to the bottom right
% pixel (because we are interested in the forward prop)
wx = flipud(fliplr(wx));
wy = flipud(fliplr(wy));

ll = inf*ones(size(Uim,3), size(Uim,4));
ll(:,1) = cumsum(wx(:,1));
ll(1,:) = cumsum(wy(1,:));
for x = 2:size(Uim,3),
for y = 2:size(Uim,4),
	ll(x,y) = min(ll(x-1,y) + wx(x,y), ll(x,y-1) + wy(x,y));
end
end
ll = flipud(fliplr(ll));


% compute the bound
%% epsilons = 0:2:100;

%% C = ones(sum(ll(:) < epsilons(end)),length(epsilons));
%% 
%% for ieps = 1:length(epsilons),
%% 	% energy level
%% 	eps = epsilons(ieps);
%% 
%% 	% propagation neighborhood
%% 	[propx, propy] = ind2sub([rsz,rsz],find(ll < eps));
%% 	nprop = length(propx);
%% 
%% 	% computes C for points in the neighborhood
%% 	for ip = 1:nprop,
%% 		px = propx(ip);
%% 		py = propy(ip);
%% 		C(ip,ieps) = worst_case_trans(Uim(:,:,px,py) - ll(px,py), eps);
%% 
%% 		imagesc(C,[0,1])
%% 		pause(.01)
%% 	end
%% end

C = ones(length(ll(:)),length(epsilons));
%Cprod = ones(1,length(epsilons));

for ieps = 1:length(epsilons),
	% energy level
	eps = epsilons(ieps);

	% propagation neighborhood
	propidx = find(ll < eps);
	[propx, propy] = ind2sub(size(ll),propidx);
	nprop = length(propx);

	% computes C for points in the neighborhood
	for ip = 1:nprop,
		px = propx(ip);
		py = propy(ip);
		C(propidx(ip),ieps) = expected_case_trans(Uim(:,:,px,py) + ll(px,py), eps);
%		Cprod(ieps) = Cprod(ieps)*C(propidx(ip),ieps);

%		imagesc(C,[0,1])
%		pause(.01)
	end
end

