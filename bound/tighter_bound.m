function [C,C2,Cold,Pi,ll,llold,Uim] = compute_patchmatch_bound(a0, b, p, epsilons, prms, mepsilons)

if nargin == 5,
	mepsilons = inf*ones(size(epsilons));
end

shift = 1;

rsz = prms.rsz;    % size of region
psz = prms.psz(1); % size of patch

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computed old version of the bound (bound.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

llold = inf*ones(size(Uim,3), size(Uim,4));
llold(:,1) = cumsum(wx(:,1));
llold(1,:) = cumsum(wy(1,:));
for x = 2:size(Uim,3),
for y = 2:size(Uim,4),
	llold(x,y) = min(llold(x-1,y) + wx(x,y), llold(x,y-1) + wy(x,y));
end
end
llold = flipud(fliplr(llold));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computed old version of the bound (bound.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









ll = nan*ones(size(Uim,3), size(Uim,4), length(epsilons));
C = ones(size(ll));
Cold = ones(size(ll));
C2 = ones(size(ll));
Pi = zeros(length(epsilons),1);
for ieps = 1:length(epsilons), eps = epsilons(ieps);

%	disp(eps)

	MM = zeros(size(Uim));
	ll(end,end,ieps) = eps;
	MM(:,:,end,end) = Uim(:,:,end,end) > ll(end,end,ieps);

	y = size(Uim,4);
	for x = size(Uim,3)-1:-1:1,
		if shift,
			MM(:,:,x,y) = [Uim(2:end,:,x+1,y) >= ll(x+1,y,ieps);  ones(1,size(Uim,2))];
%			MM(:,:,x,y) = [Uim(2:end,:,x+1,y) >= ll(x+1,y,ieps); zeros(1,size(Uim,2))];
		else
			MM(:,:,x,y) = Uim(:,:,x+1,y) >= ll(x+1,y,ieps);
		end
		Uxy = Uim(:,:,x,y);
		idx = find(MM(:,:,x,y));
		if ~isempty(idx), ll(x,y,ieps) = min(Uxy(idx));
		else              ll(x,y,ieps) = inf;
		end
	end

	x = size(Uim,3);
	for y = size(Uim,4)-1:-1:1,
		if shift,
			MM(:,:,x,y) = [(Uim(:,2:end,x,y+1) >= ll(x,y+1,ieps)),  ones(size(Uim,1),1)];
%			MM(:,:,x,y) = [(Uim(:,2:end,x,y+1) >= ll(x,y+1,ieps)), zeros(size(Uim,1),1)];
		else
			MM(:,:,x,y) = Uim(:,:,x,y+1) >= ll(x,y+1,ieps);
		end
		Uxy = Uim(:,:,x,y);
		idx = find(MM(:,:,x,y));
		if ~isempty(idx), ll(x,y,ieps) = min(Uxy(idx));
		else              ll(x,y,ieps) = inf;
		end
	end

	for x = size(Uim,3)-1:-1:1,
	for y = size(Uim,4)-1:-1:1,
		if shift,
			MMx = [(Uim(2:end,:,x+1,y) >= ll(x+1,y,ieps)); ones(1,size(Uim,2))];
			MMy = [(Uim(:,2:end,x,y+1) >= ll(x,y+1,ieps)), ones(size(Uim,1),1)];
		else
			MMx = (Uim(:,:,x+1,y) >= ll(x+1,y,ieps))
			MMy = (Uim(:,:,x,y+1) >= ll(x,y+1,ieps));
		end
		MM(:,:,x,y) = MMx .* MMy; % intersection
		Uxy = Uim(:,:,x,y);

		if prms.list == 1,
			idx = find(MM(:,:,x,y));
			if ~isempty(idx), ll(x,y,ieps) = min(Uxy(idx));
			else              ll(x,y,ieps) = inf;
			end
		else
			idx = find(MMx);
			if ~isempty(idx), llx = min(Uxy(idx));
			else              llx = inf;
			end

			idx = find(MMy);
			if ~isempty(idx), lly = min(Uxy(idx));
			else              lly = inf;
			end

			% max of the mins
			ll(x,y,ieps) = max(llx,lly);
		end
	end
	end

	% compute initial probabilities
	Pi(ieps) = sum(sum(MM(:,:,end,end)))/prod(szb - psz + 1);

	disp('Computing C');

	% scan the propagation neighborhood starting from the bottom-right corner
	szll = [size(ll,1),size(ll,2)];
	for d = 0:2*max(szll),  % d is the diagonal (0 means bottom-right corner)
		for dd = -d/2:d/2, % dd is the displacement along the diagonal
			px = szll(1) - d/2 +  dd;
			py = szll(2) - d/2 + -dd;
			if px < 1 || py < 1, continue, end

			if strcmp(prms.transition_kernel, 'acceptance'),
				[C(px,py,ieps), C2(px,py,ieps)] = better_case_trans(Uim(:,:,px,py), MM(:,:,px,py), prms, mepsilons(ieps));
                Cold(px,py,ieps) =  worst_case_trans(Uim(:,:,px,py), eps - llold(px,py), prms);
			else
				C(px,py,ieps) =  worst_case_trans(Uim(:,:,px,py), ll(px,py,ieps), prms);
				Cold(px,py,ieps) =  worst_case_trans(Uim(:,:,px,py), eps - llold(px,py), prms);
			end

			imagesc(C(:,:,ieps)); drawnow
			disp(sprintf('%3d,%3d - %6.4f - %g', px, py, C(px,py,ieps), prod(prod(C(:,:,ieps)))))

			if (prod(prod(C(:,:,ieps))) <= 1e-7),
				% no need to compute anymore
				break
			end
		end
		if (prod(prod(C(:,:,ieps))) <= 1e-7),
			% no need to compute anymore
			break
		end
	end

end


