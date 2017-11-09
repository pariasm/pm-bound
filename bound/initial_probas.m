function Pi = initial_probas(a, b, p, epsilons, prms, sc0)

psz = prms.psz(1); % size of patch
hpsz = floor(psz/2);

x = p(1)-hpsz;
y = p(2)-hpsz; % coordinates of a point

% crop a smaller image 
A = a(x + [0:psz-1], y + [0:psz-1]); A = A(:);
B = im2col(b, [psz,psz], 'sliding'); % patches of B

szb = size(b);

% L2 distances between all patches in A and B
U = real(sqrt(bsxfun(@plus, sum(B.^2,1), bsxfun(@plus, sum(A.^2,1)', -2*A'*B))/psz/psz));
Umin = min(U,[],2);
U = bsxfun(@minus, U, Umin);

% energy as a szb image
Uim = reshape(U,szb-psz+1);

Pi = zeros(length(epsilons),1);
for ieps = 1:length(epsilons), eps = epsilons(ieps);
	Pi(ieps) = sum(sum(Uim(:,:) > eps))/prod(szb - psz + 1);
end


