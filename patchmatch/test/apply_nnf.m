function u2 = apply_nnf(u1,msk1,nx_vol,ny_vol,slice)

if nargin == 4,
	slice = 1;
end

nx = nx_vol(:,:,slice);
ny = ny_vol(:,:,slice);

sz1 = size(u1);

idx1 = find(msk1 == 1);
[idx1y,idx1x] = ind2sub(sz1,idx1);

%[xx yy] = meshgrid([1:sz1(2)],[1:sz1(1)]);

idx1x = min(sz1(2), max(1, idx1x + nx(idx1) ));
idx1y = min(sz1(1), max(1, idx1y + ny(idx1) ));

idx2 = sub2ind(sz1,idx1y,idx1x);

u2 = u1;
u2(idx1) = u2(idx2);

