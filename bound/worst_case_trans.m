function C = wct(U, a)

U = (U > a);
[h,w] = size(U);

% patch match kernel is uniform
radii = fliplr(unique(round(max(h,w)./2.^[0:20])));
radii(find(radii < 1)) = [];


I = integralImage(U);

%p = zeros(h,w);
p = ones(h,w);
for r = radii,
	for i = 1:h,
	for j = 1:w,

		iM = min(h + 1,i + r); im = max(1,i - r);
		jM = min(w + 1,j + r); jm = max(1,j - r);

	%	p(i,j) = p(i,j) +  (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm); 
	%	p(i,j) = p(i,j) .* (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm); 
		pr(i,j) = (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm); 

	end
	end
	p = p .* pr; 
end

%C = max(p(:)/length(radii));
C = max(p(:));

