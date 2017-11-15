function [C,C2,p] = wct(U, Ua, prms, b)

if nargin == 2, b = inf; end

%Ua = (U > a);
[h,w] = size(Ua);

maxU = max(U(:));
minU = min(U(find(Ua)));

nvals = 21;
Ustep = (maxU - minU)/(nvals-1);
Uvals = [minU:Ustep:maxU];

I = integralImage(Ua);
Ia = zeros(h+1,w+1,nvals); % to compute acceptance probability
Ir = zeros(h+1,w+1,nvals); % to compute rejection  probability
for n = 1:nvals, u = Uvals(n);
	Ia(:,:,n) = integralImage(   Ua  .* (U < u));
	Ir(:,:,n) = integralImage((1-Ua) .* (U > u));
end

% patch match kernel is uniform
if prms.nradii == 1,
	radii = max(h,w);
else
	radii = fliplr(unique(round(max(h,w)./2.^[0:20])));
	radii(find(radii < 1)) = [];
end

p = ones(h,w); % uncomment to draw probability map
C = 0;
for i = 1:h,
for j = 1:w,
if U(i,j) >= minU && U(i,j) < b,

	% bin energy value
	n = floor((U(i,j)-minU)/Ustep) + 1;
	du = (U(i,j) - Uvals(n))/Ustep;

	pij = 1;
	if ~Ua(i,j),
		% the point is outside Ua -> a sample inside Ua must be accepted
		for r = radii,
			iM = min(h + 1,i + r); im = max(1,i - r);
			jM = min(w + 1,j + r); jm = max(1,j - r);

			% linear interpolation
			if n < nvals,
				tmp = (Ia(iM,jM,[n,n+1])+Ia(im,jm,[n,n+1])-Ia(im,jM,[n,n+1])-Ia(iM,jm,[n,n+1]))/(iM-im)/(jM-jm); 
				pij = pij * ( du * tmp(2) + (1 - du) * tmp(1) );
			else
				pij = pij * (Ia(iM,jM,n)+Ia(im,jm,n)-Ia(im,jM,n)-Ia(iM,jm,n))/(iM-im)/(jM-jm); 
			end
			if pij <= C, break, end % comment to draw probability map
		end
	else
		% the point is inside Ua -> a sample outside Ua must be rejected
		for r = radii,
			iM = min(h + 1,i + r); im = max(1,i - r);
			jM = min(w + 1,j + r); jm = max(1,j - r);

			% linear interpolation
			if n < nvals,
				tmp = (Ir(iM,jM,[n,n+1])+Ir(im,jm,[n,n+1])-Ir(im,jM,[n,n+1])-Ir(iM,jm,[n,n+1]))/(iM-im)/(jM-jm); 
				pr = ( du * tmp(2) + (1 - du) * tmp(1) );
			else
				pr = (Ir(iM,jM,n)+Ir(im,jm,n)-Ir(im,jM,n)-Ir(iM,jm,n))/(iM-im)/(jM-jm); 
			end
			pij = pij * (pr + (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm)); 
			if pij <= C, break, end % comment to draw probability map
		end
	end
	p(i,j) = pij; % uncomment to draw probability map
	C = max(pij,C);
else             % uncomment to draw probability map
	p(i,j) = 0;   % uncomment to draw probability map
end
end
end

C2 = mean(p(:));
