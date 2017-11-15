function C = wct(U, a, prms)

if nargin == 2,
	prms.nradii = 0;
	prms.list = 1;
end

U = (U >= a);
[h,w] = size(U);

if prms.nradii == 1,
	radii = max(h,w);
else
	radii = fliplr(unique(round(max(h,w)./2.^[0:20])));
	radii(find(radii < 1)) = [];
end
nradii = length(radii);


I = integralImage(U);

%% p = ones(h,w);
%% for r = radii,
%% 	for i = 1:h,
%% 	for j = 1:w,
%% 
%% 		iM = min(h + 1,i + r); im = max(1,i - r);
%% 		jM = min(w + 1,j + r); jm = max(1,j - r);
%% 
%% 	%	p(i,j) = p(i,j) +  (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm); 
%% 	%	p(i,j) = p(i,j) .* (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm); 
%% 		pr(i,j) = (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm); 
%% 
%% 	end
%% 	end
%% 	p = p .* pr; 
%% end


pradii = ones(h,w,nradii);
for i = 1:h,
for j = 1:w,

	%% real patchmatch bound!
%	pradii = ones(nradii,1);
	for ir = 1:nradii, r = radii(ir);
		iM = min(h + 1,i + r); im = max(1,i - r);
		jM = min(w + 1,j + r); jm = max(1,j - r);
		pradii(i,j,ir) = (I(iM,jM)+I(im,jm)-I(im,jM)-I(iM,jm))/(iM-im)/(jM-jm); 
	end

%	pr = sort(pradii, 1, 'descend');
%	p(i,j) = prod(pr(1:nradii - prms.list + 1)); 

end
end

% combinatorics
n = nradii;
inmin = n - prms.list + 1;
p = prod(pradii,3);
if inmin < n, for i0 = 1:n,
  idxi0 = 1:n; idxi0(find(idxi0 == i0)) = [];
  idxo0 = [];  idxo0 = [idxo0 i0];
  p = p + prod(pradii(:,:,idxi0),3).*prod(1-pradii(:,:,idxo0),3);
  if inmin < n-1, for i1 = i0+1:n,
    idxi1 = idxi0; idxi1(find(idxi1 == i1)) = [];
    idxo1 = idxo0; idxo1 = [idxo1 i1];
    p = p + prod(pradii(:,:,idxi1),3).*prod(1-pradii(:,:,idxo1),3);
    if inmin < n-2, for i2 = i1+1:n,
      idxi2 = idxi1; idxi2(find(idxi2 == i2)) = [];
      idxo2 = idxo1; idxo2 = [idxo2 i2];
      p = p + prod(pradii(:,:,idxi2),3).*prod(1-pradii(:,:,idxo2),3);
      if inmin < n-3, for i3 = i2+1:n,
        idxi3 = idxi2; idxi3(find(idxi3 == i3)) = [];
        idxo3 = idxo2; idxo3 = [idxo3 i3];
        p = p + prod(pradii(:,:,idxi3),3).*prod(1-pradii(:,:,idxo3),3);
        if inmin < n-4, for i4 = i3+1:n,
          idxi4 = idxi3; idxi4(find(idxi4 == i4)) = [];
          idxo4 = idxo3; idxo4 = [idxo4 i4];
          p = p + prod(pradii(:,:,idxi4),3).*prod(1-pradii(:,:,idxo4),3);
          if inmin < n-5, for i5 = i4+1:n,
            idxi5 = idxi4; idxi5(find(idxi5 == i5)) = [];
            idxo5 = idxo4; idxo5 = [idxo5 i5];
            p = p + prod(pradii(:,:,idxi5),3).*prod(1-pradii(:,:,idxo5),3);
            if inmin < n-6, for i6 = i5+1:n,
              idxi6 = idxi5; idxi6(find(idxi6 == i6)) = [];
              idxo6 = idxo5; idxo6 = [idxo6 i6];
              p = p + prod(pradii(:,:,idxi6),3).*prod(1-pradii(:,:,idxo6),3);
              if inmin < n-7, for i7 = i6+1:n,
                idxi7 = idxi6; idxi7(find(idxi7 == i7)) = [];
                idxo7 = idxo6; idxo7 = [idxo7 i7];
                p = p + prod(pradii(:,:,idxi7),3).*prod(1-pradii(:,:,idxo7),3);
                if inmin < n-8, for i8 = i7+1:n,
                  idxi8 = idxi7; idxi8(find(idxi8 == i8)) = [];
                  idxo8 = idxo7; idxo8 = [idxo8 i8];
                  p = p + prod(pradii(:,:,idxi8),3).*prod(1-pradii(:,:,idxo8),3);
                  if inmin < n-9, for i9 = i8+1:n,
                    idxi9 = idxi8; idxi9(find(idxi9 == i9)) = [];
                    idxo9 = idxo8; idxo9 = [idxo9 i9];
                    p = p + prod(pradii(:,:,idxi9),3).*prod(1-pradii(:,:,idxo9),3);
                  end, end %9
                end, end %8
              end, end %7
            end, end %6
          end, end %5
        end, end %4
      end, end %3
    end, end %2
  end, end %1
end, end %0

%C = max(p(:)/length(radii));
C = max(p(:));


