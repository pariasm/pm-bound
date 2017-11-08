function G = gauss_kernel(sz,a);

hsz = floor(sz/2);

[Gx,Gy] = meshgrid([-hsz(1):hsz(1)],[-hsz(1):hsz(1)]);

G = exp(-(Gx.^2 + Gy.^2)/a);
G = G/sum(G(:));

