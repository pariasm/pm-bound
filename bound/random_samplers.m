
% the random samplers
m = 100;

% the sampling domain
n = 100*100;
u2 = ones(n,1);
u2(randi(n)) = 0; % we look for this position

iters = 251;

trials = 1000;
P = zeros(iters,1);
for t = 1:trials,
	U = 1;
	Us = zeros(iters,1);
	for i = 1:iters,
		U = min(U, min(u2(randi(n,m,1))));
		Us(i) = U;
		if U == 0,
			break;
		end
	end
	P = ((t-1)*P + Us)/t;
	plot(1:iters,P)
end
