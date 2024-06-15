function [x] = PCG(Amap,b,D,toler)
%Solve for Ax = b, A is symmetric and PD, until ||Ap-b||<= toler
%D is the preconditioner

if (size(D,1) == 1)
	D = D';
end
if (isscalar(D))
	Dmap = @(x) D*x;
elseif (size(D,2) == 1)
	Dmap = @(x) D.*x;
else
	Dmap = @(x) D*x;
end

if (~isa(Amap,'function_handle'))
	Amap = @(x) Amap*x;
end
x = zeros(size(b));
r = b;
z = Dmap(r);
p = z;
nr = norm(r);
if nr < toler
	j = 0;
	return;
end

rz = r' * z;
for j = 1:length(x)
	Ap = Amap(p);
	kap = p' * Ap;
	alpha = rz / kap;
	x = x + alpha * p;
	r = r - alpha * Ap;
	nrplus = norm(r);
	if nrplus < toler
		return;
	end
	z = Dmap(r);
	newrz = r' * z;
	beta = newrz / rz;
	rz = newrz;
	p = z + beta * p;
end
end

