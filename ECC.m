% Implement ECDH
% elliptic curve coefficients
p = 17;			% finite field GF(p) = Fp
N = 17;			% order of curve; number of points on the curve
n = 17;			% order of subgroup; number of points in subgroup start from G
h = N/n;		% curve cofactor
a = 0;			% curve coefficient
b = 7;			% curve coefficient
G = [15,13];	% start point

% generate finite field 2-D map
[X, Y] = meshgrid((0:p-1), (0:p-1));

% elliptic curve equation
field = mod(Y.^2 - X.^3 - a.*X - b, p);

% list valid points on the curve in the finite field
ECpoints = [X(field == 0) Y(field == 0)];

figure
hold on;
plot(ECpoints(:, 1), ECpoints(:, 2), 'o');
hold on;
ax = gca;
ax.XTick = 0:17;
ax.YTick = 0:17;
axis([0 17 0 17]);
grid on;

% Crypto
disp('cyclic subgroup elements: ------------------------------')
linecolor_base = [0, 0.1,  1];
linecolor_grad = [1, 0, -1] ./ n;

P = G;
xP = G(1);
yP = G(2);

% First, derive R=2P
[gcd, x, y] = ExtEuclAlgo(2*yP, p);
mult_inv = mod(x, p);
s = mod((3 * xP^2 + a) * mult_inv, p);	% diff with P + Q

xR = mod(s^2 - 2 * xP, p);
yR = mod(s * (xP - xR) - yP, p);
R = [xR, yR];

% draw line from P to R=2P
line = [P; R];
plot(line(:,1), line(:,2), 'color', linecolor_base + linecolor_grad * 1)
drawnow

% iteratively traverse each points in subgroup
for r = 2:n+1
	
	% update Q
	Q = R;
	xQ = R(1);
	yQ = R(2);
	
	% calculate P + Q = R i.e. P + rP = (r+1)P
	if xQ ~= xP
		[gcd, x, y] = ExtEuclAlgo((xQ - xP), p);	% find the mult inverse
		mult_inv = mod(x, p);						% of (xQ - xP)
		s = mod((yQ - yP) * mult_inv, p);
		xR = mod(s^2 - xP - xQ, p);
		yR = mod(s * (xP - xR) - yP, p);
		R = [xR, yR]
	else
		R = [G(1), 0]
	end
	
	% draw line from Q (rP) to R ((r+1)P)
	line = [Q; R]; 
	plot(line(:,1), line(:,2), 'color', linecolor_base + linecolor_grad * r)
	drawnow
	pause(0.1);
end
%%
Q = ECpoint_Scale(G, 6, p)
%%
x = 5;
mod(x^(p-2), p)

%% 
P1 = [15, 13]
P2 = [15, 4]
% Q = ECpoint_Double([P1(1), P1(2), 1]);
Q = ECpoint_Addition([P1(1), P1(2), 1], [P2(1), P2(2), 1])
[gcd, x, y] = ExtEuclAlgo(Q(3), p);	% find the mult inverse
Z_inv = mod(x, p)
P2 = [mod(Q(1) * Z_inv^2, p), mod(Q(2) * Z_inv^3, p)]

% need deal with the infinite point O.
% how to make O + P = P?

%%
function [Q] = ECpoint_Scale(P, k, p)
	% to jacobian coordinate
	Pjac = [P(1), P(2), 1];
	Qjac = [0, 0, 0];
	l = 8;
	kbit = dec2bin(int16(k), l);
	for j = (0:l-1)
		Qjac = ECpoint_Double(Qjac);
		if kbit(j + 1) == '1'
			if (Qjac == [0, 0, 0])
				Qjac = Pjac;
			else
				Qjac = ECpoint_Addition(Pjac, Qjac);
			end
		end
	end
	[gcd, x, y] = ExtEuclAlgo(Qjac(3), p);	% find the mult inverse
	Z_inv = mod(x, p);
	Q = [mod(Qjac(1) * Z_inv^2, p), mod(Qjac(2) * Z_inv^3, p)];
end

function [Q] = ECpoint_Double(P)
	a = 0; p = 17;
	X1 = P(1); Y1 = P(2); Z1 = P(3);
	M = mod(3 * X1^2 + a * Z1^4, p);
	S = mod(4 * X1 * Y1^2, p);
	T = mod(8 * Y1^4, p);
	X2 = mod(M^2 - 2*S, p);
	Y2 = mod(M * (S-X2) - T, p);
	Z2 = mod(2 * Y1 * Z1, p);
	Q = [X2, Y2, Z2];
end

function [Q] = ECpoint_Addition(P0, P1)
	p = 17;
	X0 = P0(1);
	Y0 = P0(2); 
	Z0 = P0(3);	
	X1 = P1(1); Y1 = P1(2); Z1 = P1(3);
	M = mod(Y0 * Z1^3 + Y1 * Z0^3, p);
	R = mod(Y0 * Z1^3 - Y1 * Z0^3, p);
	T = mod(X0 * Z1^2 + X1 * Z0^2, p);
	W = mod(X0 * Z1^2 - X1 * Z0^2, p);
	X2 = mod(R^2 - T * W^2, p);
	V = mod(T * W^2 - 2 * X2, p);
	Y2 = mod((V * R - M * W^3) * 9, p);	% <-- 2Y2?
	Z2 = mod(Z0 * Z1 * W, p);
	Q = [X2, Y2, Z2];
end

%%
function [gcd, s_, t_] = ExtEuclAlgo(a, b)
% Description
%	  Extended Euclidean Algorithm
%     Returns a three-tuple (gcd, x, y) such that
%     a * x + b * y == gcd, where gcd is the greatest
%     common divisor of a and b.
% 
%     This function implements the extended Euclidean
%     algorithm and runs in O(log b) in the worst case.

	% [now, old] 
	s = [0, 1];
	t = [1, 0];
	r = [b, a];

	while r(1) ~= 0
		quotient = floor(r(2) / r(1));	% get quotient
		r = [r(2) - quotient * r(1), r(1)];
		s = [s(2) - quotient * s(1), s(1)];
		t = [t(2) - quotient * t(1), t(1)];
	end
	
	gcd = r(2);
	s_ = s(2);
	t_ = t(2);
end

%% 
% https://cryptobook.nakov.com/asymmetric-key-ciphers/elliptic-curve-cryptography-ecc
% https://zh.wikipedia.org/wiki/%E6%A4%AD%E5%9C%86%E6%9B%B2%E7%BA%BF

