# ECDH implementation
import numpy as np
import time
import sys

##### Curve parameters -------------------------------------------------------------
Keyclass = 0
if Keyclass == 1:
	# Edu curve
	p = 17;			# finite field GF(p) = Fp
	N = 17;			# order of curve; number of points on the curve
	n = 17;			# order of subgroup; number of points in subgroup start from G
	h = N/n;			# curve cofactor
	a = 0;			# curve coefficient
	b = 7;			# curve coefficient
	G = [15,13];	# start point
	L = 5			# bit length
else:
	# secp192r1 
	p = int(0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff)
	n = int(0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551)
	h = 1
	a = int(0xffffffff00000001000000000000000000000000fffffffffffffffffffffffc)
	b = int(0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b)
	G = [int(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296), \
		 int(0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)]
	L = 256
##### ------------------------------------------------------------------------------

max_iter = 0
outputfile = 1

if outputfile:
    f_base = open('./TP/basepoint.dat', 'w')
    f_end = open('./TP/endpoint.dat', 'w')
    f_scalar = open('./TP/scalar.dat', 'w')
    data_file = open("./TP/reg_record.dat", 'w')
    f_pointadd = open("./TP/pointadd.dat", 'w')
    f_pointdbl = open("./TP/pointdbl.dat", 'w')

def ECpoint_Scale(P, k):
	# Q = [k]P
	# to jacobian coordinate
	global p
	global L
	Pjac = [P[0], P[1], 1];
	Qjac = [0, 0, 0];
	kbitarray = list(np.binary_repr(k, width=L))
								# '1', '0', '1', '1' (index: 0 ~ 3)
	kbitarray.reverse()			# '1', '1', '0', '1' (index: 0 ~ 3)
	for j in reversed(range(L)):		# for L-1, L-2, ..., 1, 0
		Qjac = ECpoint_Doubling(Qjac);	# double Q: Q := 2Q
		# print(f'{j}: k[j]={kbitarray[j]}  {Qjac}')
		if int(kbitarray[j]) == 1:
			if (Qjac == [0, 0, 0]):		# first visit
				Qjac = Pjac;
			else:	
				# add P: Q := Q + P
				Qjac = ECpoint_Addition(Qjac, Pjac);

	# print(f'end {Qjac}')
	# transform back to normal domain
	Z_inv = MultInv(Qjac[2])	# cal: 1/Z mod p
    # [x y] = [X / Z^2,  Y / Z^3]
	Q = [mod(Qjac[0] * Z_inv**2, p), mod(Qjac[1] * Z_inv**3, p)];

	return Q

def ECpoint_Doubling(P):
	# Q = 2P
	global p
	global a
	[X1, Y1, Z1] = P;
	# print(P)
	M = mod(3 * X1**2 + a * Z1**4, p);
	S = mod(4 * X1 * Y1**2, p);
	T = mod(8 * Y1**4, p);
	X2 = mod(M**2 - 2*S, p);
	Y2 = mod(M * (S-X2) - T, p);
	Z2 = mod(2 * Y1 * Z1, p);
	Q = [X2, Y2, Z2];
	# print(Q)
	return Q

def ECpoint_Addition(P0, P1):
	# Q = P0 + P1
	global p
	[X0, Y0, Z0] = P0;	
	[X1, Y1, Z1] = P1;
	M = mod(Y0 * Z1**3 + Y1 * Z0**3, p);
	R = mod(Y0 * Z1**3 - Y1 * Z0**3, p);
	T = mod(X0 * Z1**2 + X1 * Z0**2, p);
	W = mod(X0 * Z1**2 - X1 * Z0**2, p);
	X2 = mod(R**2 - T * W**2, p);
	V = mod(T * W**2 - 2 * X2, p);
	Y2 = mod((V * R - M * W**3) * (p//2 + 1), p);
	Z2 = mod(Z0 * Z1 * W, p);
	Q = [X2, Y2, Z2];
	# print(Q)
	return Q

def mod(a, p):
	# b = a mod p
	b = a % p
	return b

def ECpoint_Scale_HDL(P, k):
	# Q = [k]P
	# to jacobian coordinate
	global p
	global L
	[X0, Y0, Z0] = [0, 0, 0];			# current point	P
	[X1, Y1, Z1] = [P[0], P[1], 1];		# base point	Q

    # only need to initialize G, H, I
	[A, B, C, D, E, F, G, H, I] = [0, 0, 0, 0, 0, 0, 0, 0, 0] # [X0, Y0, Z0, X0^2, 0, 0, X1, Y1, 1]
	is_first = 1

	kbitarray = list(np.binary_repr(k, width=L))
	kbitarray.reverse()
	for j in reversed(range(L)):		# for L-1, L-2, ..., 1, 0
		if (not is_first):
			[A, B, C, D, E, F, G, H, I] = \
				ECpoint_Doubling_HDL(A, B, C, D, E, F, G, H, I)	# double Q: Q := 2Q
			data_file.write(output256bits(A))
			data_file.write(output256bits(B))
			data_file.write(output256bits(C))
			data_file.write(output256bits(D))
		if int(kbitarray[j]) == 1:
			if (is_first == 1):		# first visit
				is_first = 0
				A, B, C = X1, Y1, 1
				D = mod256(A * A, p)	# prepare X0^2
			else:	
				# add P: Q := Q + P
				G, H, I = X1, Y1, 1		# refresh
				[A, B, C, D, E, F, G, H, I] = \
					ECpoint_Addition_HDL(A, B, C, D, E, F, G, H, I)
				data_file.write(output256bits(A))
				data_file.write(output256bits(B))
				data_file.write(output256bits(C))
				data_file.write(output256bits(D))
		else:
			if (not is_first):
				D = mod256(A * A, p)
	
	# transform back to normal domain
	Z_inv = MultInv(C)	# cal: 1/Z mod p
    # [x y] = [X / Z^2,  Y / Z^3]
	Z_inv2 = mod256(Z_inv * Z_inv, p)
	Z_inv3 = mod256(Z_inv2 * Z_inv, p)
	# print('Z^(-2): ' + output256bits(Z_inv2))
	# print('Z^(-3): ' + output256bits(Z_inv3))
	Q = [mod256(A * Z_inv2, p), mod256(B * Z_inv3, p)];
	# print('Qx: ' + output256bits(Q[0]))
	# print('Qy: ' + output256bits(Q[1]))

	return Q

def ECpoint_Doubling_HDL(A, B, C, D, E, F, G, H, I):
	# input: A = X1, B = Y1, C = Z1, D = X1^2
	global p
	global a
	# print('EC point DBL')
	# print([A, B, C])
	E = mod256(3 * D, p)	# 3 * X1^2
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t1 ---------------------------
	D = mod256(C * C, p)	# Z1^2
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t2 ---------------------------
	D = mod256(D * D, p)	# Z1^4
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t3 ---------------------------
	D = mod256(a * D, p)	# a * Z1^4
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t4 ---------------------------
	E = mod256(D + E, p)	# M = 3X1^2+aZ1^4
	D = mod256(B * C, p)	# Y1*Z1
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t5 ---------------------------
	B = mod256(B * B, p)	# Y1^2
	C = mod256(D + D, p)	# Z2 $$
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t6 ---------------------------
	D = mod256(B * A, p)	# X1*Y1^2
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t7 ---------------------------
	D = mod256(4 * D, p)	# S = 4X1Y1^2
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t8 ---------------------------
	A = mod256(D + D, p)	# 2S
	F = mod256(E * E, p)	# M^2
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t9 ---------------------------
	B = mod256(B * B, p)	# Y1^4
	A = mod256(F + p - A, p)	# X2 $$
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t10 --------------------------
	B = mod256(8 * B, p)	# T
	D = mod256(D + p - A, p)	# S - X2
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t11 --------------------------
	D = mod256(D * E, p)	# M(S-X2)
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t12 --------------------------
	B = mod256(D + p - B, p)	# Y2
	D = mod256(C * C, p)	# Z2^2
	outputDBL_Regfile(A, B, C, D, E, F, G, H, I)
	# t13 --------------------------
	# output: A = X2, B = Y2, C = Z2, D = Z2^2
	# did not use G, H, I
	return [A, B, C, D, E, F, G, H, I]

def ECpoint_Addition_HDL(A, B, C, D, E, F, G, H, I):
	# input: A = X0, B = Y0, C = Z0, D = Z0^2, 
	# 		 G = X1, H = Y1, I = Z1
	# P0 = Q_cur, P1 = base point
	global p
	# print('EC point ADD')
	E = mod256(G * D, p)	# U1 = X1*Z0^2
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t01 ---------------------------
	D = mod256(D * C, p)	# Z0^3
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t02 ---------------------------
	I = mod256(A + p - E, p)	# W = U0 - U1
	D = mod256(D * H, p)	# S1 = Y1 * Z0^3
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t03 ---------------------------
	F = mod256(B + p - D, p)	# R = S0 - S1
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t04 ---------------------------
	G = mod256(F * F, p)	# R^2
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t05 ---------------------------
	A = mod256(A + E, p)	# T = U0 + U1
	E = mod256(I * I, p)	# W^2
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t06 ---------------------------
	D = mod256(D + B, p)	# M
	B = mod256(A * E, p)	# TW^2
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t07 ---------------------------
	A = mod256(G + p - B, p)	# X2 $$
	E = mod256(I * E, p)	# W^3
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t08 ---------------------------
	G = mod256(E * D, p)	# MW^3
	D = mod256(A + A, p)	# 2X2
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t09 ---------------------------
	B = mod256(B + p - D, p)	# V
	# print(f'V = {B}')
	C = mod256(C * I, p)	# Z2 $$
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t10 ---------------------------
	F = mod256(F * B, p)	# VR
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t11 ---------------------------
	D = mod256(A * A, p)	# X2^2
	F = mod256(F + p - G, p)	# 2Y2
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t12 ---------------------------
	B = mod256(F * (p//2 + 1), p)	# Y2 $$
	outputADD_Regfile(A, B, C, D, E, F, G, H, I)
	# t13 ---------------------------
	# output: A = X2, B = Y2, C = Z2, D = X2^2
	# remenber to fresh G and H with P1
	return [A, B, C, D, E, F, G, H, I]

def mod256(a, p):
	# b = a mod p
	# b = a % p
	# return b
	if a >= (p*p):
		print('\nout of bound')
		print(output256bits(a))
		print('\n')
		sys.exit()
	elif a < 0:
		print('\na < 0\n')
        
	A = [0]*16
	for i in range(16):
		A[i] = (a >> (32*i)) % (2**32)	# A[i] = a[i*64 +: 64]
	T  = (A[ 7] << 32*7) + (A[ 6] << 32*6) + (A[ 5] << 32*5) + (A[ 4] << 32*4) + (A[ 3] << 32*3) + (A[ 2] << 32*2) + (A[ 1] << 32*1) + A[ 0]	
	S1 = (A[15] << 32*7) + (A[14] << 32*6) + (A[13] << 32*5) + (A[12] << 32*4) + (A[11] << 32*3)
	S2 = 				   (A[15] << 32*6) + (A[14] << 32*5) + (A[13] << 32*4) + (A[12] << 32*3)
	S3 = (A[15] << 32*7) + (A[14] << 32*6) 														 + (A[10] << 32*2) + (A[ 9] << 32*1) + A[ 8]
	S4 = (A[ 8] << 32*7) + (A[13] << 32*6) + (A[15] << 32*5) + (A[14] << 32*4) + (A[13] << 32*3) + (A[11] << 32*2) + (A[10] << 32*1) + A[ 9]
	D1 = (A[10] << 32*7) + (A[ 8] << 32*6) 									   					 + (A[13] << 32*2) + (A[12] << 32*1) + A[11]
	D2 = (A[11] << 32*7) + (A[ 9] << 32*6) 									   + (A[15] << 32*3) + (A[14] << 32*2) + (A[13] << 32*1) + A[12]
	D3 = (A[12] << 32*7) 				   + (A[10] << 32*5) + (A[ 9] << 32*4) + (A[ 8] << 32*3) + (A[15] << 32*2) + (A[14] << 32*1) + A[13]
	D4 = (A[13] << 32*7)				   + (A[11] << 32*5) + (A[10] << 32*4) + (A[ 9] << 32*3)				   + (A[15] << 32*1) + A[14]

	bs = [0]*4
	bs[0] = T - D1
	bs[1] = 2 * S1 - D2 + p
	bs[2] = 2 * S2 - D3 + p
	bs[3] = S3 + S4 - D4
    
	b1 = bs[0] + bs[1] 
	b2 = bs[2] + bs[3]
	if (b1 < 0):
		b1 = b1 + p
	elif b1 > p:
		b1 = b1 - p
		
	if (b2 < 0):
		b2 = b2 + p
	elif b2 > p:
		b2 = b2 - p

	b = b1 + b2

	# exception detection
	if (b < 0):
		print(f'mod256: b < 0: {b:X}')
	quotient = b // p		# quotient is usually 0, 1, 2, 3
	if (quotient > 4):
		print(f'mod256: b > 4p, q={quotient} !!!!!!!!!!!!!!!!!!!!!!!!!!!')
    
	if (b < p):
		return b
	elif (b - p) < p:
		return b - p
	elif (b - 2*p) < p:
		return b - 2*p
	elif (b - 3*p) < p:
		return b - 3*p
	else:
		return b - 4*p

	b = b % p
	return b

def MultInv(a):
	# b = a^(-1) (mod p)
	global p
	global max_iter	# test the maximum of k
	# print('EC Mult Inv')
	# print(output256bits(a))
	
	u = -p
	v = a
	r = 0
	s = 1
	k = 0
	y = 0
	x = a
	while x != 0:
		# print(x)
		if u%2 == 0:	# u_LSB == 0
			u = u >> 1
			s = s << 1
		elif v%2 == 0:	# v_LSB == 0
			v = v >> 1
			r = r << 1
		else:
			x = u + v
			y = r + s
			if (x < 0):
				u = x >> 1
				r = y
				s = s << 1
			else:
				v = x >> 1
				s = y
				r = r << 1
		k = k + 1
		# if (u >= 0):
		# 	print('u > 0')
		# if v < 0:
		# 	print('v < 0')
		# if r < 0:
		# 	print('r < 0')
		# elif r > p:
		# 	print('r > p')
		# if s < 0:
		# 	print('s < 0')
		# elif s > p:
		# 	print('s > p')

	if k > max_iter:
		max_iter = k
		
	# calculate (2^(-1) mod p)^k mod p
	# That is, in python: adjust = (p//2 + 1)**k mod p
	adjust = 1
	factor = (p//2 + 1)
	bitlen_k = 9							# empirical value
	for i in reversed(range(bitlen_k)):
		adjust = mod256(adjust * adjust, p)	# adjust^2 mod p
		if (k >> (i)) % 2:					# k[i] == 1?
			adjust = mod256(factor * adjust, p)
			# print(f'2^(-{i}): ' + output256bits(adjust))

	# print(f'adjust \t{output256bits(adjust)}')
	# print('2^(-k): ' + output256bits(adjust))
	# r may larger than p, when a is really large
	# r = mod256(r, p)
	if (r//p > 1):
		print(f'In MultInv r > 2p')
	if (r > p): 
		r = r - p
	# print(f'r\t{output256bits(r)}')
		
	b = mod256(r * adjust, p)
	# print('r * adj: ' + output256bits(b))

	b = p - b
	# print('Z^(-1): ' + output256bits(b))

	return b

#### output formatting function ####
def output256bits(Bin):
	out = f'{(Bin >> (7 * 32)) % (2**32):08X}'
	# out = ''
	for i in range(6, -1, -1):
		out = out + f'_{(Bin >> (i*32)) % (2**32):08X}'
	out = out + '\n'
	return out

def outputDBL_Regfile(A, B, C, D, E, F, G, H, I):
	if True:
		f_pointdbl.write(output256bits(A))
		f_pointdbl.write(output256bits(B))
		f_pointdbl.write(output256bits(C))
		f_pointdbl.write(output256bits(D))
		f_pointdbl.write(output256bits(E))
		f_pointdbl.write(output256bits(F))
		f_pointdbl.write(output256bits(G))
		f_pointdbl.write(output256bits(H))
		f_pointdbl.write(output256bits(I))

def outputADD_Regfile(A, B, C, D, E, F, G, H, I):
	if True:
		f_pointadd.write(output256bits(A))
		f_pointadd.write(output256bits(B))
		f_pointadd.write(output256bits(C))
		f_pointadd.write(output256bits(D))
		f_pointadd.write(output256bits(E))
		f_pointadd.write(output256bits(F))
		f_pointadd.write(output256bits(G))
		f_pointadd.write(output256bits(H))
		f_pointadd.write(output256bits(I))

####### Verify Multiplicative Inverse ########
# p = 97
max_iter = 0
if False:
	errornum = 0
	rdnum = np.random.randint((2**32), dtype=np.uint64)
	for i in range(7):
		rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))

	base = rdnum
	for a in range(base, base + 1):
		b = MultInv(a)
		# using Fermat's little law to verify
		# b_true = mod(a**(p-2), p)
		b_true = 1
		for exp in reversed(range(L)):
			b_true = b_true**2 % p
			if ((p-2) >> (exp)) % 2 :
				b_true = (a * b_true) % p
		 
		if b != b_true:
			errornum += 1
			print(f'\n!! Wrong !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
			print(f'{a:X} (mod p) \n\tMulInv: {output256bits(b)} \n\tFermat: {output256bits(b_true)}')
			break
		print(f'{a:X} (mod p) \n\tMulInv: {output256bits(b)} \n\tFermat: {output256bits(b_true)}')
	if errornum == 0:
		print('Nice! There is no error!')
	print('Max MultInv Iteration = ', max_iter)

###### Verify Modulo ######
if False:
	for i in range(1000):
		rdnum = np.random.randint((2**32), dtype=np.uint64)
		for i in range(7):
			rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
		n1 = rdnum % p

		rdnum = np.random.randint((2**32), dtype=np.uint64)
		for i in range(7):
			rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
		n2 = rdnum % p

		# print(output256bits(n1))
		# print(output256bits(n2))
		remainder = mod256(n1 * n2, p)
		true_remain = (n1 * n2) % p
		print(remainder == true_remain)
		if remainder != true_remain:
			print(output256bits(remainder))
			print(output256bits(true_remain))
			break

##### Key Exchange Simulation ---------------------------------
# Alice to Bob
if False:
	privateKey1 = 2
	privateKey2 = p//2 + 1
	dA = privateKey1
	Q_A = ECpoint_Scale(G, dA)
	print(f"dA: {dA} x P: \n\t[{Q_A[0]:48X},\n\t {Q_A[1]:48X}]")

	dB = privateKey2
	Q_B = ECpoint_Scale(Q_A, dB)
	print(f"dB: {dB} x QA: \n\t[{Q_B[0]:48X},\n\t {Q_B[1]:48X}]")

	key_A = Q_B

	# Bob to Alice
	dA = privateKey2
	Q_A = ECpoint_Scale(G, dA)
	print(f"dA: {dA} x P: \n\t[{Q_A[0]:48X},\n\t {Q_A[1]:48X}]")

	dB = privateKey1
	Q_B = ECpoint_Scale(Q_A, dB)
	print(f"dB: {dB} x QA: \n\t[{Q_B[0]:48X},\n\t {Q_B[1]:48X}]")

	key_B = Q_B

	## Verify
	if key_A == key_B:
		print('\n---------------------------------')
		print(f' Alice and Bob get the same key!')
		print('---------------------------------\n')
##### Key Exchange Simulation ---------------------------------

##### Verify HDL-like code -----------------------------------
start = time.time()
rounds = 1
if False:
	for i in range(rounds):
		print(f'\nround: {i:=5} --------------------------------------------')
		rdnum = np.random.randint((2**32), dtype=np.uint64)
		for i in range(7):
			rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
		key1 = rdnum
		print(f'key1 = {output256bits(key1)}')
		
		PubKey_pseu = ECpoint_Scale(G, key1)
		PubKey_veri = ECpoint_Scale_HDL(G, key1)
		
		if PubKey_pseu[:] == PubKey_veri[0:2+1]:
			print(f'identical PubKey? {PubKey_pseu[:] == PubKey_veri[0:2+1]}')
		else:
			print('PubKey inequal')
			break
		
		# key2 = np.random.randint((2**31)) 
		rdnum = np.random.randint((2**32), dtype=np.uint64)
		for i in range(7):
			rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
		key2 = rdnum
		print(f'key2 = {output256bits(key2)}')
		
		Ans_pseu = ECpoint_Scale(PubKey_pseu, key2)
		Ans_veri = ECpoint_Scale_HDL(PubKey_veri, key2)
		
		# print(f"True: Q:[{Ans_pseu[0]:48X},\n         {Ans_pseu[1]:48X}]")
		# print(f"Veri: Q:[{Ans_veri[0]:48X},\n         {Ans_veri[1]:48X}]")
		
		if Ans_pseu[:] == Ans_veri[0:2+1]:
			print(f'identical answer? {Ans_pseu[:] == Ans_veri[0:2+1]}')
		else:
			print(f'Answer mismatch')
			break

		# Check the cyclic property
		key_equiv = key1*key2 % n
		print(f'Equiv key = {key_equiv:X}')
		Equiv_pseu = ECpoint_Scale(G, key_equiv)
		Equiv_veri = ECpoint_Scale_HDL(G, key_equiv)
		# print(f"True: Q:[{Equiv_pseu[0]:48X},\n         {Equiv_pseu[1]:48X}]")
		# print(f"Veri: Q:[{Equiv_veri[0]:48X},\n         {Equiv_veri[1]:48X}]")
		print(f'equal to answer? {Equiv_veri[:] == Ans_pseu[0:2+1]}')
		if Equiv_veri[:] != Ans_pseu[0:2+1]:
			print('equivalent Key != (key1 * key2) mod n')
			print('WRONG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
			break
		print(f'---------------------------------------------------------')
	end = time.time()
	print(f'\nAvg elapsed time: {(end-start)/rounds * 1000:3.3f} ms/round')

### Verify Doubling
if False:
	P = [0x1CE62661480AF4D8F518AF5EEDBA7C26A4284361AF6CB45C, \
		 0x728798031F02328DE92BF9461E0D06BFB30A98DC09F2C492]
	Q_true = ECpoint_Doubling([P[0], P[1], 1])
	[A, B, C, D] = [P[0], P[1], 1, mod(P[0] * P[0], p)]
	Q_hdl = ECpoint_Doubling_HDL(A, B, C, D, 0, 0, 0, 0, 0)
	print(f"True: Q:[{Q_true[0]:48X},\n         {Q_true[1]:48X}\n         {Q_true[2]:48X}]")
	print(f"Veri: Q:[{Q_hdl[0]:48X},\n         {Q_hdl[1]:48X}\n         {Q_hdl[2]:48X}]")
	print(Q_true[:] == Q_hdl[0:2+1])
	Z_inv = MultInv(Q_true[2])	# cal: 1/Z mod p
	Q = [mod(Q_true[0] * Z_inv**2, p), mod(Q_true[1] * Z_inv**3, p)];
	print(f'Q: [{Q[0]:48X}\n    {Q[1]:48X}]')
	
### Verify Addition
if False:
	G2 = [5518421071991463839100833553767507511692256595899225500678, \
		   623358079454615111697064773677714472814391500098424913986, \
		   348100664587244062809715104560438820728046977854773301282, \
		  5432392459578811477078447929571628010766186319035089069573, \
		  3861615048567063030691355956960560842237002122943048358239, \
		  1432760983132731149600004535962129452791679972009845012997, \
		   602046282375688656758213480587526111916698976636884684818, \
		   174050332293622031404857552280219410364023488927386650641, \
																   1]
	Q_true = ECpoint_Addition([G2[0], G2[1], G2[2]], [G2[6], G2[7], G2[8]])
	Q_hdl = ECpoint_Addition_HDL(G2[0], G2[1], G2[2], G2[3], G2[4], G2[5], G2[6], G2[7], G2[8])
	print(f"True: Q:[{Q_true[0]:48X},\n         {Q_true[1]:48X}\n         {Q_true[2]:48X}]")
	print(f"Veri: Q:[{Q_hdl[0]:48X},\n         {Q_hdl[1]:48X}\n         {Q_hdl[2]:48X}]")
	print(Q_true[:] == Q_hdl[0:2+1])

### ECDH batch testbench file generator ---------------------------------------
testnum = 0

for rnd in range(testnum):
	rdnum = np.random.randint((2**32), dtype=np.uint64)
	for i in range(7):
		rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
	prvkeyA = rdnum
	# print(output256bits(prvkeyA))
	pubkeyA = ECpoint_Scale_HDL(G, prvkeyA)
	pubkeyA_golden = ECpoint_Scale(G, prvkeyA)

	# print(f"True: Q:[{pubkeyA_golden[0]:48X},\n         {pubkeyA_golden[1]:48X}]")
	# print(f"Veri: Q:[{pubkeyA[0]:48X},\n         {pubkeyA[1]:48X}]")
	if pubkeyA != pubkeyA_golden:
		print("WRONG pubkeyA !!!!!!!!!!!!!!")

	for j in range(2):
		f_base.write(output256bits(G[j]))
		f_end.write(output256bits(pubkeyA[j]))
	f_scalar.write(output256bits(prvkeyA))
	
	rdnum = np.random.randint((2**32), dtype=np.uint64)
	for i in range(5):
		rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
	prvkeyB = rdnum
	# prvkeyB = 7
	pubkeyB = ECpoint_Scale_HDL(pubkeyA, prvkeyB)
	pubkeyB_golden = ECpoint_Scale(pubkeyA, prvkeyB)
	if pubkeyB != pubkeyB_golden:
		print("WRONG pubkeyB !!!!!!!!!!!!!!")
	for j in range(2):
		f_base.write(output256bits(pubkeyA[j]))
		f_end.write(output256bits(pubkeyB[j]))
	f_scalar.write(output256bits(prvkeyB))


## Output MultInv data --------------------------------------------------------
def MultInv_test(a, file):
	# b = a^(-1) (mod p)
	global p
	u = -p
	v = a
	r = 0
	s = 1
	k = 0
	y = 0
	x = a
	while x != 0:
		# print(x)
		if u%2 == 0:	# u_LSB == 0
			u = u >> 1
			s = s << 1
		elif v%2 == 0:	# v_LSB == 0
			v = v >> 1
			r = r << 1
		else:
			x = u + v
			y = r + s
			if (x < 0):
				u = x >> 1
				r = y
				s = s << 1
			else:
				v = x >> 1
				s = y
				r = r << 1
		k = k + 1
		# print(k)
		mask = 0x1_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF
		file.write(f'{u & mask:065X}\n')
		file.write(f'{v & mask:065X}\n')
		file.write(f'{s & mask:065X}\n')
		file.write(f'{r & mask:065X}\n')
		file.write(f'{x & mask:065X}\n')
		file.write(f'{y & mask:065X}\n')
	print(output256bits(r))
# innum = 0x62C151D9_4EE0B260_83265B9F_49837211_A840A5E2_3D96598C_89A66817_5A3AB570
# f_multinv = open('./pat256/multinv.dat', 'w')
# MultInv_test(innum, f_multinv)
# # print(output256bits(out))
# f_multinv.close()

## 
# key = 0xB84FDB87_06381800_3AA8780D_EB517EB7_C258F039_0098BE8E_13110263_B4B575AC
# out = ECpoint_Scale_HDL(G, key)
# print(output256bits(out[0]))
# print(output256bits(out[1]))

### Top Module Testbench ------------------------------------------------------------
def output256to128bits(Bin):
	out = f'{(Bin >> (7 * 32)) % (2**32):08X}'
	for i in range(6, 3, -1):
		out = out + f'_{(Bin >> (i*32)) % (2**32):08X}'
	out = out + '\n'
    
	out = out + f'{(Bin >> (3 * 32)) % (2**32):08X}'
	for i in range(2, -1, -1):
		out = out + f'_{(Bin >> (i*32)) % (2**32):08X}'
	out = out + '\n'
	return out

f_top_inpoint = open('./TP/ECDH_test_point.dat', 'w')
f_top_outpoint = open('./TP/ECDH_out_point.dat', 'w')
f_top_key = open('./TP/ECDH_test_key.dat', 'w')

# rdnum = np.random.randint((2**32), dtype=np.uint64)
# for i in range(7):
#     rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
# prvkeyA = rdnum
prvkeyA = 0x782D1AAD_45ED0C00_CDBAC6CE_031F30BD_32369913_4251E6F6_8739170B_D9CB931F
pubkeyA = ECpoint_Scale_HDL(G, prvkeyA)
pubkeyA_golden = ECpoint_Scale(G, prvkeyA)

print(f"True: Q:[{pubkeyA_golden[0]:48X},\n         {pubkeyA_golden[1]:48X}]")
print(f"Veri: Q:[{pubkeyA[0]:48X},\n         {pubkeyA[1]:48X}]")
if pubkeyA != pubkeyA_golden:
    print("WRONG pubkeyA !!!!!!!!!!!!!!")

for j in range(2):
    f_top_inpoint.write(output256to128bits(G[j]))
    f_top_outpoint.write(output256to128bits(pubkeyA[j]))
    f_base.write(output256bits(G[j]))
    f_end.write(output256bits(pubkeyA[j]))
f_top_key.write(output256to128bits(prvkeyA))
f_scalar.write(output256bits(prvkeyA))

# rdnum = np.random.randint((2**32), dtype=np.uint64)
# for i in range(7):
#     rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
# prvkeyB = rdnum
prvkeyB = 0xBEB55F6B_F3249800_D17B6718_E94CAD3F_7A837038_D2E501D2_1311FFE1_B9174865
pubkeyB = ECpoint_Scale_HDL(pubkeyA, prvkeyB)
pubkeyB_golden = ECpoint_Scale(pubkeyA, prvkeyB)

if pubkeyB != pubkeyB_golden:
    print("WRONG pubkeyB !!!!!!!!!!!!!!")
for j in range(2):
    f_top_inpoint.write(output256to128bits(pubkeyA[j]))
    f_top_outpoint.write(output256to128bits(pubkeyB[j]))
    f_base.write(output256bits(pubkeyA[j]))
    f_end.write(output256bits(pubkeyB[j]))
f_top_key.write(output256to128bits(prvkeyB))
f_scalar.write(output256bits(prvkeyB))


### File close
if outputfile:
	f_base.close()
	f_end.close()
	f_scalar.close()
	data_file.close()
	f_pointadd.close()
	