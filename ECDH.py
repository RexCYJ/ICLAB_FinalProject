# ECDH implementation
import numpy as np
import time

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
	p = int('FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFE_FFFFFFFF_FFFFFFFF', base=16)
	n = int('FFFFFFFF_FFFFFFFF_FFFFFFFF_99DEF836_146BC9B1_B4D22831', base=16)
	h = 1
	a = int('FFFFFFFF_FFFFFFFF_FFFFFFFF_FFFFFFFE_FFFFFFFF_FFFFFFFC', base=16)
	b = int('64210519_E59C80E7_0FA7E9AB_72243049_FEB8DEEC_C146B9B1', base=16)
	G = [int('188DA80E_B03090F6_7CBF20EB_43A18800_F4FF0AFD_82FF1012', base=16), \
		 int('07192B95_FFC8DA78_631011ED_6B24CDD5_73F977A1_1E794811', base=16)]
	L = 192
##### ------------------------------------------------------------------------------

max_iter = 0

data_file = open("./pat/reg_record.dat", 'w')
f_pointadd = open("./pat/pointadd.dat", 'w')

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
	# print(f'M = {M}')
	# print(f'S = {S}')
	# print(f'T = {T}')
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

def ECpoint_Scale_HDL(P, k):
	# Q = [k]P
	# to jacobian coordinate
	global p
	global L
	[X0, Y0, Z0] = [0, 0, 0];			# current point	P
	[X1, Y1, Z1] = [P[0], P[1], 1];		# base point	Q

    # only need to initialize G, H, I
	# [A, B, C, D, E, F, G, H, I] = [0, 0, 0, 0, 0, 0, X1, Y1, 1] # [X0, Y0, Z0, X0^2, 0, 0, X1, Y1, 1]
	[A, B, C, D, E, F, G, H, I] = [0, 0, 0, 0, 0, 0, 0, 0, 0] # [X0, Y0, Z0, X0^2, 0, 0, X1, Y1, 1]
	is_first = 1

	kbitarray = list(np.binary_repr(k, width=L))
	kbitarray.reverse()
	for j in reversed(range(L)):		# for L-1, L-2, ..., 1, 0
		if (not is_first):
			[A, B, C, D, E, F, G, H, I] = \
				ECpoint_Doubling_HDL(A, B, C, D, E, F, G, H, I)	# double Q: Q := 2Q
			# data_file.write(output192bits(A))
			# data_file.write(output192bits(B))
			# data_file.write(output192bits(C))
			# data_file.write(output192bits(D))
		if int(kbitarray[j]) == 1:
			if (is_first == 1):		# first visit
				is_first = 0
				A, B, C = X1, Y1, 1
				D = mod192(A * A, p)	# prepare X0^2
				# print(output192bits(D))
			else:	
				# add P: Q := Q + P
				G, H, I = X1, Y1, 1		# refresh
				[A, B, C, D, E, F, G, H, I] = \
					ECpoint_Addition_HDL(A, B, C, D, E, F, G, H, I)
				# data_file.write(output192bits(A))
				# data_file.write(output192bits(B))
				# data_file.write(output192bits(C))
				# data_file.write(output192bits(D))
		else:
			if (not is_first):
				# print(f'{A * A:96X}')
				D = mod192(A * A, p)
				# print(output192bits(D))
	
	# data_file.write(output192bits(A))
	# data_file.write(output192bits(B))
	# data_file.write(output192bits(C))

	# transform back to normal domain
	Z_inv = MultInv(C)	# cal: 1/Z mod p
    # [x y] = [X / Z^2,  Y / Z^3]
	Z_inv2 = mod192(Z_inv * Z_inv, p)
	Z_inv3 = mod192(Z_inv2 * Z_inv, p)
	# print('Z^(-2): ' + output192bits(Z_inv2))
	# print('Z^(-3): ' + output192bits(Z_inv3))
	Q = [mod192(A * Z_inv2, p), mod192(B * Z_inv3, p)];
	# print('Qx: ' + output192bits(Q[0]))
	# print('Qy: ' + output192bits(Q[1]))

	return Q

def ECpoint_Doubling_HDL(A, B, C, D, E, F, G, H, I):
	# input: A = X1, B = Y1, C = Z1, D = X1^2
	global p
	global a
	# print([A, B, C])
	E = mod192(3 * D, p)	# 3 * X1^2
	# t1 ---------------------------
	D = mod192(C * C, p)	# Z1^2
	# t2 ---------------------------
	D = mod192(D * D, p)	# Z1^4
	# t3 ---------------------------
	D = mod192(a * D, p)	# a * Z1^4
	# t4 ---------------------------
	E = mod192(D + E, p)	# M = 3X1^2+aZ1^4
	D = mod192(B * C, p)	# Y1*Z1
	# t5 ---------------------------
	B = mod192(B * B, p)	# Y1^2
	C = mod192(D + D, p)	# Z2 $$
	# t6 ---------------------------
	D = mod192(B * A, p)	# X1*Y1^2
	# t7 ---------------------------
	D = mod192(4 * D, p)	# S = 4X1Y1^2
	# t8 ---------------------------
	A = mod192(D + D, p)	# 2S
	F = mod192(E * E, p)	# M^2
	# t9 ---------------------------
	B = mod192(B * B, p)	# Y1^4
	A = mod192(F + p - A, p)	# X2 $$
	# t10 --------------------------
	B = mod192(8 * B, p)	# T
	D = mod192(D + p - A, p)	# S - X2
	# t11 --------------------------
	D = mod192(D * E, p)	# M(S-X2)
	# t12 --------------------------
	B = mod192(D + p - B, p)	# Y2
	D = mod192(C * C, p)	# Z2^2
	# t13 --------------------------
	# output: A = X2, B = Y2, C = Z2, D = Z2^2
	# did not use G, H, I
	return [A, B, C, D, E, F, G, H, I]

def ECpoint_Addition_HDL(A, B, C, D, E, F, G, H, I):
	# input: A = X0, B = Y0, C = Z0, D = Z0^2, 
	# 		 G = X1, H = Y1, I = Z1
	# P0 = Q_cur, P1 = base point
	global p
	E = mod192(G * D, p)	# U1 = X1*Z0^2
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t01 ---------------------------
	D = mod192(D * C, p)	# Z0^3
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t02 ---------------------------
	I = mod192(A + p - E, p)	# W = U0 - U1
	D = mod192(D * H, p)	# S1 = Y1 * Z0^3
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t03 ---------------------------
	F = mod192(B + p - D, p)	# R = S0 - S1
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t04 ---------------------------
	G = mod192(F * F, p)	# R^2
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t05 ---------------------------
	A = mod192(A + E, p)	# T = U0 + U1
	E = mod192(I * I, p)	# W^2
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t06 ---------------------------
	D = mod192(D + B, p)	# M
	B = mod192(A * E, p)	# TW^2
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t07 ---------------------------
	A = mod192(G + p - B, p)	# X2 $$
	E = mod192(I * E, p)	# W^3
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t08 ---------------------------
	G = mod192(E * D, p)	# MW^3
	D = mod192(A + A, p)	# 2X2
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t09 ---------------------------
	B = mod192(B + p - D, p)	# V
	# print(f'V = {B}')
	C = mod192(C * I, p)	# Z2 $$
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t10 ---------------------------
	F = mod192(F * B, p)	# VR
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t11 ---------------------------
	D = mod192(A * A, p)	# X2^2
	F = mod192(F + p - G, p)	# 2Y2
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t12 ---------------------------
	B = mod192(F * (p//2 + 1), p)	# Y2 $$
	outputRegfile(A, B, C, D, E, F, G, H, I)
	# t13 ---------------------------
	# output: A = X2, B = Y2, C = Z2, D = X2^2
	# remenber to fresh G and H with P1
	return [A, B, C, D, E, F, G, H, I]

def mod(a, p):
	# b = a mod p
	b = a % p
	return b

def mod192(a, p):
	# b = a mod p
	# b = a % p
	# return b
	if a >= (p*p):
		print('\nout of bound\n')
	elif a < 0:
		print('\na < 0\n')
	A = [0, 0, 0, 0, 0, 0]
	# A[0] = a % (2**64)	# A[0] = a[0:63]
	for i in range(6):
		A[i] = (a >> (64*i)) % (2**64)	# A[i] = a[i*64 +: 64]
	T  = A[2] * (2**(64*2)) + A[1] * (2**(64*1)) + A[0]	# T  = {A2, A1, A0};
	S1 = 					  A[3] * (2**(64*1)) + A[3]	# S1 = 	   {A3, A3};
	S2 = A[4] * (2**(64*2)) + A[4] * (2**(64*1))		# S2 = {A4, A4,  0};
	S3 = A[5] * (2**(64*2)) + A[5] * (2**(64*1)) + A[5] # S3 = {A5, A5, A5};

	b = (T + S1 + S2 + S3)
	# print('b: ' + output192bits(b))
	# print(f'{b:X}')
	quotient = (T + S1 + S2 + S3) // p		# quotient is usually 0, 1, 2; 
											# theoritically, might be up to 3
											# but I dont want to consider
	# print(quotient)
	if (b < p):
		return b
	elif (b - p) < p:
		return b - p
	else:
		if quotient > 2:
			print('quotient > 2')
		return b - 2*p

	b = (T + S1 + S2 + S3) % p
	return b

def MultInv(a):
	# b = a^(-1) (mod p)
	global p
	global max_iter	# test the maximum of k
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
		# print(f'u\t{u:X}')
		# print(f'v\t{v:X}')
		# print(f's\t{s:X}')
		# print(f'r\t{r:X}')
		# print(f'x\t{x:X}')
		# print(f'y\t{y:X}')

	if k > max_iter:
		max_iter = k
	
	# print(k)
	# print(f'u\t{u:X}')
	# print(f'v\t{v:X}')
	# print(f's\t{s:X}')
	# print(f'r\t{r:X}')
	# print(f'x\t{x:X}')
	# print(f'y\t{y:X}')

	# calculate (2^(-1) mod p)^k mod p
	# That is, in python: adjust = (p//2 + 1)**k mod p
	adjust = 1
	factor = (p//2 + 1)
	bitlen_k = 9							# empirical value
	for i in reversed(range(bitlen_k)):
		adjust = mod192(adjust * adjust, p)	# adjust^2 mod p
		if (k >> (i)) % 2:					# k[i] == 1?
			adjust = mod192(factor * adjust, p)
			# print(f'2^(-{i}): ' + output192bits(adjust))

	# print('2^(-k): ' + output192bits(adjust))
	# r may larger than p, when a is really large
	# r = mod192(r, p)
	if (r > p): 
		r = r - p
	b = mod192(r * adjust, p)
	# print('r * adj: ' + output192bits(b))

	b = p - b
	# print('Z^(-1): ' + output192bits(b))

	return b

#### output formatting function ####
def output192bits(Bin):
	out = f'{(Bin >> (11* 16)) % (2**16):04X}'
	for i in range(10, -1, -1):
		out = out + f'_{(Bin >> (i*16)) % (2**16):04X}'
	out = out + '\n'
	return out

def outputRegfile(A, B, C, D, E, F, G, H, I):
	if False:
		f_pointadd.write(output192bits(A))
		f_pointadd.write(output192bits(B))
		f_pointadd.write(output192bits(C))
		f_pointadd.write(output192bits(D))
		f_pointadd.write(output192bits(E))
		f_pointadd.write(output192bits(F))
		f_pointadd.write(output192bits(G))
		f_pointadd.write(output192bits(H))
		f_pointadd.write(output192bits(I))

####### Verify Multiplicative Inverse ########
# p = 97
max_iter = 0
if False:
	errornum = 0
	base = 0xFFFF_FFFF_1240_3121_1105_1059_5155_6666_9892_1241_4623
	for a in range(base, base + 10000):
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
			print(f'{a:X} (mod p) \n\tMulInv: {b:48X} \n\tFermat: {b_true:48X}')
			break
		print(f'{a:X} (mod p) \n\tMulInv: {b:48X} \n\tFermat: {b_true:48X}')
	if errornum == 0:
		print('Nice! There is no error!')
	print('Max MultInv Iteration = ', max_iter)

###### Verify Modulo ######
if False:
	for i in range(100):
		n1 = (np.random.randint(2**31) * (2**33) + np.random.randint(2**31)) % p
		n2 = (np.random.randint(2**31) * (2**33) + np.random.randint(2**31)) % p
		# n2 = 0
		remainder = mod192(n1 * n2, p)
		true_remain = (n1 * n2) % p
		print(remainder)
		print(true_remain)
		print(remainder == true_remain)
		if remainder != true_remain:
			print(False)
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
rounds = 5
if False:
    for i in range(rounds):
        print(f'\nround: {i:=5} --------------------------------------------')
        rdnum2 = np.random.randint((2**63), dtype=np.uint64)
        rdnum1 = np.random.randint((2**63), dtype=np.uint64)
        rdnum0 = np.random.randint((2**63), dtype=np.uint64) 
        key1 = rdnum2 * (2**128) + rdnum1 * (2**64) + rdnum0
        print(f'key1 = {rdnum2:016X}_{rdnum1:016X}_{rdnum0:016X}')
        
        PubKey_pseu = ECpoint_Scale(G, key1)
        PubKey_veri = ECpoint_Scale_HDL(G, key1)
        
        if PubKey_pseu[:] == PubKey_veri[0:2+1]:
            print(f'identical PubKey? {PubKey_pseu[:] == PubKey_veri[0:2+1]}')
        else:
            print('PubKey inequal')
            break
        
        # key2 = np.random.randint((2**31)) 
        rdnum2 = np.random.randint((2**63), dtype=np.uint64)
        rdnum1 = np.random.randint((2**63), dtype=np.uint64)
        rdnum0 = np.random.randint((2**63), dtype=np.uint64) 
        key2 = rdnum2 * (2**128) + rdnum1 * (2**64) + rdnum0
        print(f'key2 = {rdnum2:016X}_{rdnum1:016X}_{rdnum0:016X}')
        
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

### Testbench file generator -------------------------------------------------------------------
# rdnum2 = np.random.randint((2**63), dtype=np.uint64)
# rdnum1 = np.random.randint((2**63), dtype=np.uint64)
# rdnum0 = np.random.randint((2**63), dtype=np.uint64) 
# print(f'key1 = {rdnum2:016X}_{rdnum1:016X}_{rdnum0:016X}')
f_base = open('./pat/basepoint.dat', 'w')
f_end = open('./pat/endpoint.dat', 'w')
f_scalar = open('./pat/scalar.dat', 'w')

testnum = 500

for rnd in range(testnum):
	rdnum = np.random.randint((2**32), dtype=np.uint64)
	for i in range(5):
		rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))

	prvkeyA = rdnum
	# prvkeyA = 0xB7D4_F2DE_8EC4_C000_9B51_200A_D3D8_FF2A_A32E_A163_BA78_7018
	# print(output192bits(prvkeyA))
	pubkeyA = ECpoint_Scale_HDL(G, prvkeyA)
	# pubkeyA_golden = ECpoint_Scale(G, prvkeyA)

	# print(f"True: Q:[{pubkeyA_golden[0]:48X},\n         {pubkeyA_golden[1]:48X}]")
	# print(f"Veri: Q:[{pubkeyA[0]:48X},\n         {pubkeyA[1]:48X}]")

	for j in range(2):
		f_base.write(output192bits(G[j]))
		f_end.write(output192bits(pubkeyA[j]))
	f_scalar.write(output192bits(prvkeyA))
	
	rdnum = np.random.randint((2**32), dtype=np.uint64)
	for i in range(5):
		rdnum = int(rdnum*(2**32) + np.random.randint((2**32), dtype=np.uint64))
	prvkeyB = rdnum
	pubkeyB = ECpoint_Scale_HDL(pubkeyA, prvkeyB)
	for j in range(2):
		f_base.write(output192bits(pubkeyA[j]))
		f_end.write(output192bits(pubkeyB[j]))
	f_scalar.write(output192bits(prvkeyB))

f_base.close()
f_end.close()
f_scalar.close()

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
		mask = 0x1_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF_FFFF
		file.write(f'{u & mask:049X}\n')
		file.write(f'{v & mask:049X}\n')
		file.write(f'{s & mask:049X}\n')
		file.write(f'{r & mask:049X}\n')
		file.write(f'{x & mask:049X}\n')
		file.write(f'{y & mask:049X}\n')
	print(output192bits(r))
# innum = 0xFFFF
# f_multinv = open('./pat/multinv.dat', 'w')
# MultInv_test(innum, f_multinv)
# # print(output192bits(out))
# f_multinv.close()

data_file.close()