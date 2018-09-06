# A simple toy ISD
# PK - 2018


#TODO parallelizm
#TODO Leon
#TODO Stern w/ limited memory (liek < 200 MB? Does that make sense??)

def find_lo_weight(C,maxiter,wtarget=1,wcomb=2):
	from sage.all import sample
	from sage.combinat.gray_codes import combinations

	if GF(2) != C.base_field():
		raise NotImplementedError

	n,k = C.length(), C.dimension()
	G = C.generator_matrix()
	supp = range(n)
	minwd = G.row(0) # arbitrary
	minw = minwd.hamming_weight()
	it = 0

	while it < maxiter:
		Iset = sample(supp,k)
		Gis = G.matrix_from_columns(Iset)
		try:
			Gis_inv = Gis.inverse()
		except ZeroDivisionError:
			continue
		it += 1
		Glw = Gis_inv * G
		# weight one is simple
		for c in range(k):
			cc = Glw.row(c)
			cw = cc.hamming_weight()
			if cw < minw:
				print "Found a new codeword of weight "+str(cw)
				minwd = cc
				minw = cw
		# further weights wiz Gray codes
		for w in range(wcomb)[1:]:
			base = [1]*w
			base.extend([0]*(k-w))
			v = vector(base)
			cc = v*Glw
			cw = cc.hamming_weight() #ugly repeat 1
			if cw < minw:
				print "Found a new codeword of weight "+str(cw)
				minwd = cc
				minw = cw
			for i,j in combinations(k,w):
				cc = cc + Glw.row(i)
				cc = cc + Glw.row(j)
				cw = cc.hamming_weight()
				if cw < minw: #ugly repeat 2
					print "Found a new codeword of weight "+str(cw)
					minwd = cc
					minw = cw
		# done once per iteration, as not expected to be successful many times?
		if minw <= wtarget:
			return minwd
	return minwd
