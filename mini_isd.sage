# A simple toy ISD
# PK - 2018

#TODO Leon
#TODO Stern w/ limited memory (liek < 200 MB? Does that make sense??)

def find_lo_weight(C,maxiter=1000000,wtarget=1,wcomb=2,nthreads=1):
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

	@parallel(nthreads)
	def do_Iset(curminw):
		set_random_seed()
		while True:
			Iset = sample(supp,k)
			Gis = G.matrix_from_columns(Iset)
			if (Gis.is_invertible()):
				break
		Gis_inv = Gis.inverse()
		Glw = Gis_inv * G
		curminwd = None
		# weight one is simple
		for c in range(k):
			cc = Glw.row(c)
			cw = cc.hamming_weight()
			if cw < curminw:
				print "Found a new codeword of weight "+str(cw)
				curminwd = cc
				curminw = cw
		# further weights wiz Gray codes
		for w in range(wcomb)[1:]:
			base = [1]*w
			base.extend([0]*(k-w))
			v = vector(base)
			cc = v*Glw
			cw = cc.hamming_weight() #ugly repeat 1
			if cw < curminw:
				print "Found a new codeword of weight "+str(cw)
				curminwd = cc
				curminw = cw
			for i,j in combinations(k,w):
				cc = cc + Glw.row(i)
				cc = cc + Glw.row(j)
				cw = cc.hamming_weight()
				if cw < curminw: #ugly repeat 2
					print "Found a new codeword of weight "+str(cw)
					curminwd = cc
					curminw = cw
		return (curminw,curminwd)

	while it < maxiter:
		if 1 == nthreads:
			(cw,cc) = do_Iset(minw)
		else:
			run = do_Iset([minw]*nthreads)
			res = list(run)
			(cw,cc) = min(res)[1]
		if (cw < minw):
			minw = cw
			minwd = cc
		# done once per iteration, as not expected to be successful many times?
		if minw <= wtarget:
			return minwd
		it += nthreads
	return minwd
