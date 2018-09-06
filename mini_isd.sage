# A simple toy ISD
# PK - 2018

#TODO decoding wrapper (dimension extend? error sampler?)

def find_lo_weight(C,maxiter=1000000,wtarget=1,wcomb=2,nthreads=1):
	from sage.combinat.gray_codes import combinations

	if GF(2) != C.base_field():
		raise NotImplementedError

	n,k = C.length(), C.dimension()
	G = C.generator_matrix()
	supp = range(n)
	minwd = G.row(0) # arbitrary
	minw = minwd.hamming_weight()
	it = 0
	nthreads = min(maxiter, nthreads)

	@parallel(nthreads)
	def do_Iset(curminw):
		set_random_seed()
		while True: # TODO kind of optim?
			Iset = sample(supp,k)
			Gis = G.matrix_from_columns(Iset)
			if Gis.is_invertible():
				break
		Gis_inv = Gis.inverse()
		Glw = Gis_inv * G
		curminwd = None
		# weight one is simple
		for c in range(k):
			cc = Glw.row(c)
			cw = cc.hamming_weight()
			if cw < curminw:
				curminwd,curminw = cc,cw
		# further weights wiz Gray codes
		for w in range(wcomb)[1:]:
			base = [1]*w
			base.extend([0]*(k-w))
			v = vector(base)
			cc = v*Glw
			cw = cc.hamming_weight() #ugly repeat 1
			if cw < curminw:
				curminwd,curminw = cc,cw
			for i,j in combinations(k,w):
				cc += Glw.row(i)
				cc += Glw.row(j)
				cw = cc.hamming_weight()
				if cw < curminw: #ugly repeat 2
					curminwd,curminw = cc,cw
		return (curminw,curminwd)

	while it < maxiter:
		if 1 == nthreads:
			(cw,cc) = do_Iset(minw)
		else:
			run = do_Iset([minw]*nthreads)
			res = list(run)
			(cw,cc) = min(res)[1]
		if cw < minw:
			minwd,minw = cc,cw
			print "Found a new codeword of weight "+str(cw)
		# done once per iteration, as not expected to be successful many times?
		if minw <= wtarget:
			return minwd
		it += nthreads
	return minwd



# bof
def find_lo_weight_birthday(C,matchsize,maxiter=1000000,wtarget=1,wcomb=2,nthreads=1):
	from sage.combinat.gray_codes import combinations

	if GF(2) != C.base_field():
		raise NotImplementedError

	n,k = C.length(), C.dimension()
	G = C.generator_matrix()
	supp = range(n)
	minwd = G.row(0) # arbitrary
	minw = minwd.hamming_weight()
	it = 0

	def update_app_HT(ht, key, val):
		try:
			a = ht[key]
			a.append(val)
			ht[key] = a
		except KeyError:
			ht[key] = [val]

	@parallel(nthreads)
	def do_Iset(curminw):
		set_random_seed()
		while True:
			Iset = sample(supp,k+matchsize)
			Gis = G.matrix_from_columns(Iset[0:k])
			if (Gis.is_invertible()):
				break
		Gis_inv = Gis.inverse()
		Glw = Gis_inv * G
		curminwd = None
		t_upper = {}
		t_lower = {}
		Matchset = Iset[k:k+matchsize]
		mid = int(k/2)
		# weight one is simple
		for c in range(k):
			cc = Glw.row(c)
			cc_match = vector(GF(2),cc.list_from_positions(Matchset))
			cc_match.set_immutable()
			if c < mid:
				update_app_HT(t_upper,cc_match,cc)
			else:
				update_app_HT(t_lower,cc_match,cc)
		# further weights wiz Gray codes
		# all bundled together, regardless of the weight
		# no randomization of the split for now
		for w in range(wcomb)[1:]:
			base_upper = [1]*w
			base_upper.extend([0]*(k-w))
			base_lower = [0]*mid
			base_lower.extend([1]*w)
			base_lower.extend([0]*(k-mid-w))
			v_upper = vector(base_upper)
			v_lower = vector(base_lower)
			cc_upper = v_upper*Glw
			cc_lower = v_lower*Glw
			cc_match = vector(GF(2),cc_upper.list_from_positions(Matchset))
			cc_match.set_immutable()
			update_app_HT(t_upper,cc_match,cc)
			cc_match = vector(GF(2),cc_lower.list_from_positions(Matchset))
			cc_match.set_immutable()
			update_app_HT(t_lower,cc_match,cc)
			for i,j in combinations(mid,w):
				cc_upper += Glw.row(i)
				cc_upper += Glw.row(j)
				cc_match = vector(GF(2),cc_upper.list_from_positions(Matchset))
				cc_match.set_immutable()
				update_app_HT(t_upper,cc_match,cc)
			for i,j in combinations(k-mid,w):
				cc_lower += Glw.row(mid+i)
				cc_lower += Glw.row(mid+j)
				cc_match = vector(GF(2),cc_lower.list_from_positions(Matchset))
				cc_match.set_immutable()
				update_app_HT(t_lower,cc_match,cc)
		# now teg matchingz
		for m in iter(t_upper):
			try:
				cc_low = t_lower[m]
				cc_upp = t_upper[m]
				for ccu in cc_upp:
					for ccl in cc_low:
						cc = ccu + ccl
						cw = cc.hamming_weight()
						if (cw > 0) and (cw < curminw):
							curminwd = cc
							curminw = cw
			except KeyError:
				continue
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
			print "Found a new codeword of weight "+str(cw)
		# done once per iteration, as not expected to be successful many times?
		if minw <= wtarget:
			return minwd
		it += nthreads
	return minwd
