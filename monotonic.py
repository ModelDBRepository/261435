from neuron import h, gui
from numpy import *
from scipy import optimize as opt


FI = 0 # for figure 1c
RUN = 1 # creates t, v, w file
VNULL = 1 # creates dv/dt = 0, dw/dt = 0 files
IV = 0 # for figure 1a
JAC =0 # calculates eigenvalues of fixed points
BACKWARDS = 1 # checks for bistability in FI

VF=0 # calculates vector field on fixed grid

AMP = -10 # current value - see Fig 1,2


BOUNDARY = 0 # numerically determines quasi-threshold
ATTRACTOR = 0 # numerically determines basin of attraction for bimodal
	
if JAC:
	VNULL=1 # calculation requires this to have been done

class neuron:
	def __init__(self, cType=0, vinit=None, ninit=None):
		self.soma = h.Section()
		self.soma.L     = 100.
		self.soma.diam  = 10./pi
		self.soma.insert('type21')
		#=== neuron initialization ===#
		if cType == 1 or cType == 2:
			self.soma(0.5).type21.type21 = cType
		if type(vinit) is float or type(vinit) is int:
			self.soma(0.5).v = vinit
		if type(ninit) is float or type(ninit) is int:
			self.soma(0.5).type21.ninit = ninit
		#===     set recordings     ===#
		self.v = h.Vector()
		self.v.record(self.soma(0.5)._ref_v)
		self.n = h.Vector()
		self.n.record(self.soma(0.5).type21._ref_n)
		#===    record spike times  ===#
		self.spks	= h.Vector()
		self.sptr	= h.APCount(.5, sec=self.soma)
		self.sptr.thresh = 0.#25
		self.sptr.record(self.spks)
		#==          synapse        ===#
		self.isyn  = h.Exp2Syn(0.5,sec=self.soma)
		self.isyn.tau1 = 1.
		self.isyn.tau2 = 3.
		self.isyn.e    = -75.


test = neuron(cType=1)


h.finitialize() # sets to type 1 parameters

test.soma.type21_type21 = 0 # prevent it from reseting with finitialize()
test.soma.a_type21 = 1
test.soma.b_type21 = -1
test.soma.n0_type21 = 0
test.soma.sn_type21 = 1
test.soma.el_type21 = -58.0
test.soma.gl_type21 = 0.5
test.soma.v12_type21 = -45
test.soma.sl_type21 = 10.0
test.soma.gk_type21=24
test.soma.gna_type21 = 80
test.soma.mhalf_type21 = -28.0
test.soma.mslope_type21 = 10.0 #0.5/0.055
test.soma.ninit_type21 = -1 # use ninf instead

scale = 1.0/(test.soma.ena_type21 - test.soma.ek_type21)


t0 =7.5


test.soma.v0_type21 = test.soma.v12_type21
test.soma.sg_type21 = 2*test.soma.sl_type21
test.soma.t0_type21 = t0
test.soma.st_type21 = 0
test.soma.n_type21 = 0.1
test.soma.v = -80

vstep = 0.1

varray = arange(-80,50,vstep)


c = test.soma
el = c.el_type21
gl = c.gl_type21
b = c.b_type21
a = c.a_type21
ena = c.ena_type21
ek = c.ek_type21
gna = c.gna_type21
gk = c.gk_type21
mhalf = -c.mhalf_type21
mslope = c.mslope_type21
sn = c.sn_type21
sg = c.sg_type21
sl = c.sl_type21
v12 = c.v12_type21
v0 = c.v0_type21
t0 = c.t0_type21
s = h.area(0.5,sec=test.soma)
n0 = 0
cm = c.cm
#print cm,s,el,gl, b,a,ena,ek,gna,gk,mhalf,mslope,sn,sg,sl,v12,v0,t0
h.finitialize()
ninf = lambda v: n0 + sn/(1.+exp(-(v-v12)/sl ))
ntau = lambda v: t0/(exp((v-v0)/sg)+exp(-(v-v0)/sg))

dvdt = lambda n,v: AMP/10.0+(gl*(el-v)+gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*(b*n+a)*(ena-v)+gk*(n/1.3)**4*(ek-v))/cm #mV/ms
dndt = lambda n,v: (ninf(v)-n)/ntau(v) #1/ms



nsteps = len(varray)
nrnarray = range(nsteps)

vprime =varray[0]
fp = open('nhopf.dat','w')
for i in range(nsteps): # copy neuron object
	test.soma.v = vprime
	h.finitialize()
	fp.write("%e  %e\n" % (vprime, test.soma.n_type21))
	vprime += vstep
fp.close()

s = h.IClamp(0.5,sec=test.soma)

s1 = h.IClamp(0.5,sec=test.soma)
s2 = h.IClamp(0.5,sec=test.soma)

s.dur = 1e9
s.delay = 0
s.amp = 1e-3*AMP

s1.dur = 1
s1.delay = 150
s1.amp = 110e-3 #-s.amp - 60, 75

s2.dur = 1
s2.delay = 250
s2.amp = 100e-3 #-s.amp - 150e-3





cvode = h.CVode()
cvode.active()
#h.finitialize()
h.cvode_active(True)
result = h.ref('')
cvode.statename(3,result)
#print result[0]
#quit()

def xdvdt(pstates):
	states = h.Vector(pstates)
	dstates = h.Vector(len(pstates))
	cvode.f(0,states,dstates)
	return dstates[0]

def f(n,v):
	return xdvdt([v,0,0,n])
	
def f2(n,v):	
	c = test.soma
	el = c.el_type21
	gl = c.gl_type21
	b = c.b_type21
	a = c.a_type21
	ena = c.ena_type21
	ek = c.ek_type21
	gna = c.gna_type21
	gk = c.gk_type21
	mhalf = -c.mhalf_type21
	mslope = c.mslope_type21
	s = h.area(0.5,sec=test.soma)
	return AMP/10.0 + gl*(el-v)+gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*(b*n+a)*(ena-v)+gk*(n/1.3)**4*(ek-v)

def g(v):
	try:
		return opt.bisect(f,0,1,args=v)
	except:
		return 1

v0c = []
if VNULL:
	fp = open('vnull_%.1f_hopf.dat' % (AMP),'w')
	h.finitialize()
	for v in varray:
		if v < test.soma.ek_type21:
			continue
		#n0 = opt.bisect(f,0,1,args=v)
		try:
			n0 = opt.bisect(f,0,1,args=v)
			fp.write("%e  %e\n" % (v, n0))
			v0c.append([v,n0])
		except:
			try:
				n0 = opt.bisect(f,-1,0,args=v)
				fp.write("%e  %e\n" % (v, n0))
				v0c.append([v,n0])
			except:
				pass
	fp.close()

if RUN:
#h.finitialize()
	clamp = 0
	test.soma.v = -52.5
	test.soma.ninit_type21 = ninf(-52.5)
	h.finitialize()
	#test.soma.v = -30
	#test.soma.n_type21 = 0.5
	#h.fcurrent()
	h.t = 0
	h.dt = 0.1
	vold = test.soma.v
	#while h.t < 100:
	#	h.fadvance()
	#	vold = test.soma.v
	fp  = open('run_%.2fa_type2.dat' %AMP,'w')
	#while h.t < 0:
	#	h.fadvance()
	while h.t < 50:
		h.fadvance()
	while h.t < 300:
		h.fadvance()
		dv= dvdt(test.soma.v,test.soma.n_type21)
		dw = dndt(test.soma.v,test.soma.n_type21)
		speed = dw**2 + (dv*scale)**2
		fp.write("%e  %e  %e  %e  %e\n" %( h.t, test.soma.v, test.soma.n_type21, -dv*h.area(0.5,sec=test.soma)/100.0, speed))
	fp.close()


if FI:
	irange = arange(0,30,0.1)
	insert(irange,0,0.0)
	s1.amp = 0
	s2.amp = 0
	#print 0, 0
	h.t=0
	h.dt = 0.1
	told = 0
	tspike =0
	for i in irange:
		s.amp = i*1e-3
		test.soma.v = -53.5
		test.soma.ninit_type21 = ninf(-53.5)
		h.finitialize()
		nspike = 0
		f = 0
		while h.t < 100000:
			vold = test.soma.v
			h.fadvance()
			if vold < -20 and test.soma.v > -20:
				t0 = interp(0,[vold,test.soma.v],[h.t-h.dt,h.t])
				told = tspike
				tspike = t0
				nspike +=1
				if nspike > 10:
					break
					f += 1000.0/(tspike-told)
	
				if nspike > 11:
					f = f
					print i, f
					tspike = 0
					told = 0
					h.t = 0
					break
			if h.t - tspike > 10000:
				tspike = 0
				told = 0
				h.t = 0
				print i, 0
				break
		if nspike > 10:
			break
	if BACKWARDS:
		stop = 0
		irange = flip(irange)			
		for i in irange:
			s.amp = i*1e-3
			test.soma.v = -80
			h.finitialize()
			#test.soma.n_type21 = 0
			h.t=0
			nspike = 0
			f = 0
			while h.t < 1000000:
				vold = test.soma.v
				h.fadvance()
				if vold < -20 and test.soma.v > -20:
					t0 = interp(0,[vold,test.soma.v],[h.t-h.dt,h.t])
					told = tspike
					tspike = t0
					nspike +=1
					if nspike > 5:
						f += 1000.0/(tspike-told)
		
					if nspike > 50:
						f = f/(46.0)
						print i, f
						tspike = 0
						told = 0
						h.t = 0
						break
				if h.t - tspike > 10000:
					tspike = 0
					told = 0
					h.t = 0
					stop = 1
					break
			if stop: break

if IV:
	n0 = ninf(varray[0])
	for v in varray:
		if v < 0:
			s.amp = 1.0e-3*AMP
			test.soma.v = v
			h.finitialize()
			dV = xdvdt([test.soma.v,0,0,test.soma.n_type21])
			dV2 = xdvdt([test.soma.v,0,0,n0])
			#I = h.area(0.5, sec=test.soma)*dV + s.amp
			print v, -dV*h.area(0.5,sec=test.soma)/100, -dV2*h.area(0.5,sec=test.soma)/100
		else:
			break


	


if ATTRACTOR:
	un = 0.6
	ln = 0.25
	lv = -62
	uv = -42
	from scipy.spatial import ConvexHull
	points = []
	newpoints = []
	isconvex = False
	for i in range(int(1e6)):
		h.t = 0
		vstart = lv + (uv-lv)*random.random()
		nstart = ln+(un-ln)*random.random()
		test.soma.ninit_type21 = nstart
		test.soma.v = vstart
		h.finitialize(vstart)
		nspike = 0
		vold = vstart
		vmax = -1e9
		lenold = 0
		while h.t < 1000:
			h.fadvance()
			if vold < 0 and test.soma.v > 0:
				nspike = 1
				break
			vold = test.soma.v
			if vold > vmax:
				vmax = vold
		if not nspike:
			newpoints.append([vstart,nstart])
			points.append([vstart,nstart])
		
		#print points
		if len(newpoints)>200:
			pa = array(newpoints)
			if not isconvex:
				hull = ConvexHull(pa,incremental=True)
				isconvex = True
				newpoints = []
			else:
				hull.add_points(pa)
				newpoints = []
	
	for things in hull.vertices:
		print points[things][0], points[things][1]

if BOUNDARY:
	vstep = 1
	cstep = vstep
	search = 0
	vstart = -80
	#varray2 = flip(varray):
	for nstart in arange(0,1,0.01):
		vstart = -80
		test.soma.ninit_type21 = nstart
		test.soma.v = vstart
		h.finitialize()

		search = 0
		nspike = 0
		#nstart = 1e-6
		cstep = vstep
		nspikeold = -1
		while(1):
			if vstart > 0:
				#print nspike
				break			
			h.t = 0
			#h.dt = 0.01
			test.soma.v = vstart
			test.soma.ninit_type21 = nstart
			h.finitialize(vstart)
			nspike = 0
			#print h.t,
			while h.t < 100:
				vold = test.soma.v
				h.fadvance()
				if test.soma.v > 20 and vold < 20:
					nspike+=1
					if nspike > 1:
						break
			#if nspike > 0:
				#print v, nstart, nspike
				
			if nspike != nspikeold and not search:
				if nspikeold < 0:
					nspikeold = nspike
					nspikenew = nspike
					nspikec = nspike
					vstart += cstep
					continue
				else:
					nspikenew = nspike # nspikes at top of range
					nspikec = nspike # current expected spike
					cstep = -cstep/2.0
					vstart = vstart + cstep # go back half a step
					search = 1
					continue

				
			if not search:
				vstart += vstep
				test.soma.v = vstart
				continue

			if nspike != nspikec: # if nspike has changed from last run
				#print v, nstart
				if abs(cstep) < 1e-7:
					search = 0
					#if nspikeold == 1 and nspikenew == 0:
					print vstart, test.soma.ninit_type21, nspikeold, nspikenew
					nspikeold = nspikenew
					nspikec = nspikenew
					break
				else:
					nspikec = nspike
					cstep *= -0.5
					vstart += cstep
				#print v, nstart, cstep
			else:
			   vstart += cstep

				
				
				

		#print v, nstart, nspikeold, cstep, search
				
			
	

if JAC:

	
	
	def vv(v,n):
		temp = 0
		temp += -gl - gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*(b*n+a) - gk*(n/1.3)**4
		temp += 3.0*gna*(1./(1.+exp(-(v+mhalf)/mslope)))**4*(b*n+a)*(ena-v)*exp(-(v+mhalf)/mslope)/mslope #- from argument cancels - from 1/(f(v))^3
		temp *= 1.0/cm # x10 puts dv/dt on the same units as amp
		return temp
	
	def vn(v,n):
		temp =  (4.0/1.3)*(n/1.3)**3*gk*(ek-v) + gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*b
		temp*= 1.0/cm
		return temp
	
	def nv(v,n):
		# dninf/dv / ntau + d/dv (1/(ntau) *ninf(v)
		temp = sn/(2*sl)/(1.0+cosh((v12-v)/sl))/ntau(v)
		temp += 1.0/(t0*sg)*ninf(v)*(exp((v-v0)/sg)-exp(-(v-v0)/sg))
		temp*=1
		return temp
	
	def nn(v,n):
		temp = -1/ntau(v)
		return temp

	from scipy import interpolate, linalg
	VV0c = array(v0c)
	
	
	vint = interpolate.interp1d(VV0c[:,0],VV0c[:,1],kind='cubic')

	
	def zero(v):
		return vint(v) - ninf(v)
	

	crossings = []
	step = 1
	vmin = min(VV0c[:,0])
	vmax = vmin + step
	while(1):

		try:
			roots = opt.bisect(zero,vmin,vmax)
			#print roots, ninf(roots)
			crossings.append([roots,ninf(roots)])
		except:
			vmax += step
			if vmax > max(VV0c[:,0]):
			  break
			else:
				continue
		
		vmin = vmax
		vmax += step
		
		if len(crossings) > 2 or vmax > max(VV0c[:,0]):
			break
	
	# this only works because ncrossings is known to be odd off except with saddle
	for c in crossings:
		v = c[0]
		n = c[1]
		matrix = [[vv(v,n),nv(v,n)],[vn(v,n),nn(v,n)]]
		mat = array(matrix)
		print mat
		eig = linalg.eig(mat)
		print v, n, eig[0][0], eig[0][1]
		
if VF:
	# create 2 columns of vector field at some number of points
	# first 2 values are array dimension
	warray = [0.05+0.1*j for j in range(10)]
	varray = [j for j in range(-80,45,8)]
	wsteps = len(warray)
	vsteps = len(varray)
	fp = open("vfield_%.2f_type2.dat" % AMP,'w')
	fp.write("%d  %d\n" % (vsteps, wsteps))
	for v in varray:
		for w in warray:
			dv = 0.01*c.cm*dvdt(w,v)
			dw = dndt(w,v)
			print v,w, dw, (ninf(v)-w)/ntau(v)
			fp.write("%e  %e  %e  %e\n" %(v, w, dv,dw))
	fp.close()		

quit()
