from neuron import h, gui
from numpy import *
from scipy import optimize as opt


FI = 1 # for figure 1b,d,e
RUN = 1
VNULL = 0
IV = 0
BACKWARDS=1
VF = 0
JAC = 0
BOUNDARY = 0
GHOST2 = 0 #uses backwards integration from saddle to find seperatrix for bimodal SN example


if JAC or BOUNDARY or GHOST2:
	VNULL = 1
	
t0 =5.0 # 5.0 for SNIC, 0.4 for SN


AMP = 0.0 # saddle destruction at 2.29

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
test.soma.el_type21 = -65.0
test.soma.gl_type21 = 0.3
test.soma.v12_type21 = -45
test.soma.sl_type21 = 0.5/0.065
test.soma.gk_type21=36
test.soma.gna_type21 = 120
test.soma.ninit_type21 = -1 # use ninf instead

test.soma.mhalf_type21 = -33.0
test.soma.mslope_type21 = 0.5/0.055
test.soma.ninit_type21 = -1 # use ninf instead

#t0 = 0.4 # 0.4 5.0

test.soma.v0_type21 = test.soma.v12_type21
test.soma.sg_type21 = 2*test.soma.sl_type21
test.soma.t0_type21 = t0
test.soma.n_type21 = 0.1
test.soma.v = -80
h.finitialize()
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
#s = h.area(0.5,sec=test.soma)
n0 = 0
sn = 1

ninf = lambda v: n0 + sn/(1.+exp(-(v-v12)/sl ))
ntau = lambda v: t0/(exp((v-v0)/sg)+exp(-(v-v0)/sg))
dvdt = lambda n,v: gl*(el-v)+gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*(b*n+a)*(ena-v)+gk*(n/1.3)**4*(ek-v)+ 1000.0*s.amp*100.0/h.area(0.5,sec=c)
dndt = lambda n,v: (ninf(v)-n)/ntau(v)

vstep = 0.03

varray = arange(-80,50,vstep)

nsteps = len(varray)
nrnarray = range(nsteps)

vprime =varray[0]
fp = open('nsnic_carmen.dat','w')
narray = []
for i in range(nsteps): # copy neuron object
	test.soma.v = vprime
	h.finitialize()
	fp.write("%e  %e\n" % (vprime, ninf(vprime)))
	vprime += vstep
fp.close()

s = h.IClamp(0.5,sec=test.soma)

s1 = h.IClamp(0.5,sec=test.soma)
s2 = h.IClamp(0.5,sec=test.soma)

s.dur = 1e9
s.delay = 0
s.amp = 1.0e-3*AMP

s1.dur = 2
s1.delay = 1000
s1.amp = 0*3e-2

s2.dur = 2
s2.delay = 1100
s2.amp = 0*(-3e-2)

#s1.amp = 2.34e-3
#s2.amp = 0e-3

scale = 1.0/(test.soma.ena_type21 - test.soma.ek_type21)


cvode = h.CVode()
cvode.active()
#h.finitialize()
h.cvode_active(True)
result = h.ref('')
cvode.statename(3,result)
#print result[0]
#quit()


def hdvdt(pstates):
	states = h.Vector(pstates)
	dstates = h.Vector(len(pstates))
	cvode.f(0,states,dstates)
	return dstates[0], dstates[1]

def boltz(v,a,b):
	return 1.0/(1.0+exp(-(v-a)/b))


def f2(n,v):
	return hdvdt([v,0,0,n])

def f(n,v):
	h.finitialize()
	test.soma.v = v
	test.soma.n_type21 = n
	a = test.soma
	edvdt = a.gna_type21*(boltz(v,a.mhalf_type21, a.mslope_type21)**3)*(1.0-n)*(v-a.ena_type21)+a.gl_type21*(v-a.el_type21)-1000*s.amp*100.0/h.area(0.5,sec=a)
	if n <0 :
		edvdt -= a.gk_type21*((n/1.3)**4.0)*(v-a.ek_type21)
	else:
		edvdt += a.gk_type21*((n/1.3)**4.0)*(v-a.ek_type21)
	#return dvdt([v,0,0,n])
	return edvdt
	
def g(v):
	try:
		return opt.bisect(f,0,1,args=v)
	except:
		try:
			return opt.bisect(f,-1,0,args=v)
		except:
			return 1

v0c = []
if VNULL:
	fp = open('vsnic_%.2f.dat' % AMP,'w')
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
	h.t = 0
	h.dt = 0.1
	test.soma.v = -58
	test.soma.ninit_type21 = -1
	h.finitialize()
	fp  = open('run_%.2f_type1.dat' %AMP,'w')
	#while h.t < 0:
	#	h.fadvance()
	#while h.t < 50:
	#	h.fadvance()
	while h.t < 100:
		h.fadvance()
		dv= dvdt(test.soma.v,test.soma.n_type21)
		dw = dndt(test.soma.v,test.soma.n_type21)
		speed = dw**2 + (dv*scale)**2
		fp.write("%e  %e  %e  %e  %e\n" %( h.t, test.soma.v, test.soma.n_type21, -dv*h.area(0.5,sec=test.soma)/100.0, speed))
	fp.close()

if FI:
	irange = arange(0,5,0.01)
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
		test.soma.v = -80
		h.finitialize()
		nspike = 0
		f = 0
		while h.t < 100000:
			vold = test.soma.v
			h.fadvance()
			if vold < 0 and test.soma.v > 0:
				t0 = interp(0,[vold,test.soma.v],[h.t-h.dt,h.t])
				told = tspike
				tspike = t0
				nspike +=1
				if nspike > 10:
					break
					f += 1000.0/(tspike-told)
	
				if nspike > 11:
					f = f/(51.0)
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
		irange = flip(irange)			
		for i in irange:
			s.amp = i*1e-3
			test.soma.v = 40
			h.finitialize()
			#test.soma.n_type21 = 0
			h.t=0
			nspike = 0
			f = 0
			while h.t < 100000:
				vold = test.soma.v
				h.fadvance()
				if vold < 0 and test.soma.v > 0:
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
					print i, 0
					break

if IV:
	for v in varray:
		if v < 0:
			test.soma.v = v
			h.finitialize()
			dV,dw = dvdt([test.soma.v,0,0,test.soma.n_type21])
			dV2,dw2 = dvdt([test.soma.v,0,0,ninf(-60)])
			#I = h.area(0.5, sec=test.soma)*dV + s.amp
			print v, -dV*h.area(0.5,sec=test.soma)/100
		else:
			break
	
	
if JAC:
	
	def vv(v,n):
		temp = 0
		temp += -gl - gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*(b*n+a) - gk*(n/1.3)**4
		temp += 3*gna*(1./(1.+exp(-(v+mhalf)/mslope)))**4*(b*n+a)*(ena-v)*exp(-(v+mhalf)/mslope)/mslope
		return temp
	
	def vn(v,n):
		temp = 1.0/1.3*4.0*(n/1.3)**3*gk*(ek-v) + gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*b
		return temp
	
	def nv(v,n):
		# dninf/dv / ntau + d/dv (1/(ntau) *ninf(v)
		temp = (sn/((1.+exp(-(v-v12)/sl ))**2)*exp(-(v-v12)/sl )/sl)/ntau(v)
		temp += 1.0/(t0*sg)*ninf(v)*(exp((v-v0)/sg)-exp(-(v-v0)/sg))
		return temp
	
	def nn(v,n):
		temp = -1.0/ntau(v)
		return temp

	from scipy import interpolate, linalg
	VV0c = array(v0c)
	
	
	vint = interpolate.interp1d(VV0c[:,0],VV0c[:,1],kind='cubic')

	
	def zero(v):
		return vint(v) - ninf(v)
	


	
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
		matrix = [[vv(v,n),vn(v,n)],[nv(v,n),nn(v,n)]]
		mat = array(matrix)
		eig = linalg.eig(mat)
		print v, n, eig[0][0], eig[0][1]

if VF:
	# create 2 columns of vector field at some number of points
	# first 2 values are array dimension
	warray = [0.05+0.05*j for j in range(20)]
	varray = [j+2.5 for j in range(-80,45,4)]
	wsteps = len(warray)
	vsteps = len(varray)
	fp = open("vfield_%.2f_type1.dat" % AMP,'w')
	fp.write("%d  %d\n" % (vsteps, wsteps))
	for v in varray:
		for w in warray:
			dv = 0.01*c.cm*dvdt(w,v)
			dw = dndt(w,v)
			print v,w, dw, (ninf(v)-w)/ntau(v)
			fp.write("%e  %e  %e  %e\n" %(v, w, dv,dw))
	fp.close()

if GHOST2:
	#import numpy as np
	#from scipy import optimize
	from scipy import integrate, interpolate
	import scipy.stats as sps
	h.finitialize()
	v0c = array(v0c)
	vXn = v0c[ where(v0c[:,0] < -40.) ]
	VV0c = array(vXn)

	vint = interpolate.interp1d(VV0c[:,0],VV0c[:,1],kind='cubic')
	c = test.soma
	
	def zero(v):
		return vint(v) - ninf(v)
	
	vmin = min(VV0c[:,0])
	vmax = max(VV0c[:,0])
	res = opt.minimize(zero,-60,method='nelder-mead',)
	res2 = opt.minimize(vint,-60,method='nelder-mead',)
	#print res
	vnX = [float(res2.x), vint(float(res2.x))]

	#ntau = lambda v: t0/(exp((v-v0)/sg)+exp(-(v-v0)/sg))
	#dvdt = lambda n,v: -AMP*1e-1/s + gl*(el-v)+gna*(1./(1.+exp(-(v+mhalf)/mslope)))**3*(b*n+a)*(ena-v)+gk*(n/1.3)**4*(ek-v)
	
	def rhs(t,Y): return [-dvdt(Y[1],Y[0]),-(ninf(Y[0])-Y[1])/ntau(Y[0])]
	print vnX[0], vnX[1]
	slv = integrate.ode(rhs).set_initial_value(vnX, 0)#.set_integrator('zvode', method='bdf')
	thc = [ slv.integrate(slv.t+0.001) ]
	vbl,vbr = -80.0,vnX[0]+100
	i=0
	while slv.successful() and slv.t < 100. and vbl < thc[-1][0] < vbr and 0. < thc[-1][1] < 1:
		thc.append( slv.integrate(slv.t+0.001) )
		i+=1
		if i%10 ==0:
			print thc[-1][0], thc[-1][1]
	thc = array(thc)
	
	
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
			if vstart > -30:
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
				if test.soma.v > 0 and vold < 0:
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
				
			
				

					
	
quit()
