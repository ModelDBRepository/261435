from neuron import h, gui
from numpy import *
from scipy import optimize as opt


RUN = 0
VNULL = 0
VNULL2 =0
CNULL = 0

BOUNDARY = 1
FI = 0

AMP = 62.0

cinit = 1.73
PLOT = 0

minn = 0
sinit =-1

class neuron:
	def __init__(self, cType=0, vinit=None, ninit=None):
		self.soma = h.Section()
		self.soma.L     = 1.
		self.soma.diam  = 1./pi
		self.soma.insert('chav')
		self.soma.insert('cask')
		#self.soma.insert('type21')
		#=== neuron initialization ===#
		#===     set recordings     ===#
		self.v = h.Vector()
		self.v.record(self.soma(0.5)._ref_v)
		self.n = h.Vector()
		self.h = h.Vector()
		self.s = h.Vector()
		self.c = h.Vector()
		self.n.record(self.soma(0.5).chav._ref_n)
		self.h.record(self.soma(0.5).chav._ref_h)
		self.c.record(self.soma(0.5).cask._ref_c)
		self.s.record(self.soma(0.5).cask._ref_s)
		#===    record spike times  ===#
		self.spks	= h.Vector()
		self.sptr	= h.APCount(.5, sec=self.soma)
		self.sptr.thresh = 0.#25
		self.sptr.record(self.spks)
		#==          synapse        ===#


test = neuron(cType=1)


h.finitialize() # sets to type 1 parameters

vstep = 0.1

varray = arange(-80,40,vstep)

nsteps = len(varray)
nrnarray = range(nsteps)

s = h.IClamp(0.5,sec=test.soma)

s1 = h.IClamp(0.5,sec=test.soma)
s2 = h.IClamp(0.5,sec=test.soma)

s.dur = 1e9
s.delay = 0
s.amp = 1e-6*AMP

s1.dur = 1
s1.delay = 150
s1.amp = 0e-3 #-s.amp - 150e-3

s2.dur = 1
s2.delay = 250
s2.amp = 0e-3 #-s.amp - 150e-3


test.soma.shalf_cask = -30.0
test.soma.stau_cask = 10.0e-1
test.soma.km_cask = 2.0
test.soma.gca_cask = 0.0005e-2
test.soma.gsk_cask = 0.03e-2
test.soma.gl_chav = 0.003e-2
test.soma.sslope_cask = 2.0
test.soma.RHO_cask = 0.00012*10
test.soma.kc_cask = 0.4
test.soma.SC1_chav = 0.1 #n
test.soma.SC2_chav = 0.1 #h

test.soma.gna_chav = 3.0e-2
test.soma.gk_chav = 0.8e-2

test.soma.cm = 1

h.dt = 0.1
cvode = h.CVode()
cvode.active()
#h.finitialize()
h.cvode_active(True)
result = h.ref('')
cvode.statename(3,result)

#h.secondorder = 0
#cvode.fixed_step()
#cvode.statename(3,result)
#print result[0]
#quit()

#print h.area(0.5,sec=test.soma)
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
if RUN:
#h.finitialize()
	clamp = 0
	test.soma.v = -80
	#test.soma.ninit_type21 = 0.34645
	h.finitialize()
	#test.soma.v = -30
	#test.soma.n_type21 = 0.5
	#h.fcurrent()
	h.t = 0
	h.dt = 0.1
	vold = test.soma.v
	while h.t < 500:
		h.fadvance()
		vold = test.soma.v
	h.frecord_init()
	while h.t < 3000:
		h.fadvance()
		if not PLOT:
			print h.t, test.soma.v, test.soma.c_cask, test.soma.s_cask, test.soma.n_chav, test.soma.h_chav
		vold = test.soma.v
	
	if PLOT:
		import glob
		fig = plt.figure()
		ax = plt.axes(projection='3d')
		zline = test.v.to_python()#voltage
		yline = test.h.to_python()
		xline = test.c.to_python()
		xp = []
		yp = []
		zp = []
		xm = []
		ym = []
		zm = []
		flist = glob.glob('t3_separatrix_*b.dat')
		flist.sort()
		for file in flist:
			try:
				xs,ys,zs=loadtxt(file,unpack=True,usecols=[2,1,0])
				print zs[0], zs[-1]
				ax.plot3D(xs,ys,zs,'black')
			except:
				print file
				continue
			xp.append(xs[0])
			yp.append(ys[0])
			zp.append(zs[0])
			xm.append(xs[-1])
			ym.append(ys[-1])
			zm.append(zs[-1])
		xm = asarray(xm)
		ym = asarray(ym)
		zm = asarray(zm)
		xp = asarray(xp)
		yp = asarray(yp)
		zp = asarray(zp)
		print size(xm), size(ym), size(zp)
		ax.plot3D(xm,ym,zm,'black')
		ax.plot3D(xp,yp,zp,'black')
		ax.plot3D(xline,yline,zline,'red')
		ax.view_init(150,30)
		plt.show()


indices = []

		
def dvdt(c,v):
	states = h.Vector()

	test.soma.v = v
	test.soma.cinit_cask = c
	test.soma.sinit_cask = sinit
	h.finitialize()
	cvode.states(states)
	dstates = h.Vector(len(states))
	cvode.f(0,states,dstates)
	return dstates[0],dstates[2],dstates[4]
	
def dvdt2(hh,v):
	states = h.Vector()
	test.soma.v = v
	test.soma.hinit_chav = hh
	h.finitialize()
	#test.soma.h_chav = hh
	cvode.states(states)
	dstates = h.Vector(len(states))
	cvode.f(0,states,dstates)
	return dstates[0]

def f(c,v):
	return dvdt(c,v)[0]
	
def g(c,v):
	return dvdt(c,v)[1]
	
def p(hh,v):
	return dvdt2(hh,v)

if VNULL:
	states = h.Vector()
	cvode.states(states)
	fp = open('vnull_%.3f_burst3.dat'%sinit,'w')
	for v in varray:
		#print 0,v,f(0,v)
		#print 1,v,f(1,v)
		try:
			c0 = opt.bisect(f,0,10,args=v)
			fp.write("%e  %e   %e\n" % (v, c0, test.soma.h_chav))
			v0c.append([v,c0])
		except:
			pass
	fp.close()
	
if VNULL2:
	states = h.Vector()
	cvode.states(states)
	fp = open('vnull_h_%.3f_burst3.dat'%cinit,'w')
	test.soma.cinit_cask = cinit
	test.soma.sinit_cask = -1
	test.soma.RHO_cask = 0
	for v in varray:
		#print 0,v,p(0.1,v), test.soma.h_chav, test.soma.hinit_chav
		#print 1,v,p(0.9,v), test.soma.h_chav, test.soma.hinit_chav
		
		try:
			h0 = opt.bisect(p,0.01,0.99,args=v)
			fp.write("%e  %e   %e\n" % (v, h0, cinit))
		except:
			pass
	fp.close()
	
if CNULL:
	carray = arange(0,10,0.01)
	states = h.Vector()
	cvode.states(states)
	fp = open('cnull_burst3.dat','w')
	for v in varray:
		#print 0,v,g(0,v)
		#print 1,v,g(1,v)
		try:
			c0 = opt.bisect(g,0,10,args=v)
			fp.write("%e  %e\n" % (v, c0))
			v0c.append([v,c0])
		except:
			pass
	fp.close()
	
if FI:
	irange = arange(0,5,0.01)
	test.soma.gsk_cask = 0 # remove slow current to isolate F-I curve
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
	

if BOUNDARY:
	search = 0
	vstart = -20
	test.soma.RHO_cask = 1e-9 # fix calcium at cinit
	for c in arange(3,4,0.01):
		test.soma.cinit_cask = c # choose points near bifurcations:
		h.finitialize()
		vold = 100
		nspike = 0
		while h.t < 2000:
			h.fadvance()
			if vold < 0 and test.soma.v > 0:
				nspike+=1
				if nspike > 1:
					print c
					break
			vold = test.soma.v
		

quit()
