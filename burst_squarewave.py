from neuron import h, gui
from numpy import *
from scipy import optimize as opt


RUN = 1
VNULL = 1
CNULL = 1
BOUNDARY = 0
AMP = 0.0

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
		self.c = h.Vector()
		self.n.record(self.soma(0.5).chav._ref_n)
		self.h.record(self.soma(0.5).chav._ref_h)
		self.h.record(self.soma(0.5).cask._ref_c)
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
s.amp = 1e-3*AMP

s1.dur = 1
s1.delay = 150
s1.amp = 0e-3 #-s.amp - 150e-3

s2.dur = 1
s2.delay = 250
s2.amp = 0e-3 #-s.amp - 150e-3


 #d
test.soma.shalf_cask = -50.0
test.soma.stau_cask = 1e-9 #1 
test.soma.km_cask = 0.5
test.soma.gca_cask = 0.004e-2
test.soma.gsk_cask = 0.03e-2
test.soma.sslope_cask = 1.0/0.15
test.soma.RHO_cask = 0.00025*10e6
test.soma.kc_cask = 0.0085
test.soma.SC1_chav = 0.7e-3 #n 0.03b
test.soma.SC2_chav = 0.7e-3 #h 0.03b
test.soma.gna_chav = 2.0e-2
test.soma.gk_chav = 0.24e-2
test.soma.gl_chav = 0.003e-2

test.soma.cm = 1


cvode = h.CVode()
cvode.active()
#h.finitialize()
h.cvode_active(True)
result = h.ref('')
cvode.statename(3,result)
#print result[0]
#quit()

#print h.area(0.5,sec=test.soma)

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
	while h.t < 2000:
		h.fadvance()
		vold = test.soma.v
	while h.t < 4000:
		h.fadvance()
		print h.t, test.soma.v, test.soma.c_cask, test.soma.s_cask, test.soma.n_chav, test.soma.h_chav
		vold = test.soma.v


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
	return dstates[0],dstates[2]

def f(c,v):
	return dvdt(c,v)[0]
	
def g(c,v):
	return dvdt(c,v)[1]

if VNULL:
	states = h.Vector()
	cvode.states(states)
	fp = open('vnull_burst1.dat','w')
	for v in varray:
		#print 0,v,f(0,v)
		#print 1,v,f(1,v)
		try:
			c0 = opt.bisect(f,0,2,args=v)
			fp.write("%e  %e\n" % (v, c0))
			v0c.append([v,c0])
		except:
			pass
	fp.close()
	
if CNULL:
	carray = arange(0,10,0.01)
	states = h.Vector()
	cvode.states(states)
	fp = open('cnull_burst1.dat','w')
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
	
if BOUNDARY:
	search = 0
	vstart = -40
	test.soma.RHO_cask = 1e-9
	#test.soma.stau_cask = 1e9 # turn off slow dynamics
	for cstart in arange(0.2,1,0.005):
		vstart = -45
		test.soma.cinit_cask = cstart
		test.soma.v = vstart
		h.t = 0
		h.finitialize()
		vold = vstart
		nspike = 0
		while h.t < 1000:
			h.fadvance()
			if vold < 0 and test.soma.v > 0:
				#print h.t
				nspike += 1
			if nspike > 2:
				print nspike, cstart, test.soma.c_cask
				break
			vold = test.soma.v
		if nspike < 2:
			print nspike, cstart
			break
			

quit()
