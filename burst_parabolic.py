from neuron import h, gui
from numpy import *
from scipy import optimize as opt
from scipy import interpolate as inter


RUN = 1
VNULL = 0
CNULL = 0
VSNULL = 0

SWEEP=0
BOUNDARY=1
vstart = -50

PLOT=0
if PLOT==1:
	RUN=1
	VSNULL=0
	BOUNDARY =0
AMP = -0.025

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

varray = arange(-70,0,vstep)

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

# factor of 10 throughout compresses time such that an action potential
# has a physiologically realistic width.
test.soma.stau_cask = 9500e-1
test.soma.RHO_cask = 0.00015*10
test.soma.kc_cask = 0.00425
test.soma.SC1_chav = 0.1 #n
test.soma.SC2_chav = 0.1 #h

test.soma.gsk_cask = 0.03e-2
test.soma.gca_cask = 0.004e-2

test.soma.cinit_cask = 0.2
test.soma.sinit_cask = 0.1

test.soma.gna_chav = 0*2.5e-2
test.soma.gk_chav = 0.3e-2
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

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.tri as mtri


def dvdt2(c,v,s):
	states = h.Vector()

	test.soma.v = v
	test.soma.cinit_cask = c
	test.soma.sinit_cask = s
	dstates = h.Vector()
	h.finitialize()
	cvode.states(states)
	cvode.f(0,states,dstates)
	return dstates[0],dstates[2]

if RUN:
#h.finitialize()


	clamp = 0
	test.soma.v = vstart
	#test.soma.ninit_type21 = 0.34645
	h.finitialize()
	#test.soma.v = -30
	#test.soma.n_type21 = 0.5
	#h.fcurrent()
	h.t = 0
	h.dt = 0.1
	vold = test.soma.v
	states = h.Vector()
	dstates = h.Vector()
	while h.t < 5000:
		h.fadvance()
		vold = test.soma.v
	h.frecord_init()
	while h.t < 20000:
		h.fadvance()
		if not PLOT:
			cvode.states(states)
			cvode.f(0,states,dstates)
			print h.t, test.soma.v, test.soma.c_cask, test.soma.s_cask, test.soma.n_chav, test.soma.h_chav, dstates.x[2]
		vold = test.soma.v
	
	if PLOT:
		#X1,Y1,Z1 = loadtxt("cnull_burst2.dat",unpack=True)
		Y1,Z1 = loadtxt("boundary-30.dat",unpack=True,usecols=[0,1]) # s, c
		Y2,Z2 = loadtxt("boundary-50.dat",unpack=True,usecols=[0,1]) # s, c
		
		
		YT = []
		ZT = []
		for i in range(len(Y1)):
			if Z1[i] < 0.2:
				ZT.append(Z1[i])
				YT.append(Y1[i])
				
		YT2 = []
		ZT2 = []
		for i in range(len(Y2)):
			if Z2[i] < 0.2:
				ZT2.append(Z2[i])
				YT2.append(Y2[i])
				

		sz = len(YT)
		sz2 = len(YT2)
		X1 = zeros(sz)
		X2 = zeros(sz2)

		
		
		
		a1 = inter.interp1d(YT,ZT) # c,s
		a2 = inter.interp1d(YT2,ZT2)
	
		

		
		fig = plt.figure()
		ax = plt.axes(projection='3d')


		zline = test.v.to_python()#voltage
		yline = test.s.to_python()
		xline = test.c.to_python()
		
		dold = 0
		z1 = []
		z2 = []
		
		y1 = []
		y2 = []
		
		x1 = []
		x2 = []
		
		dold = dvdt2(xline[0],zline[0],yline[0])[1]
		print dold
		
		
		crossings = 0
		scatters = []
		for i in range(len(zline)-1):
			vtemp = min(zline[i+1],-50)
			dnew = xline[i] - a1(yline[i])
			if crossings <= 2:
				z1.append(zline[i+1])
				x1.append(xline[i+1])
				y1.append(yline[i+1])
			if dnew*dold < 0:
				if crossings > 0:
					scatters.append([xline[i],yline[i],zline[i]])
				print dnew
				crossings +=1
			dold = dnew
			if crossings >= 3:
				z2.append(zline[i+1]) # v
				x2.append(xline[i+1]) # c
				y2.append(yline[i+1]) # s
			if crossings > 5:
				break
		

		# is c(s)

		z1val = -70
		z2val = 40
		
		#a = scatters.pop()
		XS = full_like(X1, z1val)
		X2 = full_like(X2, z2val)
		X3 = full_like(X2, z1val)
		X4 = full_like(X1,z2val)
		
		print XS
		
		
		print len(ZT), len(YT), len(XS)
		
		x3 = x1+x2
		y3 = y1+y2
		z3 = z1+z2
			
		ax.set_zlim([-70,40])
		ax.set_ylim([0.0,0.4])
		ax.set_xlim([0.0,0.2])
		xs = [s[0] for s in scatters]
		ys = [s[1] for s in scatters]
		zs = [s[2] for s in scatters]

		ax.plot3D(x3,y3,z3,'red',linewidth=2)
		#ax.plot3D(x2,y2,z2,'green',linewidth=2)
		#ax.plot3D(x1,y1,z1,'red',linewidth=2)
		ax.plot3D(ZT,YT,XS,'blue',linewidth=2)
		#ax.plot3D(ZT2,YT2,X2,'green',linewidth=2)
		ax.plot3D(ZT,YT,X4,'blue',linewidth=2)
		ax.scatter(xs,ys,zs)
		#ax.plot3D(ZT2,YT2,X3,'green',linewidth=2)
		#ax.plot_surface(YMESH,ZMESH,XMESH)	


		ax.set_xlabel('c')
		ax.set_ylabel('s')
		ax.set_zlabel('v')
		from mpl_toolkits.mplot3d import proj3d
		def orthogonal_proj(zfront, zback):
				a = (zfront+zback)/(zfront-zback)
				b = -2*(zfront*zback)/(zfront-zback)
				# -0.0001 added for numerical stability as suggested in:
				# http://stackoverflow.com/questions/23840756
				return numpy.array([[1,0,0,0],
									[0,1,0,0],
									[0,0,a,b],
									[0,0,-0.0001,zback]])
	 
	# Later in your plotting code ...

		ax.view_init(60,-25)
		proj3d.persp_transformation = orthogonal_proj
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
	return dstates[0],dstates[2]
	


def f(c,v):
	return dvdt(c,v)[0]
	
def g(c,v):
	return dvdt(c,v)[1]
	
def q(c,v,s):
	return dvdt2(c,v,s)[0]

if VNULL:
	states = h.Vector()
	cvode.states(states)
	fp = open('vnull_%.3f_burst2.dat'%sinit,'w')
	for v in varray:
		#print 0,v,f(0,v)
		#print 1,v,f(1,v)
		try:
			c0 = opt.bisect(f,0,1,args=v)
			test.soma.v = v
			test.soma.cinit_cask = c0
			h.finitialize()
			s = test.soma.s_cask
			fp.write("%e  %e  %e\n" % (v, c0, s))
			v0c.append([v,c0])
		except:
			pass
	fp.close()
	
if CNULL:
	carray = arange(0,1,0.01)
	states = h.Vector()
	cvode.states(states)
	sinit = -1
	fp = open('cnull_burst2.dat','w')
	for v in varray:
		test.soma.v = v
		h.finitialize()
		#print 0,v,g(0,v)
		#print 1,v,g(1,v)
		try:
			c0 = opt.bisect(g,0,1,args=v)
			test.soma.v = v
			test.soma.cinit_cask = c0
			h.finitialize()
			s = test.soma.s_cask
			#if c0 < 0.3:
			fp.write("%e  %e  %e\n" % (v, c0, s))
			v0c.append([v,c0])
		except:
			pass
	fp.close()
	
if VSNULL:
	sarray = arange(0,0.45,0.01)
	states = h.Vector()
	cvode.states(states)
	dstates = h.Vector(len(states))
	fp = open('nullsurface_burst2.dat','w')
	for s in sarray:
		for v in varray:
			try:
				c0 = opt.bisect(q,0.02,0.25,args=(v,s))
				fp.write("%e  %e  %e\n" % (v, c0, s))
				#print v,c0,s
			except:
				pass
	fp.close()
	
if SWEEP:
	#test.soma.RHO_cask = 1e-7
	#test.soma.stau_cask = 1e9 # turn off slow dynamics
	for cstart in arange(0.01,0.5,0.01):
		for astart in arange(0.01,0.5,0.01):
			test.soma.sinit_cask = astart
			test.soma.cinit_cask = cstart
			test.soma.v = vstart
			nspike = 0
			h.finitialize()
			vold = 100
			vmax = -1e9
			#while h.t < 100:
			#	h.fadvance() # skip transients?
			while h.t < 2000:
				h.fadvance()
				vold = test.soma.v
				if vold > vmax:
					vmax = vold
				h.fadvance()
				if test.soma.v > 0 and vold < 0:
					nspike+=1
					if nspike > 0:
						print test.soma.sinit_cask, test.soma.cinit_cask, h.t, vmax
						break

		
if BOUNDARY:
	astep = -0.1
	cstep = astep
	search = 0
	vstartprime = vstart
	test.soma.RHO_cask = 1e-9
	test.soma.stau_cask = 1e9 # turn off slow dynamics
	for cstart in arange(0.01,0.5,0.005):
		vstart = vstartprime
		astart = 0.5
		test.soma.sinit_cask = astart
		test.soma.cinit_cask = cstart
		test.soma.v = vstart
		h.finitialize()
		#while h.t < 100:
		#	h.fadvance() # skip transients?


		search = 0
		nspike = 0
		#nstart = 1e-6
		cstep = astep
		nspikeold = -1
		nspikec = -1
		while(1):
			if astart > 1.0 or astart < 0:
				break
			#else:
			#	print cstart, astart, nspike, nspikec		
			h.t = 0
			#h.dt = 0.01
			test.soma.v = vstart
			test.soma.cinit_cask = astart
			test.soma.sinit_cask = cstart
			h.finitialize()
			nspike = 0
			#print h.t,
			while h.t < 100:
				h.fadvance() # skip transients?
			while h.t < 500:
				h.fadvance()
				vold = test.soma.v
				h.fadvance()
				if test.soma.v > 0 and vold < 0:
					nspike+=1
					if nspike > 0:
						break
				
			if nspike != nspikeold and not search:
				if nspikeold < 0:
					nspikeold = nspike
					nspikenew = nspike
					nspikec = nspike
					astart += cstep
					continue
				else:
					nspikenew = nspike # nspikes at top of range
					nspikec = nspike # current expected spike
					cstep = -cstep/2.0
					astart = astart + cstep # go back half a step
					search = 1
					continue

				
			if not search:
				astart += astep
				continue

			if nspike != nspikec: # if nspike has changed from last run
				#print v, nstart
				if abs(cstep) < 1e-7:
					search = 0
					if nspikeold > 0 and nspikenew == 0:
						print test.soma.sinit_cask, test.soma.cinit_cask, nspikeold, nspikenew
						break
					else:
						print test.soma.sinit_cask, test.soma.cinit_cask, nspikeold, nspikenew
						break
					nspikeold = nspikenew
					nspikec = nspikenew
					cstep = astep
					astart+=cstep
				else:
					nspikec = nspike
					cstep *= -0.5
					astart += cstep
				#print v, nstart, cstep
			else:
			   astart += cstep

quit()
