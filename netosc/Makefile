all: radau5.o decsol.o vc hvc

radau5.o: radau5.f
	gfortran -fno-align-commons -L/usr/local/g2c/lib  -c radau5.f
decsol.o: decsol.f
	gfortran  -fno-align-commons -L/usr/local/g2c/lib  -c decsol.f
vc:  vcdriver.o radau5.o decsol.o redout3.o vc.o netcon.h 
	gfortran  -fno-align-commons -L/usr/local/g2c/lib -o vc vcdriver.o radau5.o decsol.o redout3.o -lm
hvc:  hvcdriver.o radau5.o decsol.o redout3.o hvc.c newcon.h 
	gfortran  -fno-align-commons -L/usr/local/g2c/lib -o hvc hvcdriver.o radau5.o decsol.o redout3.o -lm
vc.o : netcon.h vc.c
	gcc -c vc.c  
hvc.o : newcon.h hvc.c
	gcc -c hvc.c 
hvcdriver.o: newcon.h hvcdriver.c hvc.c
	gcc -c hvcdriver.c
vcdriver.o : netcon.h vcdriver.c vc.c
	gcc -c vcdriver.c
redout3.o : netcon.h redout3.c
	gcc -c redout3.c



