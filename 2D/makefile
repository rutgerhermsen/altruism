# Start of the makefile
# Defining variables
objects = multiscale.o price.o state.o g_of_r.o kernels.o fft.o utils.o data_structures.o parameters.o shorthands.o pairwise.o ziggurat.o fftpack5.1d.o altruism_default_2D.o
f90comp = gfortran
switch = -O3 -march=native -Wall -Wextra
# Makefile
altr2D.out: $(objects)
	$(f90comp) -o altr2D.out $(switch) $(objects)
altruism_default_2D.o: multiscale.o price.o state.o g_of_r.o kernels.o fft.o utils.o data_structures.o parameters.o shorthands.o pairwise.o ziggurat.o fftpack5.1d.o altruism_default_2D.f90
	$(f90comp) -c $(switch) altruism_default_2D.f90
fftpack5.1d.mod: fftpack5.1d.o fftpack5.1d.f90
	$(f90comp) -c $(switch) fftpack5.1d.f90
fftpack5.1d.o: fftpack5.1d.f90
	$(f90comp) -c $(switch) fftpack5.1d.f90
ziggurat.mod: ziggurat.o ziggurat.f90
	$(f90comp) -c $(switch) ziggurat.f90
ziggurat.o: ziggurat.f90
	$(f90comp) -c $(switch) ziggurat.f90
pairwise.mod: pairwise.o pairwise.f90
	$(f90comp) -c $(switch) pairwise.f90
pairwise.o: pairwise.f90
	$(f90comp) -c $(switch) pairwise.f90
parameters.mod: parameters.o parameters.f90
	$(f90comp) -c $(switch) parameters.f90
parameters.o: parameters.f90
	$(f90comp) -c $(switch) parameters.f90
shorthands.mod: shorthands.o shorthands.f90
	$(f90comp) -c $(switch) shorthands.f90
shorthands.o: parameters.o shorthands.f90
	$(f90comp) -c $(switch) shorthands.f90
data_structures.mod: data_structures.o data_structures.f90
	$(f90comp) -c $(switch) data_structures.f90
data_structures.o: parameters.o data_structures.f90
	$(f90comp) -c $(switch) data_structures.f90
utils.mod: utils.o utils.f90
	$(f90comp) -c $(switch) utils.f90
utils.o: pairwise.o parameters.o utils.f90
	$(f90comp) -c $(switch) utils.f90
fft.mod: fft.o fft.f90
	$(f90comp) -c $(switch) fft.f90
fft.o: parameters.o fft.f90
	$(f90comp) -c $(switch) fft.f90
kernels.mod: kernels.o kernels.f90
	$(f90comp) -c $(switch) kernels.f90
kernels.o: pairwise.o fft.o shorthands.o parameters.o kernels.f90
	$(f90comp) -c $(switch) kernels.f90
g_of_r.mod: g_of_r.o g_of_r.f90
	$(f90comp) -c $(switch) g_of_r.f90
g_of_r.o: fft.o shorthands.o parameters.o g_of_r.f90
	$(f90comp) -c $(switch) g_of_r.f90
state.mod: state.o state.f90
	$(f90comp) -c $(switch) state.f90
state.o: shorthands.o parameters.o data_structures.o state.f90
	$(f90comp) -c $(switch) state.f90
price.mod: price.o price.f90
	$(f90comp) -c $(switch) price.f90
price.o: utils.o state.o shorthands.o parameters.o data_structures.o price.f90
	$(f90comp) -c $(switch) price.f90
multiscale.mod: multiscale.o multiscale.f90
	$(f90comp) -c $(switch) price.f90
multiscale.o: price.o fft.o kernels.o pairwise.o shorthands.o parameters.o data_structures.o multiscale.f90
	$(f90comp) -c $(switch) multiscale.f90
		
# Cleaning everything
clean:
# End of the makefile
