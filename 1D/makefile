# Start of the makefile
# Defining variables
objects = pairwise.o ziggurat.o fftpack5.1d.o parameters.o data_structures.o utils.o shorthands.o fft.o kernels.o g_of_r.o state.o price.o hamilton.o multiscale.o collectives.o altruism_default_1D.o
f90comp = gfortran
switch = -O3 -march=native -Wall -Wextra
# Makefile
altr.out: $(objects)
	$(f90comp) -o altr.out $(switch) $(objects)
altruism_default_1D.o: pairwise.o ziggurat.o fftpack5.1d.o parameters.o data_structures.o utils.o shorthands.o fft.o kernels.o g_of_r.o state.o price.o multiscale.o collectives.o altruism_default_1D.f90
	$(f90comp) -c $(switch) altruism_default_1D.f90
collectives.mod: collectives.o collectives.f90
	$(f90comp) -c $(switch) collectives.f90
collectives.o: parameters.o shorthands.o data_structures.o fft.o kernels.o state.o utils.o collectives.f90
	$(f90comp) -c $(switch) collectives.f90
multiscale.mod: multiscale.o multiscale.f90
	$(f90comp) -c $(switch) multiscale.f90
multiscale.o: parameters.o shorthands.o pairwise.o kernels.o data_structures.o fft.o price.o multiscale.f90
	$(f90comp) -c $(switch) multiscale.f90
price.mod: price.o price.f90
	$(f90comp) -c $(switch) price.f90
price.o: parameters.o shorthands.o data_structures.o state.o utils.o price.f90
	$(f90comp) -c $(switch) price.f90
state.mod: state.o state.f90
	$(f90comp) -c $(switch) state.f90
state.o: parameters.o shorthands.o data_structures.o utils.o state.f90
	$(f90comp) -c $(switch) state.f90
g_of_r.mod: g_of_r.o g_of_r.f90
	$(f90comp) -c $(switch) g_of_r.f90
g_of_r.o: parameters.o shorthands.o fft.o g_of_r.f90
	$(f90comp) -c $(switch) g_of_r.f90
kernels.mod: kernels.o kernels.f90
	$(f90comp) -c $(switch) kernels.f90
kernels.o: parameters.o shorthands.o fft.o kernels.f90
	$(f90comp) -c $(switch) kernels.f90
fft.mod: fft.o fft.f90
	$(f90comp) -c $(switch) fft.f90
fft.o: parameters.o  fft.f90
	$(f90comp) -c $(switch) fft.f90
shorthands.mod: shorthands.o shorthands.f90
	$(f90comp) -c $(switch) shorthands.f90
shorthands.o: parameters.o  shorthands.f90
	$(f90comp) -c $(switch) shorthands.f90
utils.mod: utils.o utils.f90
	$(f90comp) -c $(switch) utils.f90
utils.o: parameters.o pairwise.o utils.f90
	$(f90comp) -c $(switch) utils.f90
data_structures.mod: data_structures.o data_structures.f90
	$(f90comp) -c $(switch) data_structures.f90
data_structures.o: parameters.o data_structures.f90
	$(f90comp) -c $(switch) data_structures.f90
parameters.mod: parameters.o parameters.f90
	$(f90comp) -c $(switch) parameters.f90
parameters.o: parameters.f90
	$(f90comp) -c $(switch) parameters.f90
hamilton.mod: hamilton.o hamilton.f90
	$(f90comp) -c $(switch) hamilton.f90
hamilton.o: parameters.o pairwise.o data_structures.o state.o utils.o shorthands.o hamilton.f90
	$(f90comp) -c $(switch) hamilton.f90
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
# Cleaning everything
clean:
# End of the makefile
