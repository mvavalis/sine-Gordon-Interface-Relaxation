CC = g++

DEALIIdir = /home/manolis/pdd/deal.II
DEALIIinclude = $(patsubst %, ${DEALIIdir}/%, deal.II/include base/include lac/include contrib/boost/include)
DEALIIlibdir = ${DEALIIdir}/lib
DEALIIlibfiles = libbase  libdeal_II_1d  libdeal_II_2d  libdeal_II_3d  liblac

## shared
default: DEALIIlib = $(patsubst %, ${DEALIIlibdir}/%.so, ${DEALIIlibfiles}) 
## static
#default: DEALIIlib = $(patsubst %, ${DEALIIlibdir}/%.a, ${DEALIIlibfiles}) 

default: CCFLAGS = -Wall
default: main.o


main.o: main.cpp shared_lib_export
	${CC} ${CCFLAGS} -o main.out main.cpp $(patsubst %, -I%, ${DEALIIinclude}) ${DEALIIlib}
	#${CC} ${CCFLAGS} -o batch_test.out batch_test.cpp

shared_lib_export:
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DEALIIlibdir}

clean:
	rm -rf *.o *~ main.out
