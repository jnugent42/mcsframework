SRCS=MCSUnfolding.cpp Collection.cpp MCSAnalysis.cpp
OBJS=$(SRCS:.cpp=.o) 
CC = g++
DEBUG  = -g
CFLAGS = -Wall $(DEBUG)
IFLAGS=-I${MAUS_ROOT_DIR}/ \
	 -I${MAUS_ROOT_DIR}/src/common_cpp \
	 -I${MAUS_ROOT_DIR}/src/legacy \
	 -I${MAUS_THIRD_PARTY}/third_party/install/include \
	 -I${MAUS_THIRD_PARTY}/third_party/install/include/libxml2 \
	 -I${MAUS_THIRD_PARTY}/third_party/install/include/python2.7 \
	 -I${ROOTSYS}/include \
	 -I/afs/phas.gla.ac.uk/user/j/jnugent/workarea/Unfolding/RooUnfold-1.1.1/src
LFLAGS=-L${MAUS_THIRD_PARTY}/third_party/build/root/lib -lPyROOT \
         -ljson \
         -L${MAUS_THIRD_PARTY}/third_party/build/root/lib -lSpectrum \
         -L${MAUS_THIRD_PARTY}/third_party/build/gcc-4.9.3/obj/x86_64-unknown-linux-gnu/libstdc++-v3/src/.libs -lstdc++ \
         -L${MAUS_THIRD_PARTY}/third_party/install/lib \
         -L${MAUS_THIRD_PARTY}/build \
         -L/afs/phas.gla.ac.uk/user/j/jnugent/workarea/Unfolding/RooUnfold-1.1.1 \
	 -lRooUnfold \
	 `root-config --ldflags` \
	 `${ROOTSYS}/bin/root-config --glibs` \
	 -lMausCpp \
	 -L${MAUS_ROOT_DIR}/third_party/lib -lxml2\
         -L/usr/lib64 \
         -L/usr/lib64/root \
	 -Wl,--no-as-needed
EXECUTABLE=MCSUnfolding reduce_tof2_trigger reduce_tof1_trigger

all : $(SRCS) $(EXECUTABLE)

.cpp.o:
	$(CC) -c $(CFLAGS) -o $@ $< $(IFLAGS) $(LFLAGS)


MCSUnfolding: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(IFLAGS) $(OBJS) $(LFLAGS)

reduce_tof2_trigger: reduce_tof2_trigger.cc
	$(CC) $(CFLAGS) reduce_tof2_trigger.cc -o $@ $(IFLAGS) $(LFLAGS)

reduce_tof1_trigger: reduce_tof1_trigger.cc
	$(CC) $(CFLAGS) reduce_tof1_trigger.cc -o $@ $(IFLAGS) $(LFLAGS)

clean:
	\rm *.o $(EXECUTABLE)


