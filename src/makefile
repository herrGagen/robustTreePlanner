CC=g++
CFLAGS=-O3 -pg -DOUTPUT_WEATHER_IN_XML -DDP -DDP_USES_ONE_NEIGHBOR_LAYER -DSUPPRESS_OUTPUT
LINKFLAGS=-pg

SOURCES=DemandProfile.cpp \
	NodeAndEdge.cpp \
	RoutingDAG.cpp \
	WeatherData.cpp \
	main.cpp \
	Quadrant.cpp \
	UserInterface.cpp \
	InputFileReader.cpp \
	TreeVerifier.cpp \
	DynProgTreeGenerator.cpp

OBJECTS=$(subst .cpp,.o,$(SOURCES) )

all: rtp

rtp: ${OBJECTS}
	${CC} ${LINKFLAGS} $^  -o $@


clean:
	rm -rf *.o rtp *.bak *~

%.o: %.cpp makefile
	emacs -batch $< -l indent.emacs -f emacs-format-function
	${CC} ${CFLAGS} $< -c