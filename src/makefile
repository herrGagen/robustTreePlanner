CC=g++
CFLAGS=-O3 -pg
LINKFLAGS=-pg

SOURCES=DemandProfile.cpp \
	NodeAndEdge.cpp \
	RoutingDAG.cpp \
	WeatherData.cpp \
	main.cpp \
	Quadrant.cpp \
	UserInterface.cpp \
	InputFileReader.cpp \
	TreeVerifier.cpp	

OBJECTS=$(subst .cpp,.o,$(SOURCES) )

all: rtp

rtp: ${OBJECTS}
	${CC} ${LINKFLAGS} $^  -o $@	

clean:
	rm -rf *.o rtp *.bak *~

%.o: %.cpp makefile
	${CC} ${CFLAGS} $< -c