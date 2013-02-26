CC=g++
CFLAGS=-O3

SOURCES=DemandProfile.cpp \
	NodeAndEdge.cpp \
	RoutingDAG.cpp \
	WeatherData.cpp \
	main.cpp \
	Quadrant.cpp \
	UserInterface.cpp

OBJECTS=$(subst .cpp,.o,$(SOURCES) )

all: rtp

rtp: ${OBJECTS}
	${CC} $^  -o $@	

clean:
	rm -rf *.o rtp *.bak *~

%.o: %.cpp makefile
	${CC} ${CFLAGS} $< -c