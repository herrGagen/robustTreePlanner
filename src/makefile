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
INDENTED=$(subst .cpp,.ind,$(SOURCES) )

all: rtp

libtinyxml.a:
	make -C tinyxml
	ar -cvq libtinyxml.a tinyxml/tinystr.o  tinyxml/tinyxmlerror.o  tinyxml/tinyxml.o  tinyxml/tinyxmlparser.o	

rtp: ${OBJECTS} libtinyxml.a
	${CC} ${LINKFLAGS} $^  -o $@

indent: ${INDENTED}
	echo ${INDENTED}

clean:
	rm -rf *.o rtp *.bak *~

%.o: %.cpp makefile
	${CC} ${CFLAGS} $< -c

%.ind: %.cpp makefile
	emacs -batch $< -l indent.emacs -f emacs-format-function

