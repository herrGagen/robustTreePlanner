all: build

build: main.o demandprofile.o nodeandedge.o quadrant.o routingdag.o userinterface.o weatherdata.o
	g++ main.o demandprofile.o nodeandedge.o quadrant.o routingdag.o userinterface.o weatherdata.o -o rtp

main.o: main.cpp
	g++ -c main.cpp

demandprofile.o: DemandProfile.cpp
	g++ -c DemandProfile.cpp

nodeandedge.o: NodeAndEdge.cpp
	g++ -c NodeAndEdge.cpp

quadrant.o: Quadrant.cpp
	g++ -c Quadrant.cpp

routingdag.o: RoutingDag.cpp
	g++ -c RoutingDag.cpp

userinterface.o: UserInterface.cpp
	g++ -c UserInterface.cpp

weatherdata.o: WeatherData.cpp
	g++ -c WeatherData.cpp

clean:
	rm -rf *o rtp