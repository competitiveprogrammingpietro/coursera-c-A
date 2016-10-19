CPP=g++

all : first second third

first: homework-one/basics.cpp
	$(CPP) $^ -o $@
second: homework-two/graph.cpp
	$(CPP) $^ -o $@
third: homework-three/graph.cpp
	$(CPP) $^ -o $@


clean:
	@rm first second third 2>/dev/null

.PHONY: clean
