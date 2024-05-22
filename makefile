CXX = c++
INCLUDES = -I/users/ziyzhang/topology-research/nexullance/IT_boost/boost
CXXFLAGS = -std=c++14 -O0 -g -fexceptions
# CXXFLAGS = -std=c++14 -O3 -ffast-math -Wall -pedantic
LDFLAGS = -L/users/ziyzhang/topology-research/nexullance/IT_boost/boost/stage/lib
LDLIBS = -lboost_graph

TARGET = main
SRCS := $(wildcard *.cpp)  # List of all .cpp files in the current directory
OBJS := $(SRCS:.cpp=.o)     # Corresponding .o files

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(TARGET)

# Compile each .cpp file into a corresponding .o file
%.o: %.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean
