CXX = g++
CXXFLAGS = -I/home/syoukera/miniconda3/envs/ct-build/include
LDFLAGS = -L/home/syoukera/miniconda3/envs/ct-build/lib
LIBS = -lcantera_shared -Wl,-rpath,/home/syoukera/miniconda3/envs/ct-build/lib

TARGET = main
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
