CXXFLAGS += -Wall -std=c++11
EXE = main
src=$(wildcard *.cpp)
OBJS = $(src:.cpp=.o)


n_x = 2
n_y = 2
x_L = -3
y_L = -2
x_R = 1.5
y_R = 2
I_max = 1000

.DEFAULT_GOAL = all

all: build

.PHONY: all build run clean distclean

build: $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

run:
	./$(EXE) $(n_x) $(n_y) $(x_L) $(y_L) $(x_R) $(y_R) $(I_max)

clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) $(EXE)
	$(RM) *~
