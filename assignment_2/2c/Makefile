CXXFLAGS += -Wall -std=c++11 -g
CFLAGS += -Wall -g
EXE = main

CPP_SRC = $(wildcard *.cpp)
C_SRC = $(wildcard *.c)
CPP_OBJS = $(CPP_SRC:.cpp=.o)
C_OBJS = $(C_SRC:.c=.o)
OBJS = $(CPP_OBJS) $(C_OBJS)

.DEFAULT_GOAL = all
n_x =  800
n_y = 600
x_L = -2.0
y_L = -1.2
x_R = 1.0
y_R = 1.2
I_max = 1000

all: build

build: $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

run:
	./$(EXE) $(n_x) $(n_y) $(x_L) $(y_L) $(x_R) $(y_R) $(I_max)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	rm -f $(OBJ) $(EXE)

distclean: clean
	rm -f *~

print-%: ; @echo $* = $($*)

