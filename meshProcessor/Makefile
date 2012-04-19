SRC=mesh.cpp
OBJ=$(SRC:.cpp=.o)
TARGET=libmesh3d.a
CXXFLAGS= -m32 -O2 -Wall

.PHONY: all clean

all: $(TARGET)

clean:: 
	rm -f $(TARGET)

$(TARGET) : $(OBJ)
	ar rcs $@ $^

clean::
	rm -f $(OBJ) 
