
CC = gcc
LD = gcc

CFLAGS = -Wall

##PDB_PATH=/global/homes/u/umansky/PUBLIC/PACT/
PDB_PATH=/global/homes/u/umansky/PUBLIC/PACT/pact/dev/lnx-2.4-o/

INCLUDE = -I$(PDB_PATH)/include
LIBS = -lm -L$(PDB_PATH)/lib -lpdb -lpml -lscore

TARGET = gato2pdb
OBJ = gato2pdb.o read_gato.o

.PHONY:all
all: $(TARGET)

$(TARGET): $(OBJ)
	$(LD) -o $(TARGET) $(OBJ) $(LIBS)

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)
