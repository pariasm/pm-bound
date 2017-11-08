# debug flags
CFLAGS=-O0 -fPIC -g -DDEBUG_STATS
MEXFLAGS=-g -DDEBUG -Icore

# release flags
CFLAGS=-O3 -fPIC #-fopenmp
MEXFLAGS=-Icore

SRC_NNF=modules/NNfield_l2.h modules/dImage.h
SRC_NNF_EP=modules/NNfield_l2_extpatches.h modules/dImage.h
OBJ=core/NNfield.o core/list.o

all: bin/nnfield bin/nnfield_exh bin/nnfield_ep

bin/nnfield: modules/nnfield_list_mex.c $(SRC_NNF) $(OBJ)
	mex ${MEXFLAGS} $< $(OBJ) -output $@

bin/nnfield_exh: modules/nnfield_list_mex.c $(SRC_NNF) $(OBJ)
	mex ${MEXFLAGS} -DEXHNNF $< $(OBJ) -output $@

bin/nnfield_ep: modules/nnfield_list_extpatches_mex.c $(SRC_NNF_EP) $(OBJ)
	mex ${MEXFLAGS} $< $(OBJ) -output $@

clean:
	-rm $(OBJ)

dist_clean: clean
	rm bin/nnfield.mex* bin/nnfield_exh.mex* bin/nnfield_ep.mex*
