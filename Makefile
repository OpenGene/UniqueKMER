DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := uniquekmer

BIN_TARGET := ${TARGET}

CXX ?= g++
CXXFLAGS := -std=c++11 -Xpreprocessor -fopenmp -g -O3 -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) ${CXXFLAGS}

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
LIBS := -lz -lpthread -lomp
else
LIBS := -lz -lpthread
endif

LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)


${BIN_TARGET}:${OBJ}
ifeq ($(UNAME_S),Darwin)
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS)
else
	$(CXX) $(OBJ) -o $@ $(LD_FLAGS) -fopenmp
endif

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp make_obj_dir
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY:clean
clean:
	rm obj/*.o
	rm $(TARGET)

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."
