CXXFLAGS =	 -g -fmessage-length=0 -I ../ -Wno-deprecated -Wno-ignored-attributes -Wno-attributes -c -std=c++11

OBJS = mprascore.o \
	../math/rand.o \
	../math/statfun.o \
	../math/variance.o \
	../system/consoleapp.o \
	../system/environment.o \
	../system/filehelpers.o \
	../system/filesystem.o \
	../system/linereader.o \
	../system/stderr.o \
	../system/stdout.o \
	../types/matrix.o \
	../types/string.o \
	../types/stringfun.o \
	../types/tablefun.o \
	../types/tableindex.o \
	../types/vector.o

LIBS =

TARGET =	mprascore

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
