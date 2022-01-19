# bamcmp Makefile

ifndef HTSLIBDIR
   HTSLIBDIR=/usr/local
endif


CC=gcc
CFLAGS=-g -Wall-fexceptions
CPP=g++
CPPFLAGS=-g -Wall -fexceptions
LDFLAG=-g
LDLIBS=hts

SRCDIR=src/
INCDIR=include/
BUILDDIR=build/
SRCS=SamReader.cpp BamRecVector.cpp HTSFileWrapper.cpp util.cpp bamcmp.cpp
OBJS=$(BUILDDIR)SamReader.o $(BUILDDIR)BamRecVector.o \
	$(BUILDDIR)HTSFileWrapper.o $(BUILDDIR)util.o $(BUILDDIR)bamcmp.o

bamcmp: ${OBJS} $(BUILDDIR)
	$(CPP) $(LDFLAGS) -o $(BUILDDIR)/bamcmp $(OBJS) -L $(HTSLIBDIR)/lib -l $(LDLIBS) -Wl,-rpath,/usr/local/lib

$(BUILDDIR)%.o: $(SRCDIR)%.cpp $(BUILDDIR)
	$(CPP) $(CPPFLAGS) -I $(INCDIR) -I $(HTSLIBDIR)/include -o $@ -c $< 

$(BUILDDIR):
	mkdir $(BUILDDIR)

clean:
	rm -f $(OBJS) $(BUILDIR)/bamcmp
