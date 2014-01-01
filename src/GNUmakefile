PROG	 	= qprop

PACKAGE  	= $(LIBDIR)/libqprop.a

C++BASICS     	= bar.cc ylm.cc wavefunction.cc hamop.cc grid.cc\
                  fluid.cc ned.cc xml_parse.cc \
		  factorial.cc winop.cc


PROGBASICS      = 
LOCALLIBS	=
ESQLBASICS    	= 
CBASICS      	= 

include         GNUmakefile.tmpl

#	User definable options

MYLIBS          = 
COPTI           = -O8 # -pg
C++OPTI         = -O8 -Wno-deprecated # -pg
C++SPECIAL	= -funroll-all-loops -finline-functions \
                  -fexpensive-optimizations -ffast-math#-DTEST
LDFLAGS         = -Xlinker -defsym -Xlinker MAIN__=main # -pg #-gstabs+

