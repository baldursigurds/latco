DESTDIR=/usr
libdir=/lib
includedir=/include

all:
	gcc -fPIC -Wall -g -c latco.c
	gcc -g -shared -Wl,-soname,liblatco.so.0 \
	-o liblatco.so.0.0 latco.o -lc

install:
	mkdir -p $(DESTDIR)$(libdir)
	install liblatco.so.0.0 $(DESTDIR)$(libdir)/liblatco.so.0.0
	ln -s $(DESTDIR)$(libdir)/liblatco.so.0.0 $(DESTDIR)$(libdir)/liblatco.so.0
	mkdir -p $(DESTDIR)$(includedir)
	install latco.h $(DESTDIR)$(includedir)/latco.h
