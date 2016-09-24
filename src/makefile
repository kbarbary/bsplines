
OS ?= $(shell sh -c 'uname -s | tr "[A-Z]" "[a-z]"')

NAME = bs
SOMAJOR = 0
SOMINOR = 1
SOBUGFIX = 0

ifeq ($(OS), darwin)
SONAME = lib$(NAME).dylib
SONAME_MAJOR = lib$(NAME).$(SOMAJOR).dylib
SONAME_MAJORMINOR = lib$(NAME).$(SOMAJOR).$(SOMINOR).dylib
SONAME_FULL = lib$(NAME).$(SOMAJOR).$(SOMINOR).$(SOBUGFIX).dylib
SONAME_FLAG = -install_name
LDPATHENV = DYLD_LIBRARY_PATH
else
ifeq ($(OS), linux)
SONAME = lib$(NAME).so
SONAME_MAJOR = lib$(NAME).so.$(SOMAJOR)
SONAME_MAJORMINOR = lib$(NAME).so.$(SOMAJOR).$(SOMINOR)
SONAME_FULL = lib$(NAME).so.$(SOMAJOR).$(SOMINOR).$(SOBUGFIX)
SONAME_FLAG = -soname
LDPATHENV = LD_LIBRARY_PATH
else
$(error OS not yet supported)
endif
endif

INSTALL ?= install
DESTDIR =
PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include

CC?=gcc
AR?=ar
CPPFLAGS ?=
LDFLAGS ?=

CPPFLAGS += -I.
CFLAGS += -Wall -Wextra -O3  -std=c99 # -Werror
CFLAGS_LIB = $(CFLAGS) -fPIC
LDFLAGS_LIB = $(LDFLAGS) -shared -Wl,$(SONAME_FLAG),$(SONAME_MAJOR)

default: all

bs.o: bs.c bs.h
	$(CC) $(CPPFLAGS) $(CFLAGS_LIB) -c bs.c -o $@

$(SONAME): bs.o
	$(CC) $(LDFLAGS_LIB) $^ -o $(SONAME_FULL)
	ln -sf $(SONAME_FULL) $(SONAME_MAJORMINOR)
	ln -sf $(SONAME_FULL) $(SONAME_MAJOR)
	ln -sf $(SONAME_FULL) $(SONAME)

lib$(NAME).a: $(OBJS)
	$(AR) rcs lib$(NAME).a $^

install: library
	$(INSTALL) -D $(NAME).h $(DESTDIR)$(INCLUDEDIR)/$(NAME).h
	$(INSTALL) -D $(SONAME_FULL) $(DESTDIR)$(LIBDIR)/$(SONAME_FULL)
	$(INSTALL) -D $(SONAME_MAJORMINOR) $(DESTDIR)$(LIBDIR)/$(SONAME_MAJORMINOR)
	$(INSTALL) -D $(SONAME_MAJOR) $(DESTDIR)$(LIBDIR)/$(SONAME_MAJOR)
	$(INSTALL) -D $(SONAME) $(DESTDIR)$(LIBDIR)/$(SONAME)
	$(INSTALL) -D lib$(NAME).a $(DESTDIR)$(LIBDIR)/lib$(NAME).a

uninstall:
	rm $(DESTDIR)$(INCLUDEDIR)/$(NAME).h
	rm $(DESTDIR)$(LIBDIR)/$(SONAME_FULL)
	rm $(DESTDIR)$(LIBDIR)/$(SONAME_MAJORMINOR)
	rm $(DESTDIR)$(LIBDIR)/$(SONAME_MAJOR)
	rm $(DESTDIR)$(LIBDIR)/$(SONAME)
	rm $(DESTDIR)$(LIBDIR)/lib$(NAME).a

test: runtests
	$(LDPATHENV)=. runtests

runtests: runtests.c $(NAME).c $(NAME).h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -L. runtests.c -o runtests

clean:
	rm -f *.o lib$(NAME).* runtests

all: $(SONAME) lib$(NAME).a

.PHONY: all test clean library install uninstall
