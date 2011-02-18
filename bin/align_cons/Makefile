CC		= gcc
OBJS		= align_cons.o
EXECUTABLE	= align_cons

#CFLAGS		= -g -Wall
CFLAGS		= -O2 -Wall

BINDIR		= /usr/local/bin

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -lm -o $@

install:
	cp -p $(EXECUTABLE) $(BINDIR)

clean:
	rm -f $(OBJS) $(EXECUTABLE)

.SUFFIXES: .o .c

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@
