CC		= gcc
OBJS		= nwalign.o matrix.o cell.o alignment.o
EXECUTABLE	= nwalign

#CFLAGS		= -g -W
CFLAGS		= -O2 -W

BINDIR		= /usr/local/bin

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

install:
	cp -p $(EXECUTABLE) $(BINDIR)

clean:
	rm -f $(OBJS) $(EXECUTABLE)

.SUFFIXES: .o .c

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@
