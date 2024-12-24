TARGET        = mpi-kmeans


MPICC         = mpicc

SRCDIR        = src
OBJDIR        = obj
BINDIR        = bin

CFLAGS        = -Wall -std=c90 -O2 -I./$(SRCDIR)

LDFLAGS       = -lm

C_SOURCES     = $(wildcard $(SRCDIR)/*.c)

C_OBJECTS     = $(C_SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

.PHONY: all
all: remove $(BINDIR)/$(TARGET) clean

$(BINDIR)/$(TARGET): $(C_OBJECTS) $(BINDIR)
	$(MPICC) $(C_OBJECTS) $(LDFLAGS) -o $@

windows: remove build_win clean

build_win: $(C_OBJECTS) $(BINDIR)
	$(MPICC) $(C_OBJECTS)  $(LDFLAGS) -o $(BINDIR)/$(TARGET).exe

$(C_OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c $(OBJDIR)
	$(MPICC) $(CFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $@

$(BINDIR):
	mkdir -p $@

clean:
	$(RM) -rf $(C_OBJECTS) $(OBJDIR)

remove:
	$(RM) -rf $(C_OBJECTS) $(OBJDIR) $(BINDIR)