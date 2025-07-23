LDFLAGS = -lm
CFLAGS = -g -Wall -Wextra -std=c17 -fPIC -fvisibility=hidden
CC = gcc

OBJDIR = Objs
DEPDIR = Deps

## all shared objects to make
OBJS = $(patsubst %.h,$(OBJDIR)/%.o,$(wildcard *.h))

# subdir hierarchy, compilation database, and all binaries
all: mkdirs compile_commands.json allBins

# all binaries
allBins: testAdjacency gbaCentrality.so

testAdjacency: $(OBJDIR)/testAdjacency.o $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

gbaCentrality.so: $(OBJS)
	$(CC) $(LDFLAGS) -shared -o $@ $^

# make subdirs if they don't exist yet
mkdirs:
	mkdir -p $(OBJDIR) $(DEPDIR)

# use bear to make a compilation database for LSP (usable by emacs and other IDEs)
# (install bear with "sudo dnf install bear" on RHEL/Alma/Fedora if you don't have it)
compile_commands.json:
	bear -- make allBins


# making each object file
$(OBJDIR)/%.o: %.c Makefile
	$(CC) $(CFLAGS) -MMD -MF $(DEPDIR)/$*.d -c -o $@ $<

# deps made thanks to -MMD above, just include them
-include $(patsubst %.c,$(DEPDIR)/%.d,$(wildcard *.c))


clean:
	rm -f $(OBJDIR)/*.o $(DEPDIR)/*.d compile_commands.json testAdjacency gbaCentrality.so

.PHONY: mkdirs all allBins clean
