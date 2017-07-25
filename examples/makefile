PROGRAMS=simple example2 materials

GCC=g++ -Wall -std=c++14
INCDIR=-I$(CATIMAPATH)/include
LIBDIR=-L$(CATIMAPATH)
LIBS=-lcatima



all: $(PROGRAMS)

%: %.cpp
	$(GCC) $< $(INCDIR) $(LIBDIR) $(LIBS) -o $@

clean:
	rm $(PROGRAMS)
