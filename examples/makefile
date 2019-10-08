PROGRAMS=simple dedx materials ls_coefficients

GCC=g++ -Wall -std=c++14
INCDIR=-I$(CATIMAPATH)/include
LIBDIR=-L$(CATIMAPATH)/lib
LIBS=-lcatima



all: $(PROGRAMS)

%: %.cpp
	$(GCC) $< $(INCDIR) $(LIBDIR) $(LIBS) -o $@

clean:
	rm $(PROGRAMS)
