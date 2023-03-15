#ifndef ABUNDANCE_DATABASE
#define ABUNDANCE_DATABASE
#include <utility>

namespace abundance{
    int get_isotopes_num(int z);
    int get_isotope_a(int z, int i);
    float get_abundance(int z, int i);
    std::pair<int,float> get_isotope(int z, int i);
}

#endif
