#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include <math.h>
#include "catima/abundance_database.h"
#include "testutils.h"

using namespace std;
using namespace abundance;

    TEST_CASE("isotope abundaces"){
        CHECK(get_isotopes_num(1) == 2);
        CHECK(get_isotopes_num(27) == 1);
        CHECK(get_isotopes_num(0) == 0);
        CHECK(get_isotopes_num(111) == 0);
        
        CHECK(get_isotope_a(0,0) == 0);
        CHECK(get_isotope_a(0,1) == 0);
        CHECK(get_isotope_a(1,0) == 1);
        CHECK(get_isotope_a(1,1) == 2);
        CHECK(get_isotope_a(1,2) == 0);
        CHECK(get_isotope_a(14,0) == 28);
        CHECK(get_isotope_a(14,1) == 29);
        CHECK(get_isotope_a(14,2) == 30);
        CHECK(get_isotope_a(14,3) == 0);
        CHECK(get_isotope_a(68,0) == 166);
        CHECK(get_isotope_a(68,1) == 168);
        CHECK(get_isotope_a(120,0) == 0);
        
        CHECK(get_abundance(1,0) == approx(0.999,0.01));
        CHECK(get_abundance(1,2) == 0.0);
        CHECK(get_abundance(14,0) == approx(0.922,0.01));
        CHECK(get_abundance(14,1) == approx(0.046,0.01));
        CHECK(get_abundance(14,2) == approx(0.0308,0.01));
        CHECK(get_abundance(14,3) == 0.0);
        CHECK(get_abundance(68,0) == approx(0.336,0.01));
        CHECK(get_abundance(68,1) == approx(0.267,0.01));
        CHECK(get_abundance(120,0) == 0);
        
        CHECK(get_isotope(0,0).first == 0);
        CHECK(get_isotope(0,1).first == 0);
        CHECK(get_isotope(1,0).first == 1);
        CHECK(get_isotope(1,1).first == 2);
        CHECK(get_isotope(1,2).first == 0);
        CHECK(get_isotope(14,0).first == 28);
        CHECK(get_isotope(14,1).first == 29);
        CHECK(get_isotope(14,2).first == 30);
        CHECK(get_isotope(14,3).first == 0);
        CHECK(get_isotope(68,0).first == 166);
        CHECK(get_isotope(68,1).first == 168);
        CHECK(get_isotope(120,0).first == 0);
        
        CHECK(get_isotope(1,0).second == approx(0.999,0.01));
        CHECK(get_isotope(1,2).second == 0.0);
        CHECK(get_isotope(14,0).second == approx(0.922,0.01));
        CHECK(get_isotope(14,1).second == approx(0.046,0.01));
        CHECK(get_isotope(14,2).second == approx(0.0308,0.01));
        CHECK(get_isotope(14,3).second == 0.0);
        CHECK(get_isotope(68,0).second == approx(0.336,0.01));
        CHECK(get_isotope(68,1).second == approx(0.267,0.01));
        CHECK(get_isotope(120,0).second == 0);
        
        
    }
