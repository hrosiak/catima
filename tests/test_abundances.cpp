#include "lest.hpp"
#include <math.h>
#include "catima/abundance_database.h"
#include "testutils.h"

using namespace std;
using catima::approx;
using namespace abundance;

const lest::test specification[] =
{
    CASE("isotope abundaces"){
        EXPECT(get_isotopes_num(1) == 2);
        EXPECT(get_isotopes_num(27) == 1);
        EXPECT(get_isotopes_num(0) == 0);
        EXPECT(get_isotopes_num(111) == 0);
        
        EXPECT(get_isotope_a(0,0) == 0);
        EXPECT(get_isotope_a(0,1) == 0);
        EXPECT(get_isotope_a(1,0) == 1);
        EXPECT(get_isotope_a(1,1) == 2);
        EXPECT(get_isotope_a(1,2) == 0);
        EXPECT(get_isotope_a(14,0) == 28);
        EXPECT(get_isotope_a(14,1) == 29);
        EXPECT(get_isotope_a(14,2) == 30);
        EXPECT(get_isotope_a(14,3) == 0);
        EXPECT(get_isotope_a(68,0) == 166);
        EXPECT(get_isotope_a(68,1) == 168);
        EXPECT(get_isotope_a(120,0) == 0);
        
        EXPECT(get_abundance(1,0) == approx(0.999,0.01));
        EXPECT(get_abundance(1,2) == 0.0);
        EXPECT(get_abundance(14,0) == approx(0.922,0.01));
        EXPECT(get_abundance(14,1) == approx(0.046,0.01));
        EXPECT(get_abundance(14,2) == approx(0.0308,0.01));
        EXPECT(get_abundance(14,3) == 0.0);
        EXPECT(get_abundance(68,0) == approx(0.336,0.01));
        EXPECT(get_abundance(68,1) == approx(0.267,0.01));
        EXPECT(get_abundance(120,0) == 0);
        
        EXPECT(get_isotope(0,0).first == 0);
        EXPECT(get_isotope(0,1).first == 0);
        EXPECT(get_isotope(1,0).first == 1);
        EXPECT(get_isotope(1,1).first == 2);
        EXPECT(get_isotope(1,2).first == 0);
        EXPECT(get_isotope(14,0).first == 28);
        EXPECT(get_isotope(14,1).first == 29);
        EXPECT(get_isotope(14,2).first == 30);
        EXPECT(get_isotope(14,3).first == 0);
        EXPECT(get_isotope(68,0).first == 166);
        EXPECT(get_isotope(68,1).first == 168);
        EXPECT(get_isotope(120,0).first == 0);
        
        EXPECT(get_isotope(1,0).second == approx(0.999,0.01));
        EXPECT(get_isotope(1,2).second == 0.0);
        EXPECT(get_isotope(14,0).second == approx(0.922,0.01));
        EXPECT(get_isotope(14,1).second == approx(0.046,0.01));
        EXPECT(get_isotope(14,2).second == approx(0.0308,0.01));
        EXPECT(get_isotope(14,3).second == 0.0);
        EXPECT(get_isotope(68,0).second == approx(0.336,0.01));
        EXPECT(get_isotope(68,1).second == approx(0.267,0.01));
        EXPECT(get_isotope(120,0).second == 0);
        
        
    }
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}


