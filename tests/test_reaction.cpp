#include "lest.hpp"
#include <math.h>
#include "catima/catima.h"
#include "catima/reactions.h"
#include "testutils.h"
using namespace std;
using catima::approx;
using catima::reaction_rate;
using catima::nonreaction_rate;
const lest::test specification[] =
{
    CASE("reaction"){
        catima::Projectile proj{12,6,6,1000};
        auto c = catima::get_material(6);
        c.thickness(4.0);
        double r = reaction_rate(1000,c.number_density_cm2());
        std::cout<<r<<"\n";

        EXPECT(reaction_rate(0,10.0) == 0.0);
        EXPECT(reaction_rate(1000,0.0) == 0.0);
        
        EXPECT(nonreaction_rate(0,10.0) == 1.0);
        EXPECT(nonreaction_rate(1000,0.0) == 1.0);

        auto fcc = [](double x){return 1000.0;};
        r = reaction_rate(fcc, c.number_density_cm2());
        std::cout<<r<<"\n";

        catima::reaction_rate1(proj, c);
    }
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}

