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
        catima::Projectile proj{12,6,6,870};
        auto c = catima::get_material(6);
        c.thickness(2.0);
        double r;

        r= catima::nonreaction_rate(proj, c);
        EXPECT(r == approx(0.92,0.02));

        catima::Layers l;
        l.add(c);
        l.add(c);
        l.add(c);

        auto res = catima::calculate(proj,l);
        EXPECT(res.total_result.sp == approx(r*r*r,0.02));

        c.thickness(6.0);
        r= catima::nonreaction_rate(proj, c);
        EXPECT(res.total_result.sp == approx(r,0.001));
        
        c.thickness(0.0);
        r= catima::nonreaction_rate(proj, c);
        EXPECT(r == 1.0);
        proj.T = 0.0;
        c.thickness(6);
        r= catima::nonreaction_rate(proj, c);
        EXPECT(r == -1.0);

    }
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}

