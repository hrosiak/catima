/*
 *  Copyright(C) 2017, Andrej Prochazka
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef GKQ_INTEGRATOR_H
#define GKQ_INTEGRATOR_H
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace integrators{

template<int order>
struct GK_data{
};

/**
 * @brief adaptive integrator 
 * following orders are supported: 7,15,21,31,61
 * @tparam order Kronrod order
 */
template<int order>
class GaussKronrodIntegration{
public:
    template<typename F>
    static double integrate(F& f, double a, double b, double eps = 1e-3, double reps=1e-6, int N = 1);

    template<typename F>
    static double integrate_intervals(F& f, const std::vector<std::pair<double,double>>& intervals, double eps = 1e-3, double reps=1e-6);

    template<typename F>
    static std::pair<double,double> integrate_nonadaptive(F& f, double a, double b);

    template<typename F>
    static double integrate_adaptive(F& f, double a, double b, double eps = 1e-3, double reps=1e-6, int level=49);

    static double w(int i) {return GK_data<order>::w()[i];}
    static double wg(int i) {return GK_data<order>::wg()[i];}
    static double x(int i) {return GK_data<order>::x()[i];}
    static int n() {return order;}
    static std::array<double,order> get_points(double a = -1.0, double b = 1.0);

private:
    static const size_t ngauss = (order-1)/2;
    static const size_t xsize = order/2 + 1;
    static const bool oddgauss = ((order-1)/2)&1;
};

template<int order>
template<typename F>
std::pair<double,double> GaussKronrodIntegration<order>::integrate_nonadaptive(F& f, double a, double b){
    double res;
    double gres = 0.0;
    double val;
    double p = 0.5*(b-a);
    double q = 0.5*(b+a);
    unsigned int gi = 1; //by default assume odd-order gauss
    unsigned int ki = 2;

    // 1st kronrod point
    val = f(p*x(0) + q);
    res = w(0) * val;
    if(oddgauss){ // 1st gaus point if odd order gauss
        gres = wg(0) * val;
        gi = 2;
        ki = 1;
        }

    for(unsigned int i=gi;i<xsize;i+=2){
        val = f(p*x(i) + q);
        res += w(i) * val;
        gres += wg(i/2)*val;

        val = f(-p*x(i) + q);
        res += w(i) * val;
        gres += wg(i/2)*val;
    }

    for(unsigned int i=ki;i<xsize;i+=2){
        res += w(i) * (f(p*x(i) + q) + f(-p*x(i) + q));
    }

    double err = std::max( std::abs(gres-res), std::numeric_limits<double>::epsilon());

    return std::make_pair(p*res, p*err);
}

template<int order>
template<typename F>
double GaussKronrodIntegration<order>::integrate_adaptive(F& f, double a, double b, double eps, double rprec, int level){
    double result = 0.0;
    double err = 0.0;
    const double numlimit = 10*std::numeric_limits<double>::epsilon();

    auto r = integrate_nonadaptive(f, a, b);
    result = r.first;
    err = r.second;
    //printf("level = %d, I=%lf, e=%lf, %lf - %lf\n",level,result,err,a,b);
    if( (std::abs(result) < numlimit) || ((b-a)<numlimit)){
        return result;
    }

    double aeps = std::max(rprec*std::abs(result), eps);
    if((aeps<numlimit) || (std::abs(result)<aeps)){
	    return result;
    }

    if( level && (err > aeps)){
	    double mid = 0.5*(a+b);
        result  = integrate_adaptive(f, a, mid, 0.707*aeps, 0, level-1);
        result += integrate_adaptive(f, mid, b, 0.707*aeps, 0, level-1);
        }

    return result;
    }

template<int order>
template<typename F>
double GaussKronrodIntegration<order>::integrate(F& f, double a, double b, double eps, double reps, int N){
    double step = (b-a)/N;
    double result = 0;
    for(int i=0; i<N;i++){
        double m =a+(i*step);
        result+=integrate_adaptive(f,m,m+step, eps/N, reps);
    }

    return result;
}

template<int order>
template<typename F>
double GaussKronrodIntegration<order>::integrate_intervals(F& f, const std::vector<std::pair<double,double>>& intervals, double eps, double reps){
    double result = 0;
    for(auto& i:intervals){
        result+=integrate_adaptive(f,i.first,i.second, eps/intervals.size(), reps);
    }

    return result;
}

template<int order>
std::array<double,order> GaussKronrodIntegration<order>::get_points(double a,  double b){
    std::array<double,order> points;
    double p = 0.5*(b-a);
    double q = 0.5*(b+a);

    int num = (order/2);
    for(int i=0;i< num;i++){
        points[num-i-1] = -p*x(i) + q;
        points[num+i] = p*x(i) + q;
    }
    return points;
}


/// weights and abscissas
// order = 7
template<>
struct GK_data<7>{
    static std::array<double,4> const & x(){
        static const std::array<double,4> _x =
            {0.0, 0.4342437493468025580021, 0.774596669241483377036,0.9604912687080202834235};
        return _x;
    }
    static std::array<double,4> const & w(){
        static const std::array<double,4> _w =
            {0.4509165386584741423451, 0.4013974147759622229051, 0.2684880898683334407286, 0.1046562260264672651938};
        return _w;
        }
    static std::array<double,4> const & wg(){
        static const std::array<double,4> _wg =
            {0.8888888888888888888889, 0.555555555555555555556};
        return _wg;
        }
};

//order = 15
template<>
struct GK_data<15>{
    static std::array<double,8> const & x(){
        static const std::array<double,8> _x =
            {0.0,
            0.207784955007898467600689403773244913,
            0.405845151377397166906606412076961463,
            0.586087235467691130294144838258729598,
            0.741531185599394439863864773280788407,
            0.864864423359769072789712788640926201,
            0.949107912342758524526189684047851262,
            0.991455371120812639206854697526328517};
            return _x;
    }
    static std::array<double,8> const & w(){
        static const std::array<double,8> _w =
            {
         0.209482141084727828012999174891714264,
         0.204432940075298892414161999234649085,
         0.190350578064785409913256402421013683,
         0.169004726639267902826583426598550284,
         0.140653259715525918745189590510237920,
         0.104790010322250183839876322541518017,
         0.0630920926299785532907006631892042867,
         0.0229353220105292249637320080589695920};
        return _w;
        }
    static std::array<double,4> const & wg(){
        static const std::array<double,4> _wg =
        { 0.417959183673469387755102040816326531,
          0.381830050505118944950369775488975134,
          0.279705391489276667901467771423779582,
          0.129484966168869693270611432679082018};
        return _wg;
        }
};

//order = 21
template<>
struct GK_data<21>{
    static std::array<double,11> const & x(){
        static const std::array<double,11> _x = {
            0.00000000000000000e+00,
            1.48874338981631211e-01,
            2.94392862701460198e-01,
            4.33395394129247191e-01,
            5.62757134668604683e-01,
            6.79409568299024406e-01,
            7.80817726586416897e-01,
            8.65063366688984511e-01,
            9.30157491355708226e-01,
            9.73906528517171720e-01,
            9.95657163025808081e-01};
            return _x;
    }
    static std::array<double,11> const & w(){
        static const std::array<double,11> _w =
            {
            1.49445554002916906e-01,
            1.47739104901338491e-01,
            1.42775938577060081e-01,
            1.34709217311473326e-01,
            1.23491976262065851e-01,
            1.09387158802297642e-01,
            9.31254545836976055e-02,
            7.50396748109199528e-02,
            5.47558965743519960e-02,
            3.25581623079647275e-02,
            1.16946388673718743e-02,
            };
        return _w;
        }
    static std::array<double,5> const & wg(){
        static const std::array<double,5> _wg =
        {
            2.95524224714752870e-01,
            2.69266719309996355e-01,
            2.19086362515982044e-01,
            1.49451349150580593e-01,
            6.66713443086881376e-02,};
        return _wg;
        }
};



//order = 31
template<>
struct GK_data<31>{
    static std::array<double,16> const & x(){
        static const std::array<double,16> _x = {
            0.00000000000000000e+00,
            1.01142066918717499e-01,
            2.01194093997434522e-01,
            2.99180007153168812e-01,
            3.94151347077563370e-01,
            4.85081863640239681e-01,
            5.70972172608538848e-01,
            6.50996741297416971e-01,
            7.24417731360170047e-01,
            7.90418501442465933e-01,
            8.48206583410427216e-01,
            8.97264532344081901e-01,
            9.37273392400705904e-01,
            9.67739075679139134e-01,
            9.87992518020485428e-01,
            9.98002298693397060e-01    
            };
            return _x;
    }
    static std::array<double,16> const & w(){
        static const std::array<double,16> _w =
            {
            1.01330007014791549e-01,
            1.00769845523875595e-01,
            9.91735987217919593e-02,
            9.66427269836236785e-02,
            9.31265981708253212e-02,
            8.85644430562117706e-02,
            8.30805028231330210e-02,
            7.68496807577203789e-02,
            6.98541213187282587e-02,
            6.20095678006706403e-02,
            5.34815246909280873e-02,
            4.45897513247648766e-02,
            3.53463607913758462e-02,
            2.54608473267153202e-02,
            1.50079473293161225e-02,
            5.37747987292334899e-03
            };
        return _w;
        }
    static std::array<double,8> const & wg(){
        static const std::array<double,8> _wg =
        {
        2.02578241925561273e-01,
        1.98431485327111576e-01,
        1.86161000015562211e-01,
        1.66269205816993934e-01,
        1.39570677926154314e-01,
        1.07159220467171935e-01,
        7.03660474881081247e-02,
        3.07532419961172684e-02
        };
        return _wg;
        }
};

//order = 61
template<>
struct GK_data<61>{
    static std::array<double,31> const & x(){
        static const std::array<double,31> _x = {
         0.00000000000000000e+00,
         5.14718425553176958e-02,
         1.02806937966737030e-01,
         1.53869913608583547e-01,
         2.04525116682309891e-01,
         2.54636926167889846e-01,
         3.04073202273625077e-01,
         3.52704725530878113e-01,
         4.00401254830394393e-01,
         4.47033769538089177e-01,
         4.92480467861778575e-01,
         5.36624148142019899e-01,
         5.79345235826361692e-01,
         6.20526182989242861e-01,
         6.60061064126626961e-01,
         6.97850494793315797e-01,
         7.33790062453226805e-01,
         7.67777432104826195e-01,
         7.99727835821839083e-01,
         8.29565762382768397e-01,
         8.57205233546061099e-01,
         8.82560535792052682e-01,
         9.05573307699907799e-01,
         9.26200047429274326e-01,
         9.44374444748559979e-01,
         9.60021864968307512e-01,
         9.73116322501126268e-01,
         9.83668123279747210e-01,
         9.91630996870404595e-01,
         9.96893484074649540e-01,
         9.99484410050490638e-01
            };
            return _x;
    }
    static std::array<double,31> const & w(){
        static const std::array<double,31> _w =
            {
         5.14947294294515676e-02,
         5.14261285374590259e-02,
         5.12215478492587722e-02,
         5.08817958987496065e-02,
         5.04059214027823468e-02,
         4.97956834270742064e-02,
         4.90554345550297789e-02,
         4.81858617570871291e-02,
         4.71855465692991539e-02,
         4.60592382710069881e-02,
         4.48148001331626632e-02,
         4.34525397013560693e-02,
         4.19698102151642461e-02,
         4.03745389515359591e-02,
         3.86789456247275930e-02,
         3.68823646518212292e-02,
         3.49793380280600241e-02,
         3.29814470574837260e-02,
         3.09072575623877625e-02,
         2.87540487650412928e-02,
         2.65099548823331016e-02,
         2.41911620780806014e-02,
         2.18280358216091923e-02,
         1.94141411939423812e-02,
         1.69208891890532726e-02,
         1.43697295070458048e-02,
         1.18230152534963417e-02,
         9.27327965951776343e-03,
         6.63070391593129217e-03,
         3.89046112709988405e-03,
         1.38901369867700762e-03
            };
        return _w;
        }
    static std::array<double,15> const & wg(){
        static const std::array<double,15> _wg =
        {
         1.02852652893558840e-01,
         1.01762389748405505e-01,
         9.95934205867952671e-02,
         9.63687371746442596e-02,
         9.21225222377861287e-02,
         8.68997872010829798e-02,
         8.07558952294202154e-02,
         7.37559747377052063e-02,
         6.59742298821804951e-02,
         5.74931562176190665e-02,
         4.84026728305940529e-02,
         3.87991925696270496e-02,
         2.87847078833233693e-02,
         1.84664683110909591e-02,
         7.96819249616660562e-03
        };
        return _wg;
        }
};


} // end of namespace
#endif
