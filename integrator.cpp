#include "integrator.h"

#ifdef GSL_INTEGRATION
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#endif

namespace catima{
    integrator_type integrator;
#ifdef GSL_INTEGRATION
    double funcwrapper3(double x, void *_c){
    std::function<double(double)> *f = (std::function<double(double)> *)_c;
    return (*f)(x);
    }

    IntegratorGSL::IntegratorGSL(bool adapt):adaptive(adapt){
        gsl_set_error_handler_off();
        if(adaptive){
            w=gsl_integration_workspace_alloc(100);
        }
    };
     
    IntegratorGSL::~IntegratorGSL(){
        if(adaptive){
            gsl_integration_workspace_free(w);
        }
    };
    
    double IntegratorGSL::integrate(std::function<double(double)> f, double _min, double  _max, double prec){
        gsl_function F;
        F.function = funcwrapper3;
        F.params = &f;

        min = _min;
        max = _max;
        size_t num;
        if(adaptive){
            #ifdef USE_THREADS
            std::lock_guard<std::mutex> lock(integration_mutex);
            #endif
            gsl_integration_qag(&F,_min,_max,1e-6,prec,100,6,w,&result,&error);
        }
        else{
            gsl_integration_qng(&F,_min,_max,1e-6,prec,&result,&error,&num);
            }
        return result;
    };
#endif
}
