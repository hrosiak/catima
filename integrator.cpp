#include "integrator.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"

namespace catima{
    
    IntegratorGSL integratorGSL;
    
    double funcwrapper3(double x, void *_c){
    std::function<double(double)> *f = (std::function<double(double)> *)_c;
    return (*f)(x);
    }

    IntegratorGSL::IntegratorGSL(){
        gsl_set_error_handler_off();
        w=gsl_integration_workspace_alloc(500);
    };
    
    IntegratorGSL::~IntegratorGSL(){
        gsl_integration_workspace_free(w);
    };
    
    double IntegratorGSL::Integrate(std::function<double(double)> f, double _min, double  _max, double prec, bool adaptive){
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
            gsl_integration_qag(&F,_min,_max,prec,prec,500,6,w,&result,&error);
        }
        else
            gsl_integration_qng(&F,_min,_max,0,prec,&result,&error,&num);
        return result;
    };
}
