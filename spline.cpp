/*
 * This is modification of Tino Kluge tk spline
 * https://github.com/ttk592/spline/
 *
 *  the modification is in LU caclulation, 
 * optimized for tridiagonal matrices
 */
#include "spline.h"


namespace catima{

band_matrix::band_matrix(int dim)
{
    resize(dim);
}
void band_matrix::resize(int dim)
{
    assert(dim>0);
    a.resize(dim);
    d.resize(dim);
    c.resize(dim);
}
int band_matrix::dim() const
{
    return d.size();
}


// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & band_matrix::operator () (int i, int j)
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert(k<2 && k>-2);
    if(k>0)return c[i];
    else if(k==0) return d[i];
    else return a[i];
}
double band_matrix::operator () (int i, int j) const
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    if(k>0)return c[i];
    else if(k==0) return d[i];
    else return a[i];
}



std::vector<double> band_matrix::trig_solve(const std::vector<double>& b) const
{
    assert( this->dim()==(int)b.size() );
    std::vector<double> x(this->dim());
    std::vector<double> g(this->dim());
    
    assert(d[0]!=0.0);
    x[0] = b[0]/d[0];
    double bet = d[0];
    for(int j=1;j<this->dim();j++){
        g[j] = c[j-1]/bet;
        bet = d[j] - (a[j]*g[j]);
        assert(bet != 0.0);
        x[j] = (b[j]-a[j]*x[j-1])/bet;
    }
    for(int j=this->dim()-2;j>=0;j--){
        x[j] -= g[j+1]*x[j+1];
    }
    
    return x;
}


// spline implementation
// -----------------------

void spline::set_boundary(spline::bd_type left, double left_value,
                          spline::bd_type right, double right_value,
                          bool force_linear_extrapolation)
{
    assert(n==0);          // set_points() must not have happened yet
    m_left=left;
    m_right=right;
    m_left_value=left_value;
    m_right_value=right_value;
    m_force_linear_extrapolation=force_linear_extrapolation;
}


void spline::set_points(const double *x,
                        const double *y,
                        const size_t num
                        )
{
    assert(num>2);
    m_x=x;
    m_y=y;
    n=num;
    // TODO: maybe sort x and y, rather than returning an error
    for(int i=0; i<n-1; i++) {
        assert(m_x[i]<m_x[i+1]);
    }

    
        // setting up the matrix and right hand side of the equation system
        // for the parameters b[]
        band_matrix A(n);
        std::vector<double>  rhs(n);
        for(int i=1; i<n-1; i++) {
            A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
            A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
            A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
            rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        }
        // boundary conditions
        if(m_left == spline::bd_type::second_deriv) {
            // 2*b[0] = f''
            A(0,0)=2.0;
            A(0,1)=0.0;
            rhs[0]=m_left_value;
        } else{
            // c[0] = f', needs to be re-expressed in terms of b:
            // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
            A(0,0)=2.0*(x[1]-x[0]);
            A(0,1)=1.0*(x[1]-x[0]);
            rhs[0]=3.0*((y[1]-y[0])/(x[1]-x[0])-m_left_value);
        } 

        if(m_right == spline::bd_type::second_deriv) {
            // 2*b[n-1] = f''
            A(n-1,n-1)=2.0;
            A(n-1,n-2)=0.0;
            rhs[n-1]=m_right_value;
        } else{
            // c[n-1] = f', needs to be re-expressed in terms of b:
            // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
            // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
            A(n-1,n-1)=2.0*(x[n-1]-x[n-2]);
            A(n-1,n-2)=1.0*(x[n-1]-x[n-2]);
            rhs[n-1]=3.0*(m_right_value-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
        }

        // solve the equation system to obtain the parameters b[]
        //m_b=A.lu_solve(rhs);
        m_b=A.trig_solve(rhs);

        // calculate parameters a[] and c[] based on b[]
        m_a.resize(n);
        m_c.resize(n);
        for(int i=0; i<n-1; i++) {
            m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
            m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
                   - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
        }
    
    // for left extrapolation coefficients
    m_b0 = (m_force_linear_extrapolation==false) ? m_b[0] : 0.0;
    m_c0 = m_c[0];

    // for the right extrapolation coefficients
    // f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
    double h=x[n-1]-x[n-2];
    // m_b[n-1] is determined by the boundary condition
    m_a[n-1]=0.0;
    m_c[n-1]=3.0*m_a[n-2]*h*h+2.0*m_b[n-2]*h+m_c[n-2];   // = f'_{n-2}(x_{n-1})
    if(m_force_linear_extrapolation==true)
        m_b[n-1]=0.0;
}

double spline::operator() (double x) const
{
    assert(n>2);
    // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
    auto it=std::lower_bound(m_x,m_x+n,x);
    //int idx=std::max( int(it-m_x)-1, 0);
    if(it!=m_x)it--;
    int idx = std::distance(m_x,it);
    double mx = *it;
    double h=x-mx;
    double interpol;
    if(x<m_x[0]) {
        // extrapolation to the left
        interpol=(m_b0*h + m_c0)*h + m_y[0];
    } else if(x>m_x[n-1]) {
        // extrapolation to the right
        interpol=(m_b[n-1]*h + m_c[n-1])*h + m_y[n-1];
    } else {
        // interpolation
        interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
    }
    return interpol;
}

double spline::deriv(int order, double x) const
{
    assert(order>0);
    // find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
    auto it=std::lower_bound(m_x,m_x+n,x);
    //int idx=std::max( int(it-m_x)-1, 0);
    if(it!=m_x)it--;
    int idx = std::distance(m_x,it);
    double mx = *it;
    double h=x-mx;
    double interpol;
    if(x<m_x[0]) {
        // extrapolation to the left
        switch(order) {
        case 1:
            interpol=2.0*m_b0*h + m_c0;
            break;
        case 2:
            interpol=2.0*m_b0*h;
            break;
        default:
            interpol=0.0;
            break;
        }
    } else if(x>m_x[n-1]) {
        // extrapolation to the right
        switch(order) {
        case 1:
            interpol=2.0*m_b[n-1]*h + m_c[n-1];
            break;
        case 2:
            interpol=2.0*m_b[n-1];
            break;
        default:
            interpol=0.0;
            break;
        }
    } else {
        // interpolation
        switch(order) {
        case 1:
            interpol=(3.0*m_a[idx]*h + 2.0*m_b[idx])*h + m_c[idx];
            break;
        case 2:
            interpol=6.0*m_a[idx]*h + 2.0*m_b[idx];
            break;
        case 3:
            interpol=6.0*m_a[idx];
            break;
        default:
            interpol=0.0;
            break;
        }
    }
    return interpol;
}



} // namespace tk


