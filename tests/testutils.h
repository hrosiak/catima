namespace catima{

class approx
{
public:
    explicit approx ( double magnitude, double eps=1.0)
    : epsilon_  { eps }
    , magnitude_{ magnitude } {}

    approx( approx const & other ) = default;

    approx operator()( double magnitude, double eps = 1.0 )
    {
        approx approx ( magnitude, eps);
        return approx;
    }

    double magnitude() const { return magnitude_; }

    approx & epsilon( double epsilon ) { epsilon_ = epsilon; return *this; }
    approx & R( double relative ) { epsilon_ = relative*magnitude_; return *this; }


    friend bool operator == ( double lhs, approx const & rhs )
    {
        return std::abs( lhs - rhs.magnitude_) < rhs.epsilon_;
    }

    friend bool operator == ( approx const & lhs, double rhs ) { return  operator==( rhs, lhs ); }
    friend bool operator != ( double lhs, approx const & rhs ) { return !operator==( lhs, rhs ); }
    friend bool operator != ( approx const & lhs, double rhs ) { return !operator==( rhs, lhs ); }

    friend bool operator <= ( double lhs, approx const & rhs ) { return lhs < rhs.magnitude_ || lhs == rhs; }
    friend bool operator <= ( approx const & lhs, double rhs ) { return lhs.magnitude_ < rhs || lhs == rhs; }
    friend bool operator >= ( double lhs, approx const & rhs ) { return lhs > rhs.magnitude_ || lhs == rhs; }
    friend bool operator >= ( approx const & lhs, double rhs ) { return lhs.magnitude_ > rhs || lhs == rhs; }

//private:
    double epsilon_;
    double magnitude_;
};

std::ostream & operator<<(std::ostream &os, approx const &a){
    using lest::to_string;
    return os<<to_string(a.magnitude_)<<" +- "<<a.epsilon_;
}


bool rdiff(double a, double b,double eps){
    if(fabs((a-b)/fabs(b))<eps){
      return true;
    }
    else{
      std::cout<<"\033[1;31m"<<a<<" == "<<b<<"\033[0m"<<std::endl;
      return false;
    }
      
}

bool diff(double a, double b,double eps){
    if(fabs((a-b))<eps){
      return true;
    }
    else{
      std::cout<<"\033[1;31m"<<a<<" == "<<b<<"\033[0m"<<std::endl;
      return false;
    }
      
}

}
