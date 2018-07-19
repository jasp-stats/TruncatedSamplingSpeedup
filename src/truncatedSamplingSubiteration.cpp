
#include <Rcpp.h>
#include <mpfr.h>

using namespace Rcpp;

mpfr_t BigHi, BigLo, BigA, BigLoA, BigMinusZi, BigHighBound, BigUni0, BigUni1, BigResult;

// [[Rcpp::export]]
void initBigNumbers(int bitsPrecision)
{   
    auto setPrec = [&](mpfr_t & setMyPrecision){ mpfr_init2(setMyPrecision, bitsPrecision); };
    
    setPrec(BigHi);
    setPrec(BigLo);
    setPrec(BigA);
    setPrec(BigLoA);
    setPrec(BigMinusZi);
    setPrec(BigHighBound);
    setPrec(BigUni0);
    setPrec(BigUni1);
    setPrec(BigResult);
}


// [[Rcpp::export]]
double truncatedSamplingSubiteration(double uniformSample0, double uniformSample1, double minusZi, double Lo, double ai, bool thereIsAHigherBound, double theHigherBound ) 
{
    mpfr_set_d(BigUni0,         uniformSample0, MPFR_RNDD);
    mpfr_set_d(BigUni1,         uniformSample1, MPFR_RNDD);
    mpfr_set_d(BigMinusZi,      minusZi,        MPFR_RNDD);
    mpfr_set_d(BigLo,           Lo,             MPFR_RNDD);
    mpfr_set_d(BigA,            ai,             MPFR_RNDD);
    mpfr_set_d(BigHighBound,    theHigherBound, MPFR_RNDD);

    mpfr_pow(BigLoA, BigLo, BigA, MPFR_RNDD);
    mpfr_exp(BigHi, BigMinusZi, MPFR_RNDD);
    mpfr_mul(BigHi, BigHi, BigUni0, MPFR_RNDD);
    mpfr_log(BigHi, BigHi, MPFR_RNDD);
    mpfr_mul_d(BigHi, BigHi, -1.0, MPFR_RNDD);

    if(thereIsAHigherBound)
        mpfr_min(BigHi, BigHi, BigHighBound, MPFR_RNDD);


    mpfr_pow(BigResult, BigHi, BigA, MPFR_RNDD); //Hi^a[i]
    mpfr_sub(BigResult, BigResult, BigLoA, MPFR_RNDD); //(Hi^a[i]) - (Lo^a[i])
    mpfr_mul(BigResult, BigResult, BigUni1, MPFR_RNDD); //runif(1) * ((Hi^a[i]) - (Lo^a[i]))
    mpfr_add(BigResult, BigResult, BigLoA, MPFR_RNDD); //(runif(1) * ((Hi^a[i]) - (Lo^a[i]))) + (Lo^a[i])

    mpfr_d_div(BigA, 1.0, BigA, MPFR_RNDD); // 1/a[i]
    mpfr_pow(BigResult, BigResult, BigA, MPFR_RNDD); //((runif(1) * ((Hi^a[i]) - (Lo^a[i]))) + (Lo^a[i])) ^ (1/a[i])

    double result = mpfr_get_d(BigResult, MPFR_RNDD);

    return result;
}
