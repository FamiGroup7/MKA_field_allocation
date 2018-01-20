#if !defined COCG_LLT_SMOOTH_H_INCLUDED
#define COCG_LLT_SMOOTH_H_INCLUDED

#include <cstdlib>
using namespace std;

class CGM_LLT
{
public:
    void init(size_t * gi_s, size_t * gj_s, double * di_s,
              double * gg_s, size_t n_s);
    void solve(double * solution, double * rp_s, double eps, size_t max_iter);

    CGM_LLT();
    ~CGM_LLT();
private:
    void make_LLT_decomposition();
    void mul_matrix(const double * f, double * x) const;
    void solve_L(const double * f, double * x) const;
    void solve_LT(double * f, double * x) const;
    void solve_LLT(const double * f, double * x) const;
    double dot_prod_nocj(const double * a, const double * b) const;
    double dot_prod_self(const double * a) const;
    double dot_prod_real(const double * a, const double * b) const;
    bool is_fpu_error(double x) const;

    size_t n;
    size_t * gi, * gj;
    double * di, * gg, * rp, * r, * x0, * z, * p, * s, * xs, * rs;
    double * L_di, * L_gg;
};

#endif // COCG_LLT_SMOOTH_H_INCLUDED
