/**
 * @file DGFEMSpace1D.h
 * @brief
 * @author Duan Junming, duanjm@pku.edu.cn
 * @version 1.0
 * @date 2018-04-16
 */

#ifndef DGFEMSPACE1D_H
#define DGFEMSPACE1D_H

#include "Eigen/Dense"
#include "para.h"
#include "BasFun.h"
#include "Quadrature.h"
#include "interval_crd_trs.h"
#include "engine.h"

class DGFEMSpace1D {
  private:
    u_int Nx;
    double xl, xr;
    double h;
    VEC<double> mesh;
    TemplateQuadrature TemQuad_Gauss;
    VEC<Quadrature> QUADINFO;
    VEC<double> mu;
    VEC<double> wgt;
    double sum_wm;
    EMAT M1, Mx, Fx_plus, Fx_minus, K1, Kx, F0, MM;//universal matrices
    VEC<VEC<double>> I_av;
    SOL I;//dimension: Nx*M*K
    ST_ele W, W_uw;
    BM bml, bmr;
    double BD_L, BD_R;
    EVEC FLUX;
    EMAT A;
    EVEC rhs;
    Eigen::ColPivHouseholderQR<EMAT> solver;

  public:
    DGFEMSpace1D(u_int Nx, double xl, double xr);
    void BuildQuad(u_int np);
    void st2coe(const u_int q, u_int& qx, u_int& qt);
    void BuildUniMAT();
    void Projection(u_int cell, func I0, double t, VEC<double>&);
    void Projection(const SOL& I, VEC<VEC<double>>&);
    EVEC Composition(const SOL&, u_int cell, double x);
    void init(func I0);
    double cal_dt(const SOL&, const double);
    double RAD_BE_unsteady(VEC<VEC<double>>&, SOL&,
        func_para, func_para, src, const double, const double, func, func);
    void Boundary(func BL, func BR, double& BD_L, double& BD_R, const double mu, const double t, const double dt);
    void solve_leqn(const EMAT&, const EVEC&, EVEC&);
    void STDG_predictor(const SOL&, const double dt, const u_int, const u_int, ST_ele& W);
    double flux_plus(const ST_ele& w, double dt);
    double flux_minus(const ST_ele& w, double dt);
    double source(const ST_ele& w, src q, u_int cell, u_int m, const double t, const double dt);
    void run_unsteady(func_para, func_para, src, func, func, double t_end);
    double cal_norm_I(const SOL& s1, const SOL& s2, int n);
    double cal_err(const SOL& s1, int n, double t_end);
    void print_DG_coe(std::ostream&);
    void print_solution_integral(std::ostream&);
    void print_solution_average(std::ostream&);
    void print_sol_ex_integral(std::ostream& os, func, const double t);
    int plot_I(engine*, const VEC<double>& mesh, const SOL& In);
    void scaling_limiter(const u_int, EVEC& u);
    void Reconstruct(VEC<VEC<double>>& I_av, SOL& I);
};

#endif //DGFEMSPACE1D_H

