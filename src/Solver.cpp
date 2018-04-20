#include "DGFEMSpace1D.h"

void DGFEMSpace1D::solve_leqn(const EMAT& A, const EVEC& rhs, EVEC& u) {
  solver.compute(A);
  u = solver.solve(rhs);
}

void DGFEMSpace1D::STDG_predictor(const SOL& I, const double dt, const u_int cell, const u_int m, ST_ele& W) {
  double hi = mesh[cell+1] - mesh[cell];
  double para = -1;
  A = hi*K1 + (dt*mu[m])*Kx - para*(hi*dt)*MM;//all
  //A = hi*K1 + (dt*mu[m])*Kx;//without source
  //A = hi*K1 - (hi*dt)*MM;//without advection
  rhs = hi*F0*I[cell][m];
  solve_leqn(A, rhs, W);
}

double DGFEMSpace1D::flux_plus(const ST_ele& w, double dt) {
  u_int px, pt;
  double f(0);
  VEC<double> wg = TemQuad_Gauss.weight();
  EVEC p1 = Lagrange_Poly(1);
  u_int ind(0);
  for(u_int g = 0; g < K; ++g) {
    for(u_int p = 0; p < K*K; ++p) {
      st2coe(p, px, pt);
      if(pt == g) {
        ind++;
        f += wg[g]*w[p]*p1[px];
      }
    }
  }
  return f*dt/2.0;
}

double DGFEMSpace1D::flux_minus(const ST_ele& w, double dt) {
  u_int px, pt;
  double f(0);
  VEC<double> wg = TemQuad_Gauss.weight();
  EVEC p1 = Lagrange_Poly(-1);
  for(u_int g = 0; g < K; ++g) {
    for(u_int p = 0; p < K*K; ++p) {
      st2coe(p, px, pt);
      if(pt == g) {
        f += wg[g]*w[p]*p1[px];
      }
    }
  }
  return f*dt/2.0;
}

double DGFEMSpace1D::source(const ST_ele& w, src q, u_int cell, u_int m,
    const double t, const double dt) {
  u_int px, pt;
  double s(0);
  double hi = mesh[cell+1] - mesh[cell];
  VEC<double> wg = QUADINFO[cell].weight();
  VEC<double> xg = QUADINFO[cell].points();
  VEC<double> lv(2), gv(2);
  lv[0] = -1, lv[1] = 1;
  gv[0] = t, gv[1] = t+dt;

  double tg(0);
  double para(-1);
  for(u_int p = 0; p < K*K; ++p) {
    st2coe(p, px, pt);
    local_to_global(xg[pt], lv, gv, &tg);
    s += wg[px]*wg[pt]*para*q(mu[m], w[p], xg[px], tg);
  }
  return s*hi*dt/4.0;
}

double DGFEMSpace1D::RAD_BE_unsteady(VEC<VEC<double>>& I_av, SOL& I,
    func_para sigma_t, func_para sigma_s, src q,
    const double t, const double dt, func BL, func BR) {
  VEC<ST_ele> W_list;
  W_list.resize(Nx);
  double F1, F2;
  for(u_int m = 0; m < M; ++m) {//for each light direction
    Boundary(BL, BR, BD_L, BD_R, mu[m], t, dt);
    if(mu[m] >= 0) {
      for(u_int i = 0; i < Nx; ++i) {
        STDG_predictor(I, dt, i, m, W_list[i]);

      }
      for(u_int i = 0; i < Nx; ++i) {//for each cell, update I_{m,i}^{l}
        double hi = mesh[i+1] - mesh[i];
        if(i == 0) {
          F1 = flux_plus(W_list[i], dt);
          F2 = flux_plus(W_list[Nx-1], dt);
          I_av[i][m] -= ( mu[m]*F1 - mu[m]*F2 )/hi;
        }
        else {
          F1 = flux_plus(W_list[i], dt);
          F2 = flux_plus(W_list[i-1], dt);
          I_av[i][m] -= ( mu[m]*F1 - mu[m]*F2 )/hi;
        }
        I_av[i][m] += source(W_list[i], q, i, m, t, dt)/hi;

        //double x0, x1, x2, t0, t1, t2;
        //VEC<double> lv(2), gv(2);
        //lv[0] = -1, lv[1] = 1;
        //gv[0] = mesh[0], gv[1] = mesh[1];
        //VEC<double> pg = TemQuad_Gauss.points();
        //VEC<double> wg = TemQuad_Gauss.weight();
        //local_to_global(pg[0], lv, gv, &x0);
        //local_to_global(pg[1], lv, gv, &x1);
        //local_to_global(pg[2], lv, gv, &x2);
        //gv[0] = 0, gv[1] = dt;
        //local_to_global(pg[0], lv, gv, &t0);
        //local_to_global(pg[1], lv, gv, &t1);
        //local_to_global(pg[2], lv, gv, &t2);

        //std::cout << "dt: " << dt << std::endl;
        //std::cout << "hi: " << mesh[1]-mesh[0] << std::endl;
        //std::cout << "Recon cha: " << std::endl;
        //std::cout
        //<< sin(2.0*M_PI*(x0)) - I[i][m][0] << " "
        //<< sin(2.0*M_PI*(x1)) - I[i][m][1] << " "
        //<< sin(2.0*M_PI*(x2)) - I[i][m][2] << std::endl;
        //std::cout << "I_av: " << I_av[0][0] << std::endl;
        //std::cout << "points x0, x1, t0, t1\n" << x0 << " " << x1
          //<< " " << t0 << " " << t1 << std::endl;
        //std::cout << "theorecal val\n"
          //<< sin(2.0*M_PI*(x0-t0)) << "\n"
          //<< sin(2.0*M_PI*(x1-t0)) << "\n"
          //<< sin(2.0*M_PI*(x2-t0)) << "\n"
          //<< sin(2.0*M_PI*(x0-t1)) << "\n"
          //<< sin(2.0*M_PI*(x1-t1)) << "\n"
          //<< sin(2.0*M_PI*(x2-t1)) << "\n"
          //<< sin(2.0*M_PI*(x0-t2)) << "\n"
          //<< sin(2.0*M_PI*(x1-t2)) << "\n"
          //<< sin(2.0*M_PI*(x2-t2)) << std::endl;
        //std::cout << "dg predicted val\n" <<
          //W_list[0] << std::endl;
        //std::cout << "cha\n"
          //<< sin(2.0*M_PI*(x0-t0)) - W_list[0][0] << "\n"
          //<< sin(2.0*M_PI*(x1-t0)) - W_list[0][1] << "\n"
          //<< sin(2.0*M_PI*(x2-t0)) - W_list[0][2] << "\n"
          //<< sin(2.0*M_PI*(x0-t1)) - W_list[0][3] << "\n"
          //<< sin(2.0*M_PI*(x1-t1)) - W_list[0][4] << "\n"
          //<< sin(2.0*M_PI*(x2-t1)) - W_list[0][5] << "\n"
          //<< sin(2.0*M_PI*(x0-t2)) - W_list[0][6] << "\n"
          //<< sin(2.0*M_PI*(x1-t2)) - W_list[0][7] << "\n"
          //<< sin(2.0*M_PI*(x2-t2)) - W_list[0][8] << std::endl;

        //std::cout << "F1: " << F1 << std::endl;
        //std::cout << "F2: " << F2 << std::endl;
        //std::cout << "src: " <<
          //source(W_list[i], q, i, m, dt) << std::endl;
        //std::cout << "Lagrange_Poly(1)\n" <<
          //Lagrange_Poly(1) << std::endl;
        //std::cout << "flux+\n" <<
          //((W_list[0][0]*wg[0]+W_list[0][3]*wg[1]+W_list[0][6]*wg[2])*Lagrange_Poly(1)[0] +
          //(W_list[0][1]*wg[0]+W_list[0][4]*wg[1]+W_list[0][7]*wg[2])*Lagrange_Poly(1)[1] +
          //(W_list[0][2]*wg[0]+W_list[0][5]*wg[1]+W_list[0][8]*wg[2])*Lagrange_Poly(1)[2])*dt/2.0
          //<< std::endl;
        //std::cout << "source\n" <<
          //-(W_list[0][0]*wg[0]*wg[0]+
          //W_list[0][1]*wg[0]*wg[1]+
          //W_list[0][2]*wg[0]*wg[2]+
          //W_list[0][3]*wg[1]*wg[0]+
          //W_list[0][4]*wg[1]*wg[1]+
          //W_list[0][5]*wg[1]*wg[2]+
          //W_list[0][6]*wg[2]*wg[0]+
          //W_list[0][7]*wg[2]*wg[1]+
          //W_list[0][8]*wg[2]*wg[2])*dt*hi/4.0
          //<< std::endl;
        //abort();

      }//end i
    }//end ->
    else {
      u_int i = Nx-1;
      while(i >= 0) {
        if(i-- == 0) break;
      }//end i
    }//end <-
  }//end M

  return 1;
}


/*double DGFEMSpace1D::RAD_BE_unsteady(const SOL& In, SOL& I,*/
    //func_para sigma_t, func_para sigma_s, func q,
    //const double t, const double dt, SOL& I_new, func BL, func BR) {
  //for(u_int m = 0; m < M; ++m) {//for each light direction
    //Boundary(BL, BR, BD_L, BD_R, mu[m], t, dt);
    //if(mu[m] >= 0) {
      //for(u_int i = 0; i < Nx; ++i) {//for each cell, update I_{m,i}^{l}
        //std::cout << "########i: " << i << " ########" << std::endl;
        //double hi = mesh[i+1] - mesh[i];
        //STDG_predictor(In, dt, i, m, W);

        //A = hi*M1;
        //if(i == 0) {
          //FLUX = BD_L*Lagrange_Poly(-1);
          //rhs = hi*M1*In[i][m] + (mu[m]*dt)*Mx*W - dt*Fx_plus*W + FLUX;
        //}
        //else {
          //rhs = hi*M1*In[i][m] + (mu[m]*dt)*Mx*W - dt*Fx_plus*W + dt*Fx_minus*W_uw;
        //}
        //solve_leqn(A, rhs, I_new[i][m]);

        ////do limiter
        ////scaling_limiter(i, I_new[i][m]);
        ////WENO_limiter(i, m, I_new);

        //W_uw = W;//update the upwindign W

      //}//end i
    //}//end ->
    //else {
      //u_int i = Nx-1;
      //while(i >= 0) {
        //if(i-- == 0) break;
      //}//end i
    //}//end <-
  //}//end M

  ////abort();
  //return 1;
//}

void DGFEMSpace1D::run_unsteady(func_para sigma_t, func_para sigma_s, src q, func BL, func BR, double t_end) {
  std::cout.precision(16);
  std::cout << std::showpos;
  std::cout.setf(std::ios::scientific);
  double t(0), dt(0);

  while (t < t_end) {
    dt = cal_dt(I, t);
    dt = std::min(dt, t_end-t);
    Reconstruct(I_av, I);
    RAD_BE_unsteady(I_av, I, sigma_t, sigma_s, q, t, dt, BL, BR);
    t += dt;
    std::cout << "t: " << t << ", dt: " << dt << std::endl;
  }
}

