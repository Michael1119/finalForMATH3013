#include "Parameters.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

using namespace std;

const double pi = 4 * atan(1);

class PDE{
public:
    double get_a(){return a;}
protected:
    PDE(){}
    ~PDE(){}
    double a;
};

class heat1D : public PDE{
public:
    heat1D(double constant){a = constant;}
    ~heat1D(){}
};

class heat2D : public PDE{
public:
    heat2D(double constant){a = constant;}
    ~heat2D(){}
};

class wave1D : public PDE{
public:
    wave1D(double constant){a = constant;}
    ~wave1D(){}
    double f(double x){
        return sin(2 * pi * x);
    }
    double g(double x){
        return 2 * pi * sin(2 * pi * x);
    }
};

class PDEsolver{
protected:
    PDEsolver(){}
    ~PDEsolver(){}
    Parameters par;
    int nt;
    int nx;
    int ny;
    double lambda;
    void output(fstream& file, const double& t, const vector<double>& u){
        double x = par.dx;
        for(int i = 0; i < nx - 1; ++i){
            file << t << " " << x << " " << u[i] << endl;
            x += par.dx;
        }
        file << endl << endl;
    }
    void output(fstream& file, const double& t, const vector<vector<double>>& u){
        double x = par.dx;
        double y = par.dx;
        for(int j = 0; j < ny - 1; ++j){
            x = par.dx;
            for(int i = 0; i < nx - 1; ++i){
                file << t << " " << x << " " << y << " " << u[i][j] << endl;
                x += par.dx;
            }
            y += par.dx;
        }
        file << endl << endl;
    }
};

class Explicit : public PDEsolver{
protected:
    Explicit(){}
    ~Explicit(){}
};

class Implicit : public PDEsolver{
protected:
    Implicit(){}
    ~Implicit(){}
    double diag;
    double offdiag;
    void tridiag(double A_diag, double A_offdiag, vector<double>& x, vector<double>& rhs){
        vector<double> b(x.size(), A_diag);
        vector<double> d = rhs;
        double w;
        for(int i = 1; i < x.size(); ++i){
            w = A_offdiag / b[i - 1];
            b[i] -= w * A_offdiag;
            d[i] -= w * d[i - 1];
        }
        x[x.size() - 1] = d[x.size() - 1] / b[x.size() - 1];
        for(int i = x.size() - 2; i >= 0; --i){
            x[i] = (d[i] - A_offdiag * x[i + 1]) / b[i];
        }
    }
    void conjgrad(double A_diag, double A_offdiag, vector<vector<double>>& x, vector<vector<double>>& rhs){
        vector<vector<double>> r = rhs;
        vector<vector<double>> Ax = x;
        vector<vector<double>> p = r;
        vector<vector<double>> Ap = p;
        double alpha;
        double beta;
        double rtr;
        double ptap;
        for(int j = 0; j < ny - 1; ++j){
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    Ax[i][j] = A_diag * x[i][j] + A_offdiag * x[i + 1][j];
                }else if(i == nx - 2){
                    Ax[i][j] = A_offdiag * x[i - 1][j] + A_diag * x[i][j];
                }else{
                    Ax[i][j] = A_offdiag * x[i - 1][j] + A_diag * x[i][j] + A_offdiag * x[i + 1][j];
                }
                if(j == 0){
                    Ax[i][j] += A_offdiag * x[i][j + 1]; 
                }else if(j == nx - 2){
                    Ax[i][j] += A_offdiag * x[i][j - 1];
                }else{
                    Ax[i][j] += A_offdiag * x[i][j - 1] + A_offdiag * x[i][j + 1];
                }
            }
        }
        for(int j = 0; j < ny - 1; ++j){
            for(int i = 0; i < nx - 1; ++i){
                r[i][j] -= Ax[i][j];
            }
        }
        p = r;
        while(true){
            for(int j = 0; j < ny - 1; ++j){
                for(int i = 0; i < nx - 1; ++i){
                    if(i == 0){
                        Ap[i][j] = A_diag * p[i][j] + A_offdiag * p[i + 1][j];
                    }else if(i == nx - 2){
                        Ap[i][j] = A_offdiag * p[i - 1][j] + A_diag * p[i][j];
                    }else{
                        Ap[i][j] = A_offdiag * p[i - 1][j] + A_diag * p[i][j] + A_offdiag * p[i + 1][j];
                    }
                    if(j == 0){
                        Ap[i][j] += A_offdiag * p[i][j + 1]; 
                    }else if(j == nx - 2){
                        Ap[i][j] += A_offdiag * p[i][j - 1];
                    }else{
                        Ap[i][j] += A_offdiag * p[i][j - 1] + A_offdiag * p[i][j + 1];
                    }
                }
            }
            rtr = dot(r, r);
            ptap = dot(p, Ap);
            alpha = rtr / ptap;
            for(int j = 0; j < ny - 1; ++j){
                for(int i = 0; i < nx - 1; ++i){
                    x[i][j] += alpha * p[i][j];
                    r[i][j] -= alpha * Ap[i][j];
                }
            }
            if(dot(r, r) < 1e-10){
                break;
            }
            beta = dot(r, r) / rtr;
            for(int j = 0; j < ny - 1; ++j){
                for(int i = 0; i < nx - 1; ++i){
                    p[i][j] = r[i][j] + beta * p[i][j];
                }
            }          
        }
    }
    double dot(vector<vector<double>>& A, vector<vector<double>>& B){
        double sum = 0;
        for(int j = 0; j < ny - 1; ++j){
            for(int i = 0; i < nx - 1; ++i){
                sum += A[i][j] * B[i][j];
            }
        }
        return sum;
    }
};

class ForwardEuler : public Explicit{
public:
    ForwardEuler(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
        ny = par.y_ny / par.dx;
    }
    ~ForwardEuler(){}
    void solve(heat1D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        vector<double> u(nx - 1 , par.u_t0);
        vector<double> v;
        fstream file;
        file.open("Forward Euler 1D.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    u[i] = lambda * par.u_nx + (1 - 2 * lambda) * v[i] + lambda * v[i + 1];
                }else if(i == nx - 2){
                    u[i] = lambda * v[i - 1] + (1 - 2 * lambda) * v[i] + lambda * par.u_nx;
                }else{
                    u[i] = lambda * v[i - 1] + (1 - 2 * lambda) * v[i] + lambda * v[i + 1];
                }
            }
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
    void solve(heat2D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        vector<vector<double>> u(nx - 1, vector<double>(ny - 1, par.u_t0));
        vector<vector<double>> v;
        fstream file;
        file.open("Forward Euler 2D.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            v = u;
            for(int j = 0; j < ny - 1; ++j){
                for(int i = 0; i < nx - 1; ++i){
                    if(i == 0){
                        u[i][j] = lambda * par.u_nx + (1 - 4 * lambda) * v[i][j] + lambda * v[i + 1][j];
                    }else if(i == nx - 2){
                        u[i][j] = lambda * v[i - 1][j] + (1 - 4 * lambda) * v[i][j] + lambda * par.u_nx;
                    }else{
                        u[i][j] = lambda * v[i - 1][j] + (1 - 4 * lambda) * v[i][j] + lambda * v[i + 1][j];
                    }
                    if(j == 0){
                        u[i][j] += lambda * par.u_nx + lambda * v[i][j + 1];
                    }else if(j == ny - 2){
                        u[i][j] += lambda * v[i][j - 1] + lambda * par.u_nx;
                    }else{
                        u[i][j] += lambda * v[i][j - 1] + lambda * v[i][j + 1];
                    }
                }
            }
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class BackwardEuler : public Implicit{
public:
    BackwardEuler(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
        ny = par.y_ny / par.dx;
    }
    ~BackwardEuler(){}
    void solve(heat1D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        diag = 1 + 2 * lambda;
        offdiag = - lambda;
        vector<double> u(nx - 1 , par.u_t0);
        vector<double> v;
        vector<double> b(nx - 1 , 0);
        fstream file;
        file.open("Backward Euler 1D.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    b[i] = v[i] + lambda * par.u_nx;
                }else if(i == nx - 2){
                    b[i] = v[i] + lambda * par.u_nx;
                }else{
                    b[i] = v[i];
                }
            }
            tridiag(diag, offdiag, u, b);
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
    void solve(heat2D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        diag = 1 + 4 * lambda;
        offdiag = - lambda;
        vector<vector<double>> u(nx - 1, vector<double>(ny - 1, par.u_t0));
        vector<vector<double>> v;
        vector<vector<double>> b = u;
        fstream file;
        file.open("Backward Euler 2D.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            v = u;
            for(int j = 0; j < ny - 1; ++j){
                for(int i = 0; i < nx - 1; ++i){
                    if(i == 0){
                        b[i][j] = v[i][j] + lambda * par.u_nx;
                    }else if(i == nx - 2){
                        b[i][j] = v[i][j] + lambda * par.u_nx;
                    }else{
                        b[i][j] = v[i][j];
                    }
                    if(j == 0){
                        b[i][j] += lambda * par.u_nx;
                    }else if(j == ny - 2){
                        b[i][j] += lambda * par.u_nx;
                    }
                }
            }
            conjgrad(diag, offdiag, u, b);
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class Richardson : public Explicit{
public:
    Richardson(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
    }
    ~Richardson(){}
    void solve(heat1D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        vector<double> u(nx - 1 , par.u_t0);
        vector<double> v;
        vector<double> w;
        fstream file;
        file.open("Richardson.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        v = u;
        for(int i = 0; i < nx - 1; ++i){
            if(i == 0){
                u[i] = lambda * par.u_nx + (1 - 2 * lambda) * v[i] + lambda * v[i + 1];
            }else if(i == nx - 2){
                u[i] = lambda * v[i - 1] + (1 - 2 * lambda) * v[i] + lambda * par.u_nx;
            }else{
                u[i] = lambda * v[i - 1] + (1 - 2 * lambda) * v[i] + lambda * v[i + 1];
            }
        }
        tk += par.dt;
        output(file, tk, u);
        for(int k = 1; k < nt; ++k){
            w = v;
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    u[i] = w[i] + 2 * lambda * (par.u_nx - 2 * v[i] + v[i + 1]);
                }else if(i == nx - 2){
                    u[i] = w[i] + 2 * lambda * (v[i - 1] - 2 * v[i] + par.u_nx);
                }else{
                    u[i] = w[i] + 2 * lambda * (v[i - 1] - 2 * v[i] + v[i + 1]);
                }
            }
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class DufortFrankel : public Explicit{
public:
    DufortFrankel(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
    }
    ~DufortFrankel(){}
    void solve(heat1D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        vector<double> u(nx - 1 , par.u_t0);
        vector<double> v;
        vector<double> w;
        fstream file;
        file.open("Dufort-Frankel.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        v = u;
        for(int i = 0; i < nx - 1; ++i){
            if(i == 0){
                u[i] = lambda * par.u_nx + (1 - 2 * lambda) * v[i] + lambda * v[i + 1];
            }else if(i == nx - 2){
                u[i] = lambda * v[i - 1] + (1 - 2 * lambda) * v[i] + lambda * par.u_nx;
            }else{
                u[i] = lambda * v[i - 1] + (1 - 2 * lambda) * v[i] + lambda * v[i + 1];
            }
        }
        tk += par.dt;
        output(file, tk, u);
        for(int k = 1; k < nt; ++k){
            w = v;
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    u[i] = ((1 - 2 * lambda) * w[i] + 2 * lambda * (par.u_nx + v[i + 1])) / (1 + 2 * lambda);
                }else if(i == nx - 2){
                    u[i] = ((1 - 2 * lambda) * w[i] + 2 * lambda * (v[i - 1] + par.u_nx)) / (1 + 2 * lambda);
                }else{
                    u[i] = ((1 - 2 * lambda) * w[i] + 2 * lambda * (v[i - 1] + v[i + 1])) / (1 + 2 * lambda);
                }
            }
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class CrankNicolson : public Implicit{
public:
    CrankNicolson(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
        ny = par.y_ny / par.dx;
    }
    ~CrankNicolson(){}
    void solve(heat1D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        diag = 1 + lambda;
        offdiag = - lambda / 2;
        vector<double> u(nx - 1 , par.u_t0);
        vector<double> v;
        vector<double> b(nx - 1 , 0);
        fstream file;
        file.open("Crank-Nicolson 1D.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    b[i] = lambda / 2 * par.u_nx + (1 - lambda) * v[i] + lambda / 2 * v[i + 1];
                }else if(i == nx - 2){
                    b[i] = lambda / 2 * v[i - 1] + (1 - lambda) * v[i] + lambda / 2 * par.u_nx;
                }else{
                    b[i] = lambda / 2 * v[i - 1] + (1 - lambda) * v[i] + lambda / 2 * v[i + 1];
                }
                if(i == 0){
                    b[i] += lambda / 2 * par.u_nx;
                }else if(i == nx - 2){
                    b[i] += lambda / 2 * par.u_nx;
                }
            }
            tridiag(diag, offdiag, u, b);
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
    void solve(heat2D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        diag = 1 + 2 * lambda;
        offdiag = - lambda / 2;
        vector<vector<double>> u(nx - 1, vector<double>(ny - 1, par.u_t0));
        vector<vector<double>> v;
        vector<vector<double>> b = u;
        fstream file;
        file.open("Crank-Nicolson 2D.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            v = u;
            for(int j = 0; j < ny - 1; ++j){
                for(int i = 0; i < nx - 1; ++i){
                    if(i == 0){
                        b[i][j] = lambda / 2 * par.u_nx + (1 - 2 * lambda) * v[i][j] + lambda / 2 * v[i + 1][j] + lambda / 2 * par.u_nx;
                    }else if(i == nx - 2){
                        b[i][j] = lambda / 2 * v[i - 1][j] + (1 - 2 * lambda) * v[i][j] + lambda / 2 * par.u_nx + lambda / 2 * par.u_nx;
                    }else{
                        b[i][j] = lambda / 2 * v[i - 1][j] + (1 - 2 * lambda) * v[i][j] + lambda / 2 * v[i + 1][j];
                    }
                    if(j == 0){
                        b[i][j] += lambda / 2 * par.u_nx + lambda / 2 * v[i][j + 1] + lambda / 2 * par.u_nx;
                    }else if(j == nx - 2){
                        b[i][j] += lambda / 2 * v[i][j - 1] + lambda / 2 * par.u_nx + lambda / 2 * par.u_nx;
                    }else{
                        b[i][j] += lambda / 2 * v[i][j - 1] + lambda / 2 * v[i][j + 1];
                    }
                }
            }
            conjgrad(diag, offdiag, u, b);
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class ADI : public Implicit{
public:
    ADI(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
        ny = par.y_ny / par.dx;
    }
    ~ADI(){}
    void solve(heat2D heat){
        lambda = heat.get_a() * heat.get_a() * par.dt / par.dx / par.dx;
        diag = 1 + lambda;
        offdiag = - lambda / 2;
        vector<vector<double>> u(nx - 1, vector<double>(ny - 1, par.u_t0));
        vector<vector<double>> v;
        vector<double> u_x(nx - 1, 0);
        vector<double> b_x(nx - 1, 0);
        vector<double> u_y(ny - 1, 0);
        vector<double> b_y(ny - 1, 0);
        fstream file;
        file.open("ADI.dat", std::ios::in|std::ios::out|std::ios::trunc);
        double tk = 0;
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            // implicit in x and explicit in y
            v = u;
            for(int j = 0; j < ny - 1; ++j){
                for(int i = 0; i < nx - 1; ++i){
                    if(j == 0){
                        b_x[i] = lambda / 2 * par.u_nx + (1 - lambda) * v[i][j] + lambda / 2 * v[i][j + 1];
                    }else if(j == ny - 2){
                        b_x[i] = lambda / 2 * v[i][j - 1] + (1 - lambda) * v[i][j] + lambda / 2 * par.u_nx;
                    }else{
                        b_x[i] = lambda / 2 * v[i][j - 1] + (1 - lambda) * v[i][j] + lambda / 2 * v[i][j + 1];
                    }
                    if(i == 0){
                        b_x[i] += lambda / 2 * par.u_nx;
                    }else if(i == nx - 2){
                        b_x[i] += lambda / 2 * par.u_nx;
                    }
                }
                tridiag(diag, offdiag, u_x, b_x);
                for(int i = 0; i < nx - 1; ++i){
                    u[i][j] = u_x[i];
                }
            }
            // explicit in x and implicit in y
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                for(int j = 0; j < ny - 1; ++j){
                    if(i == 0){
                        b_y[j] = lambda / 2 * par.u_nx + (1 - lambda) * v[i][j] + lambda / 2 * v[i + 1][j];
                    }else if(i == nx - 2){
                        b_y[j] = lambda / 2 * v[i - 1][j] + (1 - lambda) * v[i][j] + lambda / 2 * par.u_nx;
                    }else{
                        b_y[j] = lambda / 2 * v[i - 1][j] + (1 - lambda) * v[i][j] + lambda / 2 * v[i + 1][j];
                    }
                    if(j == 0){
                        b_y[j] += lambda / 2 * par.u_nx;
                    }else if(j == ny - 2){
                        b_y[j] += lambda / 2 * par.u_nx;
                    }
                }
                tridiag(diag, offdiag, u_y, b_y);
                for(int j = 0; j < ny - 1; ++j){
                    u[i][j] = u_y[j];
                }
            }
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class ExplicitFiniteDifference : public Explicit{
public:
    ExplicitFiniteDifference(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
    }
    ~ExplicitFiniteDifference(){}
    void solve(wave1D wave){
        lambda = wave.get_a() * par.dt / par.dx;
        double lambda2 = lambda * lambda;
        vector<double> u(nx - 1, 0);
        vector<double> v;
        vector<double> w;
        vector<double> x(nx - 1, 0);
        fstream file;
        file.open("Explicit Finite Difference.dat", std::ios::in|std::ios::out|std::ios::trunc);
        for(int i = 0; i < nx - 1; ++i){
            x[i] = (i + 1) * par.dx;
            u[i] = wave.f(x[i]);
        }
        double tk = 0;
        output(file, tk, u);
        v = u;
        for(int i = 0; i < nx - 1; ++i){
            if(i == 0){
                u[i] = (1 - lambda2) * wave.f(x[i]) + lambda2 / 2 * wave.f(x[i + 1]) + par.dt * wave.g(x[i]);
            }else if(i == nx - 2){
                u[i] = lambda2 / 2 * wave.f(x[i - 1]) + (1 - lambda2) * wave.f(x[i]) + par.dt * wave.g(x[i]);
            }else{
                u[i] = lambda2 / 2 * wave.f(x[i - 1]) + (1 - lambda2) * wave.f(x[i]) + lambda2 / 2 * wave.f(x[i + 1]) + par.dt * wave.g(x[i]);
            }
        }
        tk += par.dt;
        output(file, tk, u);
        for(int k = 1; k < nt; ++k){
            w = v;
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    u[i] = 2 * (1 - lambda2) * v[i] + lambda2 * v[i + 1] - w[i];
                }else if(i == nx - 2){
                    u[i] = lambda2 * v[i - 1] + 2 * (1 - lambda2) * v[i] - w[i];
                }else{
                    u[i] = lambda2 * v[i - 1] + 2 * (1 - lambda2) * v[i] + lambda2 * v[i + 1] - w[i];
                }
            }
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class ImplicitFiniteDifference : public Implicit{
public:
    ImplicitFiniteDifference(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
    }
    ~ImplicitFiniteDifference(){}
    void solve(wave1D wave){
        lambda = wave.get_a() * par.dt / par.dx;
        double lambda2 = lambda * lambda;
        diag = 1 + lambda2 / 2;
        offdiag = - lambda2 / 4;
        vector<double> u(nx - 1, 0);
        vector<double> v;
        vector<double> w;
        vector<double> x(nx - 1, 0);
        vector<double> b(nx - 1, 0);
        fstream file;
        file.open("Implicit Finite Difference.dat", std::ios::in|std::ios::out|std::ios::trunc);
        for(int i = 0; i < nx - 1; ++i){
            x[i] = (i + 1) * par.dx;
            u[i] = wave.f(x[i]);
        }
        double tk = 0;
        output(file, tk, u);
        v = u;
        for(int i = 0; i < nx - 1; ++i){
            if(i == 0){
                u[i] = (1 - lambda2) * wave.f(x[i]) + lambda2 / 2 * wave.f(x[i + 1]) + par.dt * wave.g(x[i]);
            }else if(i == nx - 2){
                u[i] = lambda2 / 2 * wave.f(x[i - 1]) + (1 - lambda2) * wave.f(x[i]) + par.dt * wave.g(x[i]);
            }else{
                u[i] = lambda2 / 2 * wave.f(x[i - 1]) + (1 - lambda2) * wave.f(x[i]) + lambda2 / 2 * wave.f(x[i + 1]) + par.dt * wave.g(x[i]);
            }
        }
        tk += par.dt;
        output(file, tk, u);
        for(int k = 1; k < nt; ++k){
            w = v;
            v = u;
            for(int i = 0; i < nx - 1; ++i){
                if(i == 0){
                    b[i] = lambda2 * ((- 2 * v[i] + v[i + 1]) / 2 + (- 2 * w[i] + w[i + 1]) / 4) + 2 * v[i] - w[i];
                }else if(i == nx - 2){
                    b[i] = lambda2 * ((v[i - 1] - 2 * v[i]) / 2 + (w[i - 1] - 2 * w[i]) / 4) + 2 * v[i] - w[i];
                }else{
                    b[i] = lambda2 * ((v[i - 1] - 2 * v[i] + v[i + 1]) / 2 + (w[i - 1] - 2 * w[i] + w[i + 1]) / 4) + 2 * v[i] - w[i];
                }
            }
            tridiag(diag, offdiag, u, b);
            tk += par.dt;
            output(file, tk, u);
        }
        file.close();
    }
};

class Exact : public PDEsolver{
public:
    Exact(Parameters parameters){
        par = parameters;
        nt = par.t_nt / par.dt;
        nx = par.x_nx / par.dx;
    }
    ~Exact(){}
    void solve(wave1D wave){
        vector<double> u(nx - 1, 0);
        fstream file;
        file.open("Exact solution.dat", std::ios::in|std::ios::out|std::ios::trunc);
        vector<double> x(nx - 1, 0);
        double tk = 0;
        for(int i = 0; i < nx - 1; ++i){
            x[i] = (i + 1) * par.dx;
            u[i] = sin(2 * pi * x[i]) * (cos(2 * pi * tk) + sin(2 * pi * tk));
        }
        output(file, tk, u);
        for(int k = 0; k < nt; ++k){
            tk += par.dt;
            for(int i = 0; i < nx - 1; ++i){
                u[i] = sin(2 * pi * x[i]) * (cos(2 * pi * tk) + sin(2 * pi * tk));
            }
            output(file, tk, u);
        }
        file.close();
    }
};

int main(){
    Parameters mypar;
    mypar.dt = 0.1;
    mypar.dx = 0.1;
    mypar.t_nt = 10;
    mypar.x_nx = 1;
    mypar.y_ny = 1;
    mypar.u_t0 = 100;
    mypar.u_nx = 20;
    heat1D myheat1d(0.1);
    heat2D myheat2d(0.1);
    wave1D mywave1d(1);
    ForwardEuler myFE(mypar);
    BackwardEuler myBE(mypar);
    Richardson myRI(mypar);
    DufortFrankel myDF(mypar);
    CrankNicolson myCN(mypar);
    ADI myADI(mypar);
    ExplicitFiniteDifference myEFD(mypar);
    ImplicitFiniteDifference myIFD(mypar);
    Exact myEX(mypar);
    myFE.solve(myheat1d);
    myBE.solve(myheat1d);
    myRI.solve(myheat1d);
    myDF.solve(myheat1d);
    myCN.solve(myheat1d);
    myFE.solve(myheat2d);
    myBE.solve(myheat2d);
    myCN.solve(myheat2d);
    myADI.solve(myheat2d);
    myEFD.solve(mywave1d);
    myIFD.solve(mywave1d);
    myEX.solve(mywave1d);
    return 0;
}