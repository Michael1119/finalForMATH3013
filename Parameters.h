#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters{
public:
    // step size
    double dt;
    double dx; // assume dy = dx
    // domain
    double t_nt;
    double x_nx;
    double y_ny;
    // initial condition
    double u_t0;
    // boundary condition
    double u_nx;
};

#endif // PARAMETERS_H