#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>

// https://www.gnu.org/software/gsl/manual/html_node/ODE-Example-programs.html#ODE-Example-programs
// http://stackoverflow.com/questions/27913858/how-to-consult-gsl-ode-with-many-parameters-and-harmonic-functions

typedef struct {
    double beta;
    double mu;
    double gamma;
} ode_params;

int 
SI_dt (double t, const double y[], double f[], void *params)
{
    ode_params *my_params = params;
    double beta = my_params->beta;
    double mu = my_params->mu;
    double gamma = my_params->gamma;
    double I = y[0];
    double S = y[1];
    double N = I + S;

    f[0] = beta*I*S/N - (gamma + mu)*I;
    f[1] = -beta*I*S/N + (gamma + mu)*I;
    return GSL_SUCCESS;
}

int 
SI_jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    ode_params *my_params = params;
    double beta = my_params->beta;
    double mu = my_params->mu;
    double gamma = my_params->gamma;
    double I = y[0];
    double S = y[1];
    double N = I + S;

    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix *m = &dfdy_mat.matrix;
    gsl_matrix_set(m, 0, 0, beta*S/N - gamma - mu);
    gsl_matrix_set(m, 0, 1, beta*I/N);
    gsl_matrix_set(m, 1, 0, -beta*S/N + gamma + mu);
    gsl_matrix_set(m, 1, 1, -beta*I/N);
    
    // no direct time-dependence
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}
    
int 
main (void)
{
    ode_params my_params = {0.011, 0.005, 0.005};
    gsl_odeiv2_system sys = {SI_dt, SI_jac, 2, &my_params};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys,
            gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
    int i;
    double t = 0.0, t1 = 100000.0;
    double y[2] = {1, 4999};

    for (i = 1; i <= 100; ++i) {
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

    gsl_odeiv2_driver_free(d);
    return 0;
}
