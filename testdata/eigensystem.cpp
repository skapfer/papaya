

#include "../util.h"
#include <iostream>


void rot (double phi, mat_t *m) {
    double c = cos (phi);
    double s = sin (phi);
    double a1 = (*m)(0,0);
    double b2 = (*m)(1,1);
    (*m)(0,0) = a1*c*c + b2*s*s;
    (*m)(1,1) = b2*c*c + a1*s*s;
    (*m)(0,1) = (b2-a1)*c*s;
    (*m)(1,0) = (*m)(0,1);
}

vec_t ev [] = {
    vec_t (0., 1.),
    vec_t (.1, 1.),
    vec_t (.5, 1.),
    vec_t (.99, 1.),
    vec_t (1-1e-6, 1.),
    vec_t (1-1e-8, 1.),
    vec_t (-0., 1.),
    vec_t (-.1, 1.),
    vec_t (-.5, 1.),
    vec_t (-.99, 1.),
    vec_t (-1+1e-6, 1.),
    vec_t (-1+1e-8, 1.),
};

int num_ev = sizeof (ev) / sizeof (*ev);

int main () {
    mat_t mat;
    for (int i = 0; i != num_ev; ++i) {
        for (int n = 0; n != 120; ++n) {
            double phi = -M_PI/2 + M_PI/120*n;
            mat(0,0) = ev[i][0];
            mat(1,1) = ev[i][1];
            rot (-phi, &mat);
            EigenSystem es;
            eigensystem_symm (&es, mat);
            if (fabs (es.eval[0] - ev[i][0]) > fabs (es.eval[1] - ev[i][0]))
                swap_eigenvalues (&es);
            double del_ev = std::max (
                fabs (es.eval[0]-ev[i][0]),
                fabs (es.eval[1]-ev[i][1]));
            double retphi = atan2 (es.evec[0][1], es.evec[0][0]);
            if (retphi > M_PI/2)
                retphi -= M_PI;
            if (retphi < -M_PI/2)
                retphi += M_PI;
            fprintf (stdout, "%i %.20e %.20e\n", 200*i + n,
                retphi-phi, del_ev);
        }
    }
    mat(0,0) = -2.;
    mat(0,1) = 0.;
    mat(1,0) = 0.;
    mat(1,1) = -3.;
    EigenSystem es;
    eigensystem_symm (&es, mat);
    //es.dump (std::cout);
    return 0;
}
