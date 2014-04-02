TITLE Large conductance Ca2+ activated K+ channel mslo, with channel noise.

COMMENT

Parameters from Cox et al. (1987) J Gen Physiol 110:257-81 (patch 1).

Written by Sungho Hong, Okinawa Institute of Science and Technology, March 2009.
Channel noise added by Daniele Linaro, TNB lab @ UA, January 2012.

ENDCOMMENT

VERBATIM
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

/*#define DEBUG*/

#define N_STATES      10
#define N_OPEN_STATES  5

/* transition matrix */
#define cn_A cn_A_mslo_cn
double cn_A_mslo_cn[N_STATES*N_STATES];
#define cn_tmp_A cn_tmp_A_mslo_cn
double cn_tmp_A_mslo_cn[N_STATES*N_STATES];
/* vector with the steady state values of all states */
#define cn_pinf cn_pinf_mslo_cn
double cn_pinf_mslo_cn[N_STATES];
/* vector with the indeces of the open states */
#define cn_open_states cn_open_states_mslo_cn
int cn_open_states_mslo_cn[N_OPEN_STATES] = {5,6,7,8,9};
#define cn_z cn_z_mslo_cn
double cn_z_mslo_cn[N_OPEN_STATES][N_STATES-1];
#define cn_tau cn_tau_mslo_cn
double cn_tau_mslo_cn[N_STATES-1];
#define cn_ss cn_ss_mslo_cn
double cn_ss_mslo_cn[N_OPEN_STATES][N_STATES-1];
#define cn_mu cn_mu_mslo_cn
double cn_mu_mslo_cn[N_OPEN_STATES][N_STATES-1];
#define cn_noise cn_noise_mslo_cn
double cn_noise_mslo_cn[N_OPEN_STATES][N_STATES-1];
#define cn_Z cn_Z_mslo_cn
double cn_Z_mslo_cn;

void init_cn(void);
void fill_A(void);
void compute_noise(double v, double ca, int num_channels);
int compute_sigmas(int num_channels);
int compute_taus();

ENDVERBATIM

NEURON {
  SUFFIX mslo_cn
  USEION k READ ek WRITE ik
  USEION ca READ cai
  RANGE g, gbar, ik, wca
  RANGE Nk, gamma_k, seed
}

UNITS { 
    (mV) = (millivolt)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (kilocoulombs)
    R = (k-mole) (joule/degC)
    (um) = (micrometer)
}

CONSTANT {
    q10 = 3
}

PARAMETER {
    gbar = 0.01 (S/cm2)
    wca = 1.0     : diffusion
    
    Qo = 0.73
    Qc = -0.67
    
    k1 = 1.0e3 (/mM)
    onoffrate = 1 (/ms)
    
    L0 = 1806
    Kc = 11.0e-3 (mM)
    Ko = 1.1e-3 (mM)
    
    pf0 = 2.39e-3  (/ms)
    pf1 = 7.0e-3  (/ms)
    pf2 = 40e-3   (/ms)
    pf3 = 295e-3  (/ms)
    pf4 = 557e-3  (/ms)
    
    pb0 = 3936e-3 (/ms)
    pb1 = 1152e-3 (/ms)
    pb2 = 659e-3  (/ms)
    pb3 = 486e-3  (/ms)
    pb4 = 92e-3  (/ms)
    
    gamma_k = 10  (pS)
    seed = 5061983 (1)
}

ASSIGNED {
    : rates
    c01    (/ms)
    c12    (/ms)
    c23    (/ms)
    c34    (/ms)
    o01    (/ms)
    o12    (/ms)
    o23    (/ms)
    o34    (/ms)
    f0     (/ms)
    f1     (/ms)
    f2     (/ms)
    f3     (/ms)
    f4     (/ms)

    c10    (/ms)
    c21    (/ms)
    c32    (/ms)
    c43    (/ms)
    o10    (/ms)
    o21    (/ms)
    o32    (/ms)
    o43    (/ms)
    b0     (/ms)
    b1     (/ms)
    b2     (/ms)
    b3     (/ms)
    b4     (/ms)
    
    v            (mV)
    cai          (mM)
    ek           (mV)
    ik           (milliamp/cm2)
    g            (S/cm2)
    celsius      (degC)
    Nk           (1)
    area         (um2)
}

STATE {
    C0 FROM 0 TO 1
    C1 FROM 0 TO 1
    C2 FROM 0 TO 1
    C3 FROM 0 TO 1
    C4 FROM 0 TO 1
    O0 FROM 0 TO 1
    O1 FROM 0 TO 1
    O2 FROM 0 TO 1
    O3 FROM 0 TO 1
    O4 FROM 0 TO 1
}

BREAKPOINT {
    SOLVE activation METHOD sparse
    VERBATIM
    compute_noise(v,wca*cai,Nk);
    g = gbar * (O0+O1+O2+O3+O4 + cn_Z);
    if (g < 0.) {
	g = 0.;
    }
    else if (g > gbar) {
	g = gbar;
    }
    ENDVERBATIM
    ik = g * (v - ek)
}

INITIAL {
    SOLVE activation STEADYSTATE sparse
    Nk = ceil(((1e-8)*area)*(gbar)/((1e-12)*gamma_k))
    printf("Nk = %.0f.\n", Nk)
    printf("Area = %f.\n", area)
    set_seed(seed)
    VERBATIM
    init_cn();
    ENDVERBATIM

}

KINETIC activation {
    rates(v, cai*wca)
    ~ C0 <-> C1      (c01,c10)
    ~ C1 <-> C2      (c12,c21)
    ~ C2 <-> C3      (c23,c32)
    ~ C3 <-> C4      (c34,c43)
    ~ O0 <-> O1      (o01,o10)
    ~ O1 <-> O2      (o12,o21)
    ~ O2 <-> O3      (o23,o32)
    ~ O3 <-> O4      (o34,o43)
    ~ C0 <-> O0      (f0 , b0)
    ~ C1 <-> O1      (f1 , b1)
    ~ C2 <-> O2      (f2 , b2)
    ~ C3 <-> O3      (f3 , b3)
    ~ C4 <-> O4      (f4 , b4)

CONSERVE C0 + C1 + C2 + C3 + C4 + O0 + O1 + O2 + O3 + O4 = 1
}

PROCEDURE rates(v(mV), ca (mM)) { 
    LOCAL qt, alpha, beta
    
    qt = q10^((celsius-23 (degC))/10 (degC))
    
    c01 = 4*ca*k1*onoffrate*qt
    c12 = 3*ca*k1*onoffrate*qt
    c23 = 2*ca*k1*onoffrate*qt
    c34 = 1*ca*k1*onoffrate*qt
    o01 = 4*ca*k1*onoffrate*qt
    o12 = 3*ca*k1*onoffrate*qt
    o23 = 2*ca*k1*onoffrate*qt
    o34 = 1*ca*k1*onoffrate*qt
    
    c10 = 1*Kc*k1*onoffrate*qt
    c21 = 2*Kc*k1*onoffrate*qt
    c32 = 3*Kc*k1*onoffrate*qt
    c43 = 4*Kc*k1*onoffrate*qt
    o10 = 1*Ko*k1*onoffrate*qt
    o21 = 2*Ko*k1*onoffrate*qt
    o32 = 3*Ko*k1*onoffrate*qt
    o43 = 4*Ko*k1*onoffrate*qt
    
    alpha = exp(Qo*FARADAY*v/R/(273.15 + celsius))
    beta  = exp(Qc*FARADAY*v/R/(273.15 + celsius))
    
    f0  = pf0*alpha*qt
    f1  = pf1*alpha*qt
    f2  = pf2*alpha*qt
    f3  = pf3*alpha*qt
    f4  = pf4*alpha*qt
    
    b0  = pb0*beta*qt
    b1  = pb1*beta*qt
    b2  = pb2*beta*qt
    b3  = pb3*beta*qt
    b4  = pb4*beta*qt
}

VERBATIM

void init_cn(void) {
    int i, j;
    cn_Z = 0.0;
    for (i=0; i<N_OPEN_STATES; i++) {
	for (j=0; j<N_STATES; j++) {
	    cn_z[i][j] = 0.0;
	}
    }
}

void fill_A(void) {
    int i, j;
    
    cn_A[ 0] = -(c01 + f0);
    cn_A[ 1] = c10;
    cn_A[ 2] = 0;
    cn_A[ 3] = 0;
    cn_A[ 4] = 0;
    cn_A[ 5] = b0;
    cn_A[ 6] = 0;
    cn_A[ 7] = 0;
    cn_A[ 8] = 0;
    cn_A[ 9] = 0;
    cn_A[10] = c01;
    cn_A[11] = -(c10+f1+c12);
    cn_A[12] = c21;
    cn_A[13] = 0;
    cn_A[14] = 0;
    cn_A[15] = 0;
    cn_A[16] = b1;
    cn_A[17] = 0;
    cn_A[18] = 0;
    cn_A[19] = 0;
    cn_A[20] = 0;
    cn_A[21] = c12;
    cn_A[22] = -(c21+f2+c23);
    cn_A[23] = c32;
    cn_A[24] = 0;
    cn_A[25] = 0;
    cn_A[26] = 0;
    cn_A[27] = b2;
    cn_A[28] = 0;
    cn_A[29] = 0;
    cn_A[30] = 0;
    cn_A[31] = 0;
    cn_A[32] = c23;
    cn_A[33] = -(c32+f3+c34);
    cn_A[34] = c43;
    cn_A[35] = 0;
    cn_A[36] = 0;
    cn_A[37] = 0;
    cn_A[38] = b3;
    cn_A[39] = 0;
    cn_A[40] = 0;
    cn_A[41] = 0;
    cn_A[42] = 0;
    cn_A[43] = c34;
    cn_A[44] = -(c43+f4);
    cn_A[45] = 0;
    cn_A[46] = 0;
    cn_A[47] = 0;
    cn_A[48] = 0;
    cn_A[49] = b4;
    cn_A[50] = f0;
    cn_A[51] = 0;
    cn_A[52] = 0;
    cn_A[53] = 0;
    cn_A[54] = 0;
    cn_A[55] = -(b0+o01);
    cn_A[56] = o10;
    cn_A[57] = 0;
    cn_A[58] = 0;
    cn_A[59] = 0;
    cn_A[60] = 0;
    cn_A[61] = f1;
    cn_A[62] = 0;
    cn_A[63] = 0;
    cn_A[64] = 0;
    cn_A[65] = o01;
    cn_A[66] = -(o10+b1+o12);
    cn_A[67] = o21;
    cn_A[68] = 0;
    cn_A[69] = 0;
    cn_A[70] = 0;
    cn_A[71] = 0;
    cn_A[72] = f2;
    cn_A[73] = 0;
    cn_A[74] = 0;
    cn_A[75] = 0;
    cn_A[76] = o12;
    cn_A[77] = -(o21+b2+o23);
    cn_A[78] = o32;
    cn_A[79] = 0;
    cn_A[80] = 0;
    cn_A[81] = 0;
    cn_A[82] = 0;
    cn_A[83] = f3;
    cn_A[84] = 0;
    cn_A[85] = 0;
    cn_A[86] = 0;
    cn_A[87] = o23;
    cn_A[88] = -(o32+b3+o34);
    cn_A[89] = o43;
    cn_A[90] = 0;
    cn_A[91] = 0;
    cn_A[92] = 0;
    cn_A[93] = 0;
    cn_A[94] = f4;
    cn_A[95] = 0;
    cn_A[96] = 0;
    cn_A[97] = 0;
    cn_A[98] = o34;
    cn_A[99] = -(o43+b4);
}
	
int compute_sigmas(int num_channels) {
    int i, j, k, flag, s;
    double pinf, coeff = 1.0 / num_channels;
    double c[N_STATES];
    gsl_permutation *p;
    gsl_matrix_view A_view;
    gsl_vector_view c_view, pinf_view;
    
    for (i=0; i<N_STATES; i++) {
	c[i] = 0.0;
    }
    c[cn_open_states[0]] = 1.0;
    
    for (i=0; i<N_STATES*N_STATES; i++) {
	cn_tmp_A[i] = cn_A[i];
    }
    for (i=cn_open_states[0]*N_STATES; i<(cn_open_states[0]+1)*N_STATES; i++) {
	cn_tmp_A[i] = 1.0;
    }
    
    A_view = gsl_matrix_view_array(cn_tmp_A, N_STATES, N_STATES);
    c_view = gsl_vector_view_array(c, N_STATES);
    pinf_view = gsl_vector_view_array(cn_pinf, N_STATES);
    
    p = gsl_permutation_alloc(N_STATES);
    flag = gsl_linalg_LU_decomp(&A_view.matrix, p, &s);
    if (flag != 0) {
	goto sigmas_end;
    }
    
    flag = gsl_linalg_LU_solve(&A_view.matrix, p, &c_view.vector, &pinf_view.vector);
    if (flag != 0) {
	goto sigmas_end;
    }
    
    for (i=0; i<N_OPEN_STATES; i++) {
	pinf = gsl_vector_get(&pinf_view.vector, cn_open_states[i]);
	for (j=0, k=0; j<N_STATES; j++) {
	    if (j != cn_open_states[i]) {
		cn_ss[i][k] = coeff * gsl_vector_get(&pinf_view.vector,j) * pinf;
		k++;
	    }
	}
    }
    
sigmas_end:
    
    gsl_permutation_free(p);    
    return flag;
}

int compute_taus() {
    int i, j, flag;
    double lambda;
    double eigs[N_STATES*2]; /* in general, eigenvalues of a non-symmetric matrix are complex */
    
    gsl_matrix_view A_view;
    gsl_vector_complex_view eigs_view;
    gsl_eigen_nonsymm_workspace *w;
    
    for (i=0; i<N_STATES*N_STATES; i++) {
	cn_tmp_A[i] = cn_A[i];
    }
    
    A_view = gsl_matrix_view_array(cn_tmp_A, N_STATES, N_STATES);
    eigs_view = gsl_vector_complex_view_array(eigs, N_STATES);
    w = gsl_eigen_nonsymm_alloc(N_STATES);
    gsl_eigen_nonsymm_params(0, 1, w);
    flag = gsl_eigen_nonsymm(&A_view.matrix, &eigs_view.vector, w);
    if (flag != 0) {
	goto taus_end;
    }
    
    for (i=0, j=0; i<N_STATES; i++) {
        lambda = GSL_REAL(gsl_vector_complex_get(&eigs_view.vector,i));
        if (lambda < -1e-8) {
            cn_tau[j] = - 1.0 / lambda;
	    j++;
	}
    }
    
taus_end:
    
    gsl_eigen_nonsymm_free(w);
    return flag;
}

void compute_noise(double v, double ca, int num_channels) {
    int i, j;
    
    fill_A();
    
    if (compute_taus() != 0) {
	printf("compute_taus returned an error flag @ t = %g: V = %g.\n", t, v);
	exit(1);
    }
    
    if (compute_sigmas(num_channels) != 0) {
	printf("compute_sigmas returned an error flag @ t = %g: V = %g.\n", t, v);
	exit(1);
    }
    
    cn_Z = 0.0;
    
    for (i=0; i<N_OPEN_STATES; i++) {
	for (j=0; j<N_STATES-1; j++) {
            cn_mu[i][j] = exp(-dt/cn_tau[j]);
            cn_noise[i][j] = sqrt(cn_ss[i][j]*(1-cn_mu[i][j]*cn_mu[i][j])) * normrand(0,1);
            cn_z[i][j] = cn_z[i][j]*cn_mu[i][j] + cn_noise[i][j];
            cn_Z += cn_z[i][j];
	}
    }
}

ENDVERBATIM
