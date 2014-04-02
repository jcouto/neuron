TITLE Rsg sodium channel with channel noise
: Resurgent sodium channel (with blocking particle)
: with updated kinetic parameters from Raman and Bean  

VERBATIM
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

/*#define DEBUG*/

#define ALPHA(v)   alpha*exp((v)/x1)
#define BETA(v)    beta*exp((v)/x2)
#define GAMMA(v)   gamma*exp((v)/x3)
#define DELTA(v)   delta*exp((v)/x4)
#define EPSILON(v) epsilon*exp((v)/x5)
#define ZETA(v)    zeta*exp((v)/x6)
#define NSTATES    13
#define OPEN_STATE 5

/* transition matrix, it will be factorized in A = U D S, where S = U^-1 */
#define cn_A cn_A_naRsgcn
double cn_A_naRsgcn[NSTATES*NSTATES];
#define cn_tmp_A cn_tmp_A_naRsgcn
double cn_tmp_A_naRsgcn[NSTATES*NSTATES];
/* matrix with the eigenvectors of A */
#define cn_U cn_U_naRsgcn 
double cn_U_naRsgcn[NSTATES*NSTATES*2];
/* matrix with the inverse of U */
#define cn_S cn_S_naRsgcn 
double cn_S_naRsgcn[NSTATES*NSTATES];
/* diagonal matrix with the eigenvalues */
#define cn_D cn_D_naRsgcn 
double cn_D_naRsgcn[NSTATES*NSTATES];
/* matrix with the LU decomposition of U */
#define cn_LU cn_LU_naRsgcn 
double cn_LU_naRsgcn[NSTATES*NSTATES];
/* vector with the steady state values of all states */
#define cn_pinf cn_pinf_naRsgcn
double cn_pinf_naRsgcn[NSTATES];
/* vector with the ``open'' state set to 1 */
#define cn_c cn_c_naRsgcn
double cn_c_naRsgcn[NSTATES] = {0,0,0,0,0,1,0,0,0,0,0,0,0};

#define cn_z cn_z_naRsgcn
double cn_z_naRsgcn[NSTATES-1];
#define cn_tau cn_tau_naRsgcn
double cn_tau_naRsgcn[NSTATES-1];
#define cn_ss cn_ss_naRsgcn
double cn_ss_naRsgcn[NSTATES-1];
#define cn_mu cn_mu_naRsgcn
double cn_mu_naRsgcn[NSTATES-1];
#define cn_noise cn_noise_naRsgcn
double cn_noise_naRsgcn[NSTATES-1];
#define cn_Z cn_Z_naRsgcn
double cn_Z_naRsgcn;

void init_cn(void);
void test_cn(void);
void fill_A(double v);
int compute_sigmas(double v, int num_channels);
int compute_taus(double v);
int compute_coefficients(double v, int open_state, int num_channels);
void compute_noise(double v);

ENDVERBATIM

NEURON {
    SUFFIX naRsg_cn
    USEION na READ ena WRITE ina
    RANGE g, gbar
    RANGE Nna, gamma_na, seed
}

UNITS { 
    (mV) = (millivolt)
    (S) = (siemens)
    (pS) = (picosiemens)
    (um) = (micrometer)
}

PARAMETER {
    gbar = .015			(S/cm2)
    
    : kinetic parameters
    Con = 0.005			(/ms)		: closed -> inactivated transitions
    Coff = 0.5			(/ms)		: inactivated -> closed transitions
    Oon = .75			(/ms)		: open -> Ineg transition
    Ooff = 0.005		(/ms)		: Ineg -> open transition
    alpha = 150			(/ms)		: activation
    beta = 3			(/ms)		: deactivation
    gamma = 150			(/ms)		: opening
    delta = 40			(/ms)		: closing, greater than BEAN/KUO = 0.2
    epsilon = 1.75		(/ms)		: open -> Iplus for tau = 0.3 ms at +30 with x5
    zeta = 0.03			(/ms)		: Iplus -> open for tau = 25 ms at -30 with x6
    
    : Vdep
    x1 = 20			(mV)		: Vdep of activation (alpha)
    x2 = -20			(mV)		: Vdep of deactivation (beta)
    x3 = 1e12			(mV)		: Vdep of opening (gamma)
    x4 = -1e12			(mV)		: Vdep of closing (delta)
    x5 = 1e12			(mV)		: Vdep into Ipos (epsilon)
    x6 = -25			(mV)		: Vdep out of Ipos (zeta)
    
    gamma_na = 10               (pS)
    seed = 5061983              (1)
}

ASSIGNED {
    alfac   				: microscopic reversibility factors
    btfac				
    
    : rates
    f01  		(/ms)
    f02  		(/ms)
    f03 		(/ms)
    f04			(/ms)
    f0O 		(/ms)
    fip 		(/ms)
    f11 		(/ms)
    f12 		(/ms)
    f13 		(/ms)
    f14 		(/ms)
    f1n 		(/ms)
    fi1 		(/ms)
    fi2 		(/ms)
    fi3 		(/ms)
    fi4 		(/ms)
    fi5 		(/ms)
    fin 		(/ms)
    
    b01 		(/ms)
    b02 		(/ms)
    b03 		(/ms)
    b04  		(/ms)
    b0O 		(/ms)
    bip 		(/ms)
    b11  		(/ms)
    b12 		(/ms)
    b13 		(/ms)
    b14 		(/ms)
    b1n 		(/ms)
    bi1 		(/ms)
    bi2 		(/ms)
    bi3 		(/ms)
    bi4 		(/ms)
    bi5 		(/ms)
    bin 		(/ms)
    
    v			(mV)
    ena			(mV)
    ina 		(milliamp/cm2)
    g			(S/cm2)
    
    Nna                 (1)
    area                (um2)
}

STATE {
    C1 FROM 0 TO 1
    C2 FROM 0 TO 1
    C3 FROM 0 TO 1
    C4 FROM 0 TO 1
    C5 FROM 0 TO 1
    I1 FROM 0 TO 1
    I2 FROM 0 TO 1
    I3 FROM 0 TO 1
    I4 FROM 0 TO 1
    I5 FROM 0 TO 1
    O FROM 0 TO 1
    B FROM 0 TO 1
    I6 FROM 0 TO 1
}

BREAKPOINT {
    SOLVE activation METHOD sparse
    VERBATIM
    compute_noise(v);
    g = gbar * (O + cn_Z);
    if (g < 0.) {
	g = 0.;
    }
    else if (g > gbar) {
	g = gbar;
    }
    ENDVERBATIM
    ina = g * (v - ena)
}

INITIAL {
    rates(v)
    SOLVE seqinitial
    Nna = ceil(((1e-8)*area)*(gbar)/((1e-12)*gamma_na))
    printf("Nna = %.0f.\n", Nna)
    printf("Area = %f.\n", area)
    set_seed(seed)
    VERBATIM
    init_cn();
    test_cn();
    ENDVERBATIM
}

KINETIC activation {
    rates(v)
    ~ C1 <-> C2					(f01,b01)
    ~ C2 <-> C3					(f02,b02)
    ~ C3 <-> C4					(f03,b03)
    ~ C4 <-> C5					(f04,b04)
    ~ C5 <-> O					(f0O,b0O)
    ~ O <-> B					(fip,bip)
    ~ O <-> I6					(fin,bin)
    ~ I1 <-> I2					(f11,b11)
    ~ I2 <-> I3					(f12,b12)
    ~ I3 <-> I4					(f13,b13)
    ~ I4 <-> I5					(f14,b14)
    ~ I5 <-> I6					(f1n,b1n)
    ~ C1 <-> I1					(fi1,bi1)
    ~ C2 <-> I2					(fi2,bi2)
    ~ C3 <-> I3					(fi3,bi3)
    ~ C4 <-> I4					(fi4,bi4)
    ~ C5 <-> I5					(fi5,bi5)
    
    CONSERVE C1 + C2 + C3 + C4 + C5 + O + B + I1 + I2 + I3 + I4 + I5 + I6 = 1
}

LINEAR seqinitial { : sets initial equilibrium
    ~          I1*bi1 + C2*b01 - C1*(    fi1+f01) = 0
    ~ C1*f01 + I2*bi2 + C3*b02 - C2*(b01+fi2+f02) = 0
    ~ C2*f02 + I3*bi3 + C4*b03 - C3*(b02+fi3+f03) = 0
    ~ C3*f03 + I4*bi4 + C5*b04 - C4*(b03+fi4+f04) = 0
    ~ C4*f04 + I5*bi5 + O*b0O - C5*(b04+fi5+f0O) = 0
    ~ C5*f0O + B*bip + I6*bin - O*(b0O+fip+fin) = 0
    ~ O*fip + B*bip = 0
    
    ~          C1*fi1 + I2*b11 - I1*(    bi1+f11) = 0
    ~ I1*f11 + C2*fi2 + I3*b12 - I2*(b11+bi2+f12) = 0
    ~ I2*f12 + C3*fi3 + I4*bi3 - I3*(b12+bi3+f13) = 0
    ~ I3*f13 + C4*fi4 + I5*b14 - I4*(b13+bi4+f14) = 0
    ~ I4*f14 + C5*fi5 + I6*b1n - I5*(b14+bi5+f1n) = 0
    
    ~ C1 + C2 + C3 + C4 + C5 + O + B + I1 + I2 + I3 + I4 + I5 + I6 = 1
}

PROCEDURE rates(v (mV)) {
    alfac = (Oon/Con)^(1/4)
    btfac = (Ooff/Coff)^(1/4) 
    f01 = 4 * alpha * exp(v/x1)
    f02 = 3 * alpha * exp(v/x1)
    f03 = 2 * alpha * exp(v/x1)
    f04 = 1 * alpha * exp(v/x1)
    f0O = gamma * exp(v/x3)
    fip = epsilon * exp(v/x5)
    f11 = 4 * alpha * alfac * exp(v/x1)
    f12 = 3 * alpha * alfac * exp(v/x1)
    f13 = 2 * alpha * alfac * exp(v/x1)
    f14 = 1 * alpha * alfac * exp(v/x1)
    f1n = gamma * exp(v/x3)
    fi1 = Con
    fi2 = Con * alfac
    fi3 = Con * alfac^2
    fi4 = Con * alfac^3
    fi5 = Con * alfac^4
    fin = Oon
    
    b01 = 1 * beta * exp(v/x2)
    b02 = 2 * beta * exp(v/x2)
    b03 = 3 * beta * exp(v/x2)
    b04 = 4 * beta * exp(v/x2)
    b0O = delta * exp(v/x4)
    bip = zeta * exp(v/x6)
    b11 = 1 * beta * btfac * exp(v/x2)
    b12 = 2 * beta * btfac * exp(v/x2)
    b13 = 3 * beta * btfac * exp(v/x2)
    b14 = 4 * beta * btfac * exp(v/x2)
    b1n = delta * exp(v/x4)
    bi1 = Coff
    bi2 = Coff * btfac
    bi3 = Coff * btfac^2
    bi4 = Coff * btfac^3
    bi5 = Coff * btfac^4
    bin = Ooff
}

VERBATIM

void init_cn(void) {
    int i;
    cn_Z = 0.0;
    for (i=0; i<NSTATES-1; i++) {
	cn_z[i] = 0.0;
    }
}

void test_cn(void) {
    int i, j;
    FILE *fid;
    fid = fopen("cn_coeff.dat","w");
    fprintf(stdout, "Computing coefficients... ");
    fflush(stdout);
    for (i=-100; i<=50; i++) {
	fprintf(fid, "%e", (double) i);
	fill_A((double) i);
	/*compute_coefficients((double) i, OPEN_STATE, 1);*/
	compute_sigmas((double) i, 1);
	for (j=0; j<NSTATES-1; j++) {
	    fprintf(fid, " %e", cn_ss[j]);
	}
	compute_taus((double) i);
	for (j=0; j<NSTATES-1; j++) {
	    fprintf(fid, " %e", cn_tau[j]);
	}
	fprintf(fid, "\n");
    }
    fprintf(stdout, "done!\n");
    fclose(fid);
}

void fill_A(double v) {
    cn_A[0] = -(4*ALPHA(v) + Con);
    cn_A[1] = BETA(v);
    cn_A[2] = 0;
    cn_A[3] = 0;
    cn_A[4] = 0;
    cn_A[5] = 0;
    cn_A[6] = 0;
    cn_A[7] = Coff;
    cn_A[8] = 0;
    cn_A[9] = 0;
    cn_A[10] = 0;
    cn_A[11] = 0;
    cn_A[12] = 0;
    cn_A[13] = 4*ALPHA(v);
    cn_A[14] = -(BETA(v) + Con*alfac + 3*ALPHA(v));
    cn_A[15] = 2*BETA(v);
    cn_A[16] = 0;
    cn_A[17] = 0;
    cn_A[18] = 0;
    cn_A[19] = 0;
    cn_A[20] = 0;
    cn_A[21] = Coff*btfac;
    cn_A[22] = 0;
    cn_A[23] = 0;
    cn_A[24] = 0;
    cn_A[25] = 0;
    cn_A[26] = 0;
    cn_A[27] = 3*ALPHA(v);
    cn_A[28] = -(2*BETA(v) + Con*alfac*alfac + 2*ALPHA(v));
    cn_A[29] = 3*BETA(v);
    cn_A[30] = 0;
    cn_A[31] = 0;
    cn_A[32] = 0;
    cn_A[33] = 0;
    cn_A[34] = 0;
    cn_A[35] = Coff*btfac*btfac;
    cn_A[36] = 0;
    cn_A[37] = 0;
    cn_A[38] = 0;
    cn_A[39] = 0;
    cn_A[40] = 0;
    cn_A[41] = 2*ALPHA(v);
    cn_A[42] = -(3*BETA(v) + Con*alfac*alfac*alfac + ALPHA(v));
    cn_A[43] = 4*BETA(v);
    cn_A[44] = 0;
    cn_A[45] = 0;
    cn_A[46] = 0;
    cn_A[47] = 0;
    cn_A[48] = 0;
    cn_A[49] = Coff*btfac*btfac*btfac;
    cn_A[50] = 0;
    cn_A[51] = 0;
    cn_A[52] = 0;
    cn_A[53] = 0;
    cn_A[54] = 0;
    cn_A[55] = ALPHA(v);
    cn_A[56] = -(4*BETA(v) + Con*alfac*alfac*alfac*alfac + GAMMA(v));
    cn_A[57] = DELTA(v);
    cn_A[58] = 0;
    cn_A[59] = 0;
    cn_A[60] = 0;
    cn_A[61] = 0;
    cn_A[62] = 0;
    cn_A[63] = Coff*btfac*btfac*btfac*btfac;
    cn_A[64] = 0;
    cn_A[65] = 0;
    cn_A[66] = 0;
    cn_A[67] = 0;
    cn_A[68] = 0;
    cn_A[69] = GAMMA(v);
    cn_A[70] = -(DELTA(v) + Oon + EPSILON(v));
    cn_A[71] = ZETA(v);
    cn_A[72] = 0;
    cn_A[73] = 0;
    cn_A[74] = 0;
    cn_A[75] = 0;
    cn_A[76] = 0;
    cn_A[77] = Ooff;
    cn_A[78] = 0;
    cn_A[79] = 0;
    cn_A[80] = 0;
    cn_A[81] = 0;
    cn_A[82] = 0;
    cn_A[83] = EPSILON(v);
    cn_A[84] = -ZETA(v);
    cn_A[85] = 0;
    cn_A[86] = 0;
    cn_A[87] = 0;
    cn_A[88] = 0;
    cn_A[89] = 0;
    cn_A[90] = 0;
    cn_A[91] = Con;
    cn_A[92] = 0;
    cn_A[93] = 0;
    cn_A[94] = 0;
    cn_A[95] = 0;
    cn_A[96] = 0;
    cn_A[97] = 0;
    cn_A[98] = -(Coff + 4*ALPHA(v)*alfac);
    cn_A[99] = BETA(v)*btfac;
    cn_A[100] = 0;
    cn_A[101] = 0;
    cn_A[102] = 0;
    cn_A[103] = 0;
    cn_A[104] = 0;
    cn_A[105] = Con*alfac;
    cn_A[106] = 0;
    cn_A[107] = 0;
    cn_A[108] = 0;
    cn_A[109] = 0;
    cn_A[110] = 0;
    cn_A[111] = 4*ALPHA(v)*alfac;
    cn_A[112] = -(BETA(v)*btfac + Coff*btfac + 3*ALPHA(v)*alfac);
    cn_A[113] = 2*BETA(v)*btfac;
    cn_A[114] = 0;
    cn_A[115] = 0;
    cn_A[116] = 0;
    cn_A[117] = 0;
    cn_A[118] = 0;
    cn_A[119] = Con*alfac*alfac;
    cn_A[120] = 0;
    cn_A[121] = 0;
    cn_A[122] = 0;
    cn_A[123] = 0;
    cn_A[124] = 0;
    cn_A[125] = 3*ALPHA(v)*alfac;
    cn_A[126] = -(2*BETA(v)*btfac + Coff*btfac*btfac + 2*ALPHA(v)*alfac);
    cn_A[127] = 3*BETA(v)*btfac;
    cn_A[128] = 0;
    cn_A[129] = 0;
    cn_A[130] = 0;
    cn_A[131] = 0;
    cn_A[132] = 0;
    cn_A[133] = Con*alfac*alfac*alfac;
    cn_A[134] = 0;
    cn_A[135] = 0;
    cn_A[136] = 0;
    cn_A[137] = 0;
    cn_A[138] = 0;
    cn_A[139] = 2*ALPHA(v)*alfac;
    cn_A[140] = -(3*BETA(v)*btfac + Coff*btfac*btfac*btfac + ALPHA(v)*alfac);
    cn_A[141] = 4*BETA(v)*btfac;
    cn_A[142] = 0;
    cn_A[143] = 0;
    cn_A[144] = 0;
    cn_A[145] = 0;
    cn_A[146] = 0;
    cn_A[147] = Con*alfac*alfac*alfac*alfac;
    cn_A[148] = 0;
    cn_A[149] = 0;
    cn_A[150] = 0;
    cn_A[151] = 0;
    cn_A[152] = 0;
    cn_A[153] = ALPHA(v)*alfac;
    cn_A[154] = -(4*BETA(v)*btfac + Coff*btfac*btfac*btfac*btfac + GAMMA(v));
    cn_A[155] = DELTA(v);
    cn_A[156] = 0;
    cn_A[157] = 0;
    cn_A[158] = 0;
    cn_A[159] = 0;
    cn_A[160] = 0;
    cn_A[161] = Oon;
    cn_A[162] = 0;
    cn_A[163] = 0;
    cn_A[164] = 0;
    cn_A[165] = 0;
    cn_A[166] = 0;
    cn_A[167] = GAMMA(v);
    cn_A[168] = -(DELTA(v) + Ooff);
}

int compute_sigmas(double v, int num_channels) {
    int i, j, flag, s, open_state;
    double pinf, coeff = 1.0 / num_channels;
    gsl_permutation *p;
    gsl_matrix_view A_view;
    gsl_vector_view c_view, pinf_view;
    
    for (i=0; i<NSTATES; i++) {
	if (cn_c[i]) {
	    open_state = i;
	    break;
	}
    }
    
    for (i=0; i<NSTATES*NSTATES; i++) {
	cn_tmp_A[i] = cn_A[i];
    }
    for (i=open_state*NSTATES; i<(open_state+1)*NSTATES; i++) {
	cn_tmp_A[i] = 1.0;
    }
    
    A_view = gsl_matrix_view_array(cn_tmp_A, NSTATES, NSTATES);
    c_view = gsl_vector_view_array(cn_c, NSTATES);
    pinf_view = gsl_vector_view_array(cn_pinf, NSTATES);
    
    p = gsl_permutation_alloc(NSTATES);
    flag = gsl_linalg_LU_decomp(&A_view.matrix, p, &s);
    if (flag != 0) {
	goto sigmas_end;
    }
    
    flag = gsl_linalg_LU_solve(&A_view.matrix, p, &c_view.vector, &pinf_view.vector);
    if (flag != 0) {
	goto sigmas_end;
    }
    
    pinf = gsl_vector_get(&pinf_view.vector, open_state);
    for (i=0, j=0; i<NSTATES; i++) {
	if (i != open_state) {
	    cn_ss[j] = coeff * gsl_vector_get(&pinf_view.vector,i) * pinf;
	    j++;
	}
    }
    
sigmas_end:
    
    gsl_permutation_free(p);    
    return flag;
}

int compute_taus(double v) {
    int i, j, flag;
    double lambda;
    double eigs[NSTATES*2]; /* in general, eigenvalues of a non-symmetric matrix are complex */
    
    gsl_matrix_view A_view;
    gsl_vector_complex_view eigs_view;
    gsl_eigen_nonsymm_workspace *w;
    
    for (i=0; i<NSTATES*NSTATES; i++) {
	cn_tmp_A[i] = cn_A[i];
    }
    
    A_view = gsl_matrix_view_array(cn_tmp_A, NSTATES, NSTATES);
    eigs_view = gsl_vector_complex_view_array(eigs, NSTATES);
    w = gsl_eigen_nonsymm_alloc(NSTATES);
    gsl_eigen_nonsymm_params(0, 1, w);
    flag = gsl_eigen_nonsymm(&A_view.matrix, &eigs_view.vector, w);
    if (flag != 0) {
	goto taus_end;
    }
    
    for (i=0, j=0; i<NSTATES; i++) {
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

int compute_coefficients(double v, int open_state, int num_channels) {
    int i, j, flag, s;
    double k, lambda, pinf, coeff;
    double eigs[NSTATES*2]; /* in general, eigenvalues of a non-symmetric matrix are complex */
    
    gsl_matrix_view A_view, LU_view, S_view, D_view;
    gsl_matrix_complex_view U_view;
    
    gsl_vector_complex_view eigs_view;
    
    gsl_eigen_nonsymmv_workspace *w;
    
    gsl_permutation *p;
    
    A_view = gsl_matrix_view_array(cn_A, NSTATES, NSTATES);
    U_view = gsl_matrix_complex_view_array(cn_U, NSTATES, NSTATES);
    LU_view = gsl_matrix_view_array(cn_LU, NSTATES, NSTATES);
    S_view = gsl_matrix_view_array(cn_S, NSTATES, NSTATES);
    D_view = gsl_matrix_view_array(cn_D, NSTATES, NSTATES);
    
    eigs_view = gsl_vector_complex_view_array(eigs, NSTATES);
    
    p = gsl_permutation_alloc(NSTATES);
    
    w = gsl_eigen_nonsymmv_alloc(NSTATES);
    gsl_eigen_nonsymmv_params(1, w);
    
    flag = gsl_eigen_nonsymmv(&A_view.matrix, &eigs_view.vector, &U_view.matrix, w);
    
    if (flag != 0) {
	goto coeff_end;
    }
    
    gsl_eigen_nonsymmv_sort(&eigs_view.vector, &U_view.matrix, GSL_EIGEN_SORT_ABS_ASC);
    k = 0;
    for (i=0; i<NSTATES; i++) {
        k += GSL_REAL(gsl_matrix_complex_get(&U_view.matrix, i, 0));
        lambda = GSL_REAL(gsl_vector_complex_get(&eigs_view.vector,i));
        gsl_matrix_set(&D_view.matrix, i, i, lambda);
        if (i) {
            cn_tau[i-1] = - 1.0 / lambda;
	}
        for (j=0; j<NSTATES; j++) {
            gsl_matrix_set(&LU_view.matrix, i, j, GSL_REAL(gsl_matrix_complex_get(&U_view.matrix,i,j)));
        }
    }
    pinf = GSL_REAL(gsl_matrix_complex_get(&U_view.matrix, open_state, 0))/k;
    
    flag = gsl_linalg_LU_decomp(&LU_view.matrix, p, &s);
    if (flag != 0) {
	goto coeff_end;
    }
    
    flag = gsl_linalg_LU_invert(&LU_view.matrix, p, &S_view.matrix);
    if (flag != 0) {
        goto coeff_end;
    }
    
    coeff = 1.0 / num_channels;
    for (i=1; i<NSTATES; i++) {
	cn_ss[i-1] = coeff * pinf * 
        GSL_REAL(gsl_matrix_complex_get(&U_view.matrix, open_state, i)) *
        gsl_matrix_get(&S_view.matrix, i, open_state);
    }
    
    coeff_end:
    
    gsl_eigen_nonsymmv_free(w);
    gsl_permutation_free(p);
    
    return flag;
}

void compute_noise(double v) {
    int i;
    
    fill_A(v);
    
    /*
    if (compute_coefficients(v,OPEN_STATE,Nna) != 0) {
	printf("compute_coefficients returned an error flag @ t = %g: V = %g.\n", t, v);
	exit(1);
    }
    */
    
    if (compute_taus(v) != 0) {
	printf("compute_taus returned an error flag @ t = %g: V = %g.\n", t, v);
	exit(1);
    }

    if (compute_sigmas(v,Nna) != 0) {
	printf("compute_sigmas returned an error flag @ t = %g: V = %g.\n", t, v);
	exit(1);
    }

    #ifdef DEBUG
    printf("--\nV = %g\n  taus = [ ", v);
    for (i=0; i<NSTATES-1; i++) {
	printf("%g ", cn_tau[i]);
    }
    printf("]\nsigmas = [ ");
    for (i=0; i<NSTATES-1; i++) {
	printf("%g ", cn_ss[i]);
    }
    printf("]\n");
    getchar();
    #endif
    
    cn_Z = 0.0;
    for (i=0; i<NSTATES-1; i++) {
        cn_mu[i] = exp(-dt/cn_tau[i]);
        cn_noise[i] = sqrt(cn_ss[i]*(1-cn_mu[i]*cn_mu[i])) * normrand(0,1);
        cn_z[i] = cn_z[i]*cn_mu[i] + cn_noise[i];
        cn_Z += cn_z[i];
    }
}

ENDVERBATIM

