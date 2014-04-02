: HH Slow TEA-insensitive Purkinje potassium current
: Created 8/7/02 - nwg

NEURON {
    SUFFIX kpkjslow_cn
    USEION k READ ek WRITE ik
    RANGE gkbar, ik
    GLOBAL ninf, ntau
    : channel noise - start
    RANGE gk, gamma_k
    RANGE Nk, one_over_Nk
    RANGE seed
    : channel noise - end
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    : channel noise - start
    (S) = (siemens)
    (pS) = (picosiemens)
    : channel noise - end
}

PARAMETER {
    v		(mV)
    gkbar = .004	(mho/cm2)
    
    nivh = -16.5	(mV)
    nik = 18.4
    
    ek
    : channel noise - start
    seed = 5061983 (1)
    gamma_k = 10 (pS)
    : channel noise - end
}

ASSIGNED {
    ik
    ninf
    ntau		(ms)
    : channel noise - start

    gk (S/cm2)
    Nk (1)
    one_over_Nk (1)
    
    dt (ms)
    area (um2)
    
    tau1_kpkjslow (ms) tau2_kpkjslow (ms) tau3_kpkjslow (ms) tau4_kpkjslow
    sigma1_kpkjslow (ms2) sigma2_kpkjslow (ms2) sigma3_kpkjslow (ms2) sigma4_kpkjslow (ms2)
    noise1_kpkjslow noise2_kpkjslow noise3_kpkjslow noise4_kpkjslow
    mu1_kpkjslow mu2_kpkjslow mu3_kpkjslow mu4_kpkjslow
    
    : channel noise - end
}

STATE {
    n
    : channel noise - start
    z1_kpkjslow z2_kpkjslow z3_kpkjslow z4_kpkjslow
    : channel noise - end
}

INITIAL {
    rates(v)
    n = ninf
    : channel noise - start
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))
    one_over_Nk = 1.0 / Nk
    printf("kpkjslow>> the number of channels is %.0f.\n", Nk)
    z1_kpkjslow = 0.
    z2_kpkjslow = 0.
    z3_kpkjslow = 0.
    z4_kpkjslow = 0.
    : channel noise - end
}

BREAKPOINT {
    SOLVE states
    gk = gkbar * (n*n*n*n + z1_kpkjslow+z2_kpkjslow+z3_kpkjslow+z4_kpkjslow)
    if (gk < 0) {
        gk = 0
    }
    else if (gk > gkbar) {
        gk = gkbar
    }
    ik = gk * (v - ek)
}

PROCEDURE states() {
    rates(v)
    n = n + dt * (ninf - n) / ntau
    : channel noise - start
    z1_kpkjslow = z1_kpkjslow*mu1_kpkjslow + noise1_kpkjslow
    z2_kpkjslow = z2_kpkjslow*mu2_kpkjslow + noise2_kpkjslow
    z3_kpkjslow = z3_kpkjslow*mu3_kpkjslow + noise3_kpkjslow
    z4_kpkjslow = z4_kpkjslow*mu4_kpkjslow + noise4_kpkjslow
    : channel noise - end
}

PROCEDURE rates(Vm (mV)) {
    LOCAL v,n4,one_minus_n
    v = Vm + 11	: Account for Junction Potential
    ninf = 1/(1+exp(-(v-nivh)/nik))
    ntau = 1000 * ntau_func(v)
    
    : channel noise - start
    tau1_kpkjslow = ntau
    tau2_kpkjslow = 0.5 * ntau
    tau3_kpkjslow = 0.3333333 * ntau
    tau4_kpkjslow = 0.25 * ntau
    
    mu1_kpkjslow = exp(-dt/tau1_kpkjslow)
    mu2_kpkjslow = exp(-dt/tau2_kpkjslow)
    mu3_kpkjslow = exp(-dt/tau3_kpkjslow)
    mu4_kpkjslow = exp(-dt/tau4_kpkjslow)
    
    n4 = ninf*ninf*ninf*ninf
    one_minus_n = 1. - ninf
    sigma1_kpkjslow = one_over_Nk * 4*n4*ninf*ninf*ninf * one_minus_n
    sigma2_kpkjslow = one_over_Nk * 6*n4*ninf*ninf * one_minus_n*one_minus_n
    sigma3_kpkjslow = one_over_Nk * 4*n4*n * one_minus_n*one_minus_n*one_minus_n
    sigma4_kpkjslow = one_over_Nk * n4 * one_minus_n*one_minus_n*one_minus_n*one_minus_n
    
    noise1_kpkjslow = sqrt(sigma1_kpkjslow * (1-mu1_kpkjslow*mu1_kpkjslow)) * normrand(0,1)
    noise2_kpkjslow = sqrt(sigma2_kpkjslow * (1-mu2_kpkjslow*mu2_kpkjslow)) * normrand(0,1)
    noise3_kpkjslow = sqrt(sigma3_kpkjslow * (1-mu3_kpkjslow*mu3_kpkjslow)) * normrand(0,1)
    noise4_kpkjslow = sqrt(sigma4_kpkjslow * (1-mu4_kpkjslow*mu4_kpkjslow)) * normrand(0,1)

   : channel noise - end    
}

FUNCTION ntau_func(v (mV)) {
    ntau_func = .000796 + 1/(exp((v+73.2)/11.7)+exp((v-306.7)/-74.2))
}
