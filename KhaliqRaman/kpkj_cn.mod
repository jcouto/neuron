: HH TEA-sensitive Purkinje potassium current (including channel noise)
: Created: 7th October 2011 - Daniele Linaro

NEURON {
    SUFFIX kpkj_cn
    USEION k READ ek WRITE ik
    RANGE gkbar, ik
    GLOBAL minf, hinf, mtau, htau
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
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micrometer)
    : channel noise - end
}

PARAMETER {
    gkbar = .004 (mho/cm2)
    
    mivh = -24 (mV)
    mik = 15.4 (1)
    mty0 = .00012851 	
    mtvh1 = 100.7 (mV)
    mtk1 = 12.9 (1)
    mtvh2 = -56.0 (mV)
    mtk2 = -23.1 (1)
    
    hiy0 = .31	
    hiA = .78
    hivh = -5.802 (mV)
    hik = 11.2 (1)
    
    ek (mV)
    
    : channel noise - start
    gamma_k = 10 (pS)
    seed = 5061983 (1)
    : channel noise - end
}

ASSIGNED {
    v (mV)
    
    ik (mA/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
    
    : channel noise - start
    gk (S/cm2)
    Nk (1)
    one_over_Nk (1)
   
    dt (ms)
    area (um2)
    
    tau1_kpkj (ms) tau2_kpkj (ms) tau3_kpkj (ms) tau4_kpkj (ms) tau5_kpkj (ms) tau6_kpkj (ms) tau7_kpkj (ms) 
    sigma1_kpkj (ms2) sigma2_kpkj (ms2) sigma3_kpkj (ms2) sigma4_kpkj (ms2) sigma5_kpkj (ms2) sigma6_kpkj (ms2) sigma7_kpkj (ms2)
    mu1_kpkj mu2_kpkj mu3_kpkj mu4_kpkj mu5_kpkj mu6_kpkj mu7_kpkj
    noise1_kpkj noise2_kpkj noise3_kpkj noise4_kpkj noise5_kpkj noise6_kpkj noise7_kpkj
    : channel noise - end
}

STATE {
    m h
    : channel noise  - start
    z1_kpkj z2_kpkj z3_kpkj z4_kpkj z5_kpkj z6_kpkj z7_kpkj
    : channel noise - end
}

INITIAL {
    rates(v)
    m = minf
    h = hinf
    
    : channel noise - start
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))
    one_over_Nk = 1.0 / Nk
    printf("kpkj>> the number of channels is %.0f.\n", Nk)
    z1_kpkj = 0.0
    z2_kpkj = 0.0
    z3_kpkj = 0.0
    z4_kpkj = 0.0
    z5_kpkj = 0.0
    z6_kpkj = 0.0
    z7_kpkj = 0.0
    set_seed(seed)
    : channel noise - end
}

BREAKPOINT {
    SOLVE states
    gk = gkbar * (m*m*m*h + z1_kpkj+z2_kpkj+z3_kpkj+z4_kpkj+z5_kpkj+z6_kpkj+z7_kpkj)
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
    m = m + dt * (minf - m) / mtau
    h = h + dt * (hinf - h) / htau

    : channel noise - start
    z1_kpkj = z1_kpkj*mu1_kpkj + noise1_kpkj
    z2_kpkj = z2_kpkj*mu2_kpkj + noise2_kpkj
    z3_kpkj = z3_kpkj*mu3_kpkj + noise3_kpkj
    z4_kpkj = z4_kpkj*mu4_kpkj + noise4_kpkj
    z5_kpkj = z5_kpkj*mu5_kpkj + noise5_kpkj
    z6_kpkj = z6_kpkj*mu6_kpkj + noise6_kpkj
    z7_kpkj = z7_kpkj*mu7_kpkj + noise7_kpkj
    : channel noise - end
}

PROCEDURE rates( Vm (mV)) {
    LOCAL lv,m3inf,one_minus_m,one_minus_h
    lv = Vm + 11 : Account for Junction Potential
    
    minf = 1.0 / (1+exp(-(lv-mivh)/mik))
    mtau = 1000 * mtau_func(lv)
    hinf = hiy0 + hiA/(1+exp((lv-hivh)/hik))
    htau = 1000 * htau_func(lv)
    
    : weird check that shouldn't be necessary, but at times it is...
    if (minf < 0) {
	minf = 0.
    } else if (minf > 1) {
	minf = 1.
    }
    if (hinf < 0) {
	hinf = 0.
    } else if (hinf > 1) {
	hinf = 1.
    }

    : channel noise - start
    
    m3inf = minf*minf*minf
    one_minus_m = 1. - minf
    one_minus_h = 1. - hinf
    tau1_kpkj = htau
    tau2_kpkj = mtau
    tau3_kpkj = 0.5*mtau
    tau4_kpkj = 0.3333333*mtau
    tau5_kpkj = mtau*htau/(mtau+htau)
    tau6_kpkj = mtau*htau/(mtau+2*htau)
    tau7_kpkj = mtau*htau/(mtau+3*htau)
    sigma1_kpkj = one_over_Nk * m3inf*m3inf*hinf * one_minus_h
    sigma2_kpkj = one_over_Nk * 3*m3inf*minf*minf*hinf*hinf * one_minus_m
    sigma3_kpkj = one_over_Nk * 3*m3inf*minf*hinf*hinf * one_minus_m*one_minus_m
    sigma4_kpkj = one_over_Nk * m3inf*hinf*hinf * one_minus_m*one_minus_m*one_minus_m
    sigma5_kpkj = one_over_Nk * 3*m3inf*minf*minf*hinf * one_minus_m*one_minus_h
    sigma6_kpkj = one_over_Nk * 3*m3inf*minf*hinf * one_minus_m*one_minus_m*one_minus_h
    sigma7_kpkj = one_over_Nk * m3inf*hinf * one_minus_m*one_minus_m*one_minus_m*one_minus_h
    mu1_kpkj = exp(-dt/tau1_kpkj)
    mu2_kpkj = exp(-dt/tau2_kpkj)
    mu3_kpkj = exp(-dt/tau3_kpkj)
    mu4_kpkj = exp(-dt/tau4_kpkj)
    mu5_kpkj = exp(-dt/tau5_kpkj)
    mu6_kpkj = exp(-dt/tau6_kpkj)
    mu7_kpkj = exp(-dt/tau7_kpkj)
    noise1_kpkj = sqrt(sigma1_kpkj * (1-mu1_kpkj*mu1_kpkj)) * normrand(0,1)
    noise2_kpkj = sqrt(sigma2_kpkj * (1-mu2_kpkj*mu2_kpkj)) * normrand(0,1)
    noise3_kpkj = sqrt(sigma3_kpkj * (1-mu3_kpkj*mu3_kpkj)) * normrand(0,1)
    noise4_kpkj = sqrt(sigma4_kpkj * (1-mu4_kpkj*mu4_kpkj)) * normrand(0,1)
    noise5_kpkj = sqrt(sigma5_kpkj * (1-mu5_kpkj*mu5_kpkj)) * normrand(0,1)
    noise6_kpkj = sqrt(sigma6_kpkj * (1-mu6_kpkj*mu6_kpkj)) * normrand(0,1)
    noise7_kpkj = sqrt(sigma7_kpkj * (1-mu7_kpkj*mu7_kpkj)) * normrand(0,1)
    
    : channel noise - end
}

FUNCTION mtau_func (Vm (mV)) (ms) {
    if (Vm < -35) {
	mtau_func = 3*(3.4225e-5+.00498*exp(Vm/28.29))
    } else {
	mtau_func = (mty0 + 1/(exp((Vm+mtvh1)/mtk1) + exp((Vm+mtvh2)/mtk2)))
    }
}

FUNCTION htau_func(Vm (mV)) (ms) {
    if (Vm > 0) {
	htau_func = .0012+.0023*exp(-.141*Vm)
    } else {
	htau_func = 1.2202e-05 + .012 * exp(-((Vm+56.3)/49.6)^2)
    }
}
