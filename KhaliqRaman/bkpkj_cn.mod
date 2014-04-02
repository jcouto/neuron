: BK-type Purkinje calcium-activated potassium current
: Created 8/19/02 - nwg
: Modified to include channel noise 12/12/11 - Daniele Linaro

NEURON {
    SUFFIX bkpkj_cn
    USEION k READ ek WRITE ik
    USEION ca READ cai
    RANGE gkbar, ik
    GLOBAL minf, mtau, hinf, htau, zinf, ztau
    GLOBAL m_vh, m_k, mtau_y0, mtau_vh1, mtau_vh2, mtau_k1, mtau_k2
    GLOBAL z_coef, ztau
    GLOBAL h_y0, h_vh, h_k, htau_y0, htau_vh1, htau_vh2, htau_k1, htau_k2
    : channel noise - start
    RANGE gk, gamma_k
    RANGE Nk, one_over_Nk
    RANGE seed
    : channel noise - end
}

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
    : channel noise - start
    (S) = (siemens)
    (pS) = (picosiemens)
    : channel noise - end
}

PARAMETER {
    v            (mV)
    gkbar = .007 (mho/cm2)
    
    m_vh = -28.9           (mV)
    m_k = 6.2            (mV)
    mtau_y0 = .000505     (s)
    mtau_vh1 = -33.3     (mV)
    mtau_k1 = -10         (mV)
    mtau_vh2 = 86.4       (mV)
    mtau_k2 = 10.1        (mV)
    
    z_coef = .001        (mM)
    ztau = 1              (ms)
    
    h_y0 = .085
    h_vh = -32          (mV)
    h_k = 5.8             (mV)
    htau_y0 = .0019      (s)
    htau_vh1 = -54.2       (mV)
    htau_k1 = -12.9       (mV)
    htau_vh2 = 48.5      (mV)
    htau_k2 = 5.2        (mV)
    
    ek           (mV)
    cai          (mM)
    
    : channel noise - start
    seed = 5061983 (1)
    gamma_k = 10 (pS)
    : channel noise - end
}

ASSIGNED {
    minf
    mtau          (ms)
    hinf
    htau          (ms)
    zinf
    
    ik            (mA/cm2)
    
    : channel noise - start

    gk (S/cm2)
    Nk (1)
    one_over_Nk (1)
    
    dt (ms)
    area (um2)
    
    tau1_bkpkj (ms) tau2_bkpkj (ms) tau3_bkpkj (ms) tau4_bkpkj (ms) tau5_bkpkj (ms) tau6_bkpkj (ms)
    tau7_bkpkj (ms) tau8_bkpkj (ms) tau9_bkpkj (ms) tau10_bkpkj (ms) tau11_bkpkj (ms) tau12_bkpkj (ms)
    tau13_bkpkj (ms) tau14_bkpkj (ms) tau15_bkpkj (ms) tau16_bkpkj (ms) tau17_bkpkj (ms) tau18_bkpkj (ms)
    tau19_bkpkj (ms) tau20_bkpkj (ms) tau21_bkpkj (ms) tau22_bkpkj (ms) tau23_bkpkj (ms)

    sigma1_bkpkj (ms2) sigma2_bkpkj (ms2) sigma3_bkpkj (ms2) sigma4_bkpkj (ms2) sigma5_bkpkj (ms2) sigma6_bkpkj (ms2)
    sigma7_bkpkj (ms2) sigma8_bkpkj (ms2) sigma9_bkpkj (ms2) sigma10_bkpkj (ms2) sigma11_bkpkj (ms2) sigma12_bkpkj (ms2)
    sigma13_bkpkj (ms2) sigma14_bkpkj (ms2) sigma15_bkpkj (ms2) sigma16_bkpkj (ms2) sigma17_bkpkj (ms2) sigma18_bkpkj (ms2)
    sigma19_bkpkj (ms2) sigma20_bkpkj (ms2) sigma21_bkpkj (ms2) sigma22_bkpkj (ms2) sigma23_bkpkj (ms2)

    noise1_bkpkj noise2_bkpkj noise3_bkpkj noise4_bkpkj noise5_bkpkj noise6_bkpkj
    noise7_bkpkj noise8_bkpkj noise9_bkpkj noise10_bkpkj noise11_bkpkj noise12_bkpkj
    noise13_bkpkj noise14_bkpkj noise15_bkpkj noise16_bkpkj noise17_bkpkj noise18_bkpkj
    noise19_bkpkj noise20_bkpkj noise21_bkpkj noise22_bkpkj noise23_bkpkj

    mu1_bkpkj mu2_bkpkj mu3_bkpkj mu4_bkpkj mu5_bkpkj mu6_bkpkj
    mu7_bkpkj mu8_bkpkj mu9_bkpkj mu10_bkpkj mu11_bkpkj mu12_bkpkj
    mu13_bkpkj mu14_bkpkj mu15_bkpkj mu16_bkpkj mu17_bkpkj mu18_bkpkj
    mu19_bkpkj mu20_bkpkj mu21_bkpkj mu22_bkpkj mu23_bkpkj
    
    : channel noise - end
}

STATE {
    m   FROM 0 TO 1
    z   FROM 0 TO 1
    h   FROM 0 TO 1
    : channel noise - start
    z1_bkpkj z2_bkpkj z3_bkpkj z4_bkpkj z5_bkpkj z6_bkpkj 
    z7_bkpkj z8_bkpkj z9_bkpkj z10_bkpkj z11_bkpkj z12_bkpkj 
    z13_bkpkj z14_bkpkj z15_bkpkj z16_bkpkj z17_bkpkj z18_bkpkj 
    z19_bkpkj z20_bkpkj z21_bkpkj z22_bkpkj z23_bkpkj
    : channel noise - end
}

BREAKPOINT {
    LOCAL Z
    SOLVE states
    Z =     z1_bkpkj+z2_bkpkj+z3_bkpkj+z4_bkpkj+z5_bkpkj+z6_bkpkj
    Z = Z + z7_bkpkj+z8_bkpkj+z9_bkpkj+z10_bkpkj+z11_bkpkj+z12_bkpkj
    Z = Z + z13_bkpkj+z14_bkpkj+z15_bkpkj+z16_bkpkj+z17_bkpkj+z18_bkpkj
    Z = Z + z19_bkpkj+z20_bkpkj+z21_bkpkj+z21_bkpkj+z23_bkpkj
    gk = gkbar * (m*m*m*z*z*h + Z)
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
    m = m + dt * (minf-m) / mtau
    h = h + dt * (hinf-h) / htau
    z = z + dt * (zinf-z) / ztau

    : channel noise - start
    z1_bkpkj = z1_bkpkj*mu1_bkpkj + noise1_bkpkj
    z2_bkpkj = z2_bkpkj*mu2_bkpkj + noise2_bkpkj
    z3_bkpkj = z3_bkpkj*mu3_bkpkj + noise3_bkpkj
    z4_bkpkj = z4_bkpkj*mu4_bkpkj + noise4_bkpkj
    z5_bkpkj = z5_bkpkj*mu5_bkpkj + noise5_bkpkj
    z6_bkpkj = z6_bkpkj*mu6_bkpkj + noise6_bkpkj
    z7_bkpkj = z7_bkpkj*mu7_bkpkj + noise7_bkpkj
    z8_bkpkj = z8_bkpkj*mu8_bkpkj + noise8_bkpkj
    z9_bkpkj = z9_bkpkj*mu9_bkpkj + noise9_bkpkj
    z10_bkpkj = z10_bkpkj*mu10_bkpkj + noise10_bkpkj
    z11_bkpkj = z11_bkpkj*mu11_bkpkj + noise11_bkpkj
    z12_bkpkj = z12_bkpkj*mu12_bkpkj + noise12_bkpkj
    z13_bkpkj = z13_bkpkj*mu13_bkpkj + noise13_bkpkj
    z14_bkpkj = z14_bkpkj*mu14_bkpkj + noise14_bkpkj
    z15_bkpkj = z15_bkpkj*mu15_bkpkj + noise15_bkpkj
    z16_bkpkj = z16_bkpkj*mu16_bkpkj + noise16_bkpkj
    z17_bkpkj = z17_bkpkj*mu17_bkpkj + noise17_bkpkj
    z18_bkpkj = z18_bkpkj*mu18_bkpkj + noise18_bkpkj
    z19_bkpkj = z19_bkpkj*mu19_bkpkj + noise19_bkpkj
    z20_bkpkj = z20_bkpkj*mu20_bkpkj + noise20_bkpkj
    z21_bkpkj = z21_bkpkj*mu21_bkpkj + noise21_bkpkj
    z22_bkpkj = z22_bkpkj*mu22_bkpkj + noise22_bkpkj
    z23_bkpkj = z23_bkpkj*mu23_bkpkj + noise23_bkpkj
    : channel noise - end
}

PROCEDURE rates(Vm (mV)) {
    LOCAL v,m2,m3,m6,z2,h2,one_minus_m,one_minus_z,one_minus_h

    v = Vm + 5
    minf = 1 / (1 + exp(-(v - (m_vh)) / m_k))
    mtau = (1e3) * (mtau_y0 + 1/(exp((v+ mtau_vh1)/mtau_k1) + exp((v+mtau_vh2)/mtau_k2)))
    zinf = 1/(1 + z_coef / cai)
    hinf = h_y0 + (1-h_y0) / (1+exp((v - h_vh)/h_k))
    htau = (1e3) * (htau_y0 + 1/(exp((v + htau_vh1)/htau_k1)+exp((v+htau_vh2)/htau_k2)))
    
    : channel noise - start
    tau1_bkpkj = htau
    tau2_bkpkj = 0.5*ztau
    tau3_bkpkj = (htau*ztau) / (2*htau+ztau)
    tau4_bkpkj = ztau
    tau5_bkpkj = (htau*ztau) / (htau+ztau)
    tau6_bkpkj = 0.3333333 * mtau
    tau7_bkpkj = (htau*mtau) / (3*htau+mtau)
    tau8_bkpkj = (mtau*ztau) / (mtau+3*ztau)
    tau9_bkpkj = (htau*mtau*ztau) / (htau*mtau + 3*htau*ztau + mtau*ztau)
    tau10_bkpkj = (mtau*ztau) / (2*mtau+3*ztau)
    tau11_bkpkj = (htau*mtau*ztau) / (2*htau*mtau + 3*htau*ztau + mtau*ztau)
    tau12_bkpkj = mtau
    tau13_bkpkj = (htau*mtau) / (htau+mtau)
    tau14_bkpkj = (mtau*ztau) / (mtau+ztau)
    tau15_bkpkj = (htau*mtau*ztau) / (htau*mtau + htau*ztau + mtau*ztau)
    tau16_bkpkj = (mtau*ztau) / (2*mtau+ztau)
    tau17_bkpkj = (htau*mtau*ztau) / (2*htau*mtau + htau*ztau + mtau*ztau)
    tau18_bkpkj = 0.5*mtau
    tau19_bkpkj = (htau*mtau) / (2*htau+mtau)
    tau20_bkpkj = 0.5 * (mtau*ztau) / (mtau+ztau)
    tau21_bkpkj = (htau*mtau*ztau) / (2*htau*mtau + 2*htau*ztau + mtau*ztau)
    tau22_bkpkj = (mtau*ztau) / (mtau+2*ztau)
    tau23_bkpkj = (htau*mtau*ztau) / (htau*mtau + 2*htau*ztau + mtau*ztau)
    
    m2 = minf*minf
    m3 = m2*minf
    m6 = m3*m3
    z2 = zinf*zinf
    h2 = hinf*hinf
    one_minus_m = 1. - minf
    one_minus_z = 1. - zinf
    one_minus_h = 1. - hinf
    sigma1_bkpkj = one_over_Nk * m6*z2*z2*h*one_minus_h
    sigma2_bkpkj = one_over_Nk * m6*z2*h2*one_minus_z*one_minus_z
    sigma3_bkpkj = one_over_Nk * m6*z2*h*one_minus_z*one_minus_z*one_minus_h
    sigma4_bkpkj = one_over_Nk * 2*m6*z2*z*h2*one_minus_z
    sigma5_bkpkj = one_over_Nk * 2*m6*z2*z*h*one_minus_z*one_minus_h
    sigma6_bkpkj = one_over_Nk * m3*z2*z2*h2*one_minus_m*one_minus_m*one_minus_m
    sigma7_bkpkj = one_over_Nk * m3*z2*z2*h*one_minus_m*one_minus_m*one_minus_m*one_minus_h
    sigma8_bkpkj = one_over_Nk * 2*m3*z2*z*h2*one_minus_m*one_minus_m*one_minus_m*one_minus_z
    sigma9_bkpkj = one_over_Nk * 2*m3*z2*z*h*one_minus_m*one_minus_m*one_minus_m*one_minus_z*one_minus_h
    sigma10_bkpkj = one_over_Nk * m3*z2*h2*one_minus_m*one_minus_m*one_minus_m*one_minus_z*one_minus_z
    sigma11_bkpkj = one_over_Nk * m3*z2*h*one_minus_m*one_minus_m*one_minus_m*one_minus_z*one_minus_z*one_minus_h
    sigma12_bkpkj = one_over_Nk * 3*m3*m2*z2*z2*h2*one_minus_m
    sigma13_bkpkj = one_over_Nk * 3*m3*m2*z2*z2*h*one_minus_m*one_minus_h
    sigma14_bkpkj = one_over_Nk * 6*m3*m2*z2*z*h2*one_minus_m*one_minus_z
    sigma15_bkpkj = one_over_Nk * 6*m3*m2*z2*z*h*one_minus_m*one_minus_z*one_minus_h
    sigma16_bkpkj = one_over_Nk * 3*m3*m2*z2*h2*one_minus_m*one_minus_z*one_minus_z
    sigma17_bkpkj = one_over_Nk * 3*m3*m2*z2*h*one_minus_m*one_minus_z*one_minus_z*one_minus_h
    sigma18_bkpkj = one_over_Nk * 3*m2*m2*z2*z2*h2*one_minus_m*one_minus_m
    sigma19_bkpkj = one_over_Nk * 3*m2*m2*z2*z2*h*one_minus_m*one_minus_m*one_minus_h
    sigma20_bkpkj = one_over_Nk * 3*m2*m2*z2*h2*one_minus_m*one_minus_m*one_minus_z*one_minus_z
    sigma21_bkpkj = one_over_Nk * 3*m2*m2*z2*h*one_minus_m*one_minus_m*one_minus_z*one_minus_z*one_minus_h
    sigma22_bkpkj = one_over_Nk * 6*m2*m2*z2*z*h2*one_minus_m*one_minus_m*one_minus_z
    sigma23_bkpkj = one_over_Nk * 6*m2*m2*z2*z*h*one_minus_m*one_minus_m*one_minus_z*one_minus_h
    
    mu1_bkpkj = exp(-dt/tau1_bkpkj)
    mu2_bkpkj = exp(-dt/tau2_bkpkj)
    mu3_bkpkj = exp(-dt/tau3_bkpkj)
    mu4_bkpkj = exp(-dt/tau4_bkpkj)
    mu5_bkpkj = exp(-dt/tau5_bkpkj)
    mu6_bkpkj = exp(-dt/tau6_bkpkj)
    mu7_bkpkj = exp(-dt/tau7_bkpkj)
    mu8_bkpkj = exp(-dt/tau8_bkpkj)
    mu9_bkpkj = exp(-dt/tau9_bkpkj)
    mu10_bkpkj = exp(-dt/tau10_bkpkj)
    mu11_bkpkj = exp(-dt/tau11_bkpkj)
    mu12_bkpkj = exp(-dt/tau12_bkpkj)
    mu13_bkpkj = exp(-dt/tau13_bkpkj)
    mu14_bkpkj = exp(-dt/tau14_bkpkj)
    mu15_bkpkj = exp(-dt/tau15_bkpkj)
    mu16_bkpkj = exp(-dt/tau16_bkpkj)
    mu17_bkpkj = exp(-dt/tau17_bkpkj)
    mu18_bkpkj = exp(-dt/tau18_bkpkj)
    mu19_bkpkj = exp(-dt/tau19_bkpkj)
    mu20_bkpkj = exp(-dt/tau20_bkpkj)
    mu21_bkpkj = exp(-dt/tau21_bkpkj)
    mu22_bkpkj = exp(-dt/tau22_bkpkj)
    mu23_bkpkj = exp(-dt/tau23_bkpkj)

    noise1_bkpkj = sqrt(sigma1_bkpkj * (1-mu1_bkpkj*mu1_bkpkj)) * normrand(0,1)
    noise2_bkpkj = sqrt(sigma2_bkpkj * (1-mu2_bkpkj*mu2_bkpkj)) * normrand(0,1)
    noise3_bkpkj = sqrt(sigma3_bkpkj * (1-mu3_bkpkj*mu3_bkpkj)) * normrand(0,1)
    noise4_bkpkj = sqrt(sigma4_bkpkj * (1-mu4_bkpkj*mu4_bkpkj)) * normrand(0,1)
    noise5_bkpkj = sqrt(sigma5_bkpkj * (1-mu5_bkpkj*mu5_bkpkj)) * normrand(0,1)
    noise6_bkpkj = sqrt(sigma6_bkpkj * (1-mu6_bkpkj*mu6_bkpkj)) * normrand(0,1)
    noise7_bkpkj = sqrt(sigma7_bkpkj * (1-mu7_bkpkj*mu7_bkpkj)) * normrand(0,1)
    noise8_bkpkj = sqrt(sigma8_bkpkj * (1-mu8_bkpkj*mu8_bkpkj)) * normrand(0,1)
    noise9_bkpkj = sqrt(sigma9_bkpkj * (1-mu9_bkpkj*mu9_bkpkj)) * normrand(0,1)
    noise10_bkpkj = sqrt(sigma10_bkpkj * (1-mu10_bkpkj*mu10_bkpkj)) * normrand(0,1)
    noise11_bkpkj = sqrt(sigma11_bkpkj * (1-mu11_bkpkj*mu11_bkpkj)) * normrand(0,1)
    noise12_bkpkj = sqrt(sigma12_bkpkj * (1-mu12_bkpkj*mu12_bkpkj)) * normrand(0,1)
    noise13_bkpkj = sqrt(sigma13_bkpkj * (1-mu13_bkpkj*mu13_bkpkj)) * normrand(0,1)
    noise14_bkpkj = sqrt(sigma14_bkpkj * (1-mu14_bkpkj*mu14_bkpkj)) * normrand(0,1)
    noise15_bkpkj = sqrt(sigma15_bkpkj * (1-mu15_bkpkj*mu15_bkpkj)) * normrand(0,1)
    noise16_bkpkj = sqrt(sigma16_bkpkj * (1-mu16_bkpkj*mu16_bkpkj)) * normrand(0,1)
    noise17_bkpkj = sqrt(sigma17_bkpkj * (1-mu17_bkpkj*mu17_bkpkj)) * normrand(0,1)
    noise18_bkpkj = sqrt(sigma18_bkpkj * (1-mu18_bkpkj*mu18_bkpkj)) * normrand(0,1)
    noise19_bkpkj = sqrt(sigma19_bkpkj * (1-mu19_bkpkj*mu19_bkpkj)) * normrand(0,1)
    noise20_bkpkj = sqrt(sigma20_bkpkj * (1-mu20_bkpkj*mu20_bkpkj)) * normrand(0,1)
    noise21_bkpkj = sqrt(sigma21_bkpkj * (1-mu21_bkpkj*mu21_bkpkj)) * normrand(0,1)
    noise22_bkpkj = sqrt(sigma22_bkpkj * (1-mu22_bkpkj*mu22_bkpkj)) * normrand(0,1)
    noise23_bkpkj = sqrt(sigma23_bkpkj * (1-mu23_bkpkj*mu23_bkpkj)) * normrand(0,1)
    : channel noise - end
}

INITIAL {
    rates(v)
    m = minf
    z = zinf
    h = hinf

    : channel noise - start
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))
    one_over_Nk = 1.0 / Nk
    printf("bkpkj>> the number of channels is %.0f.\n", Nk)
    z1_bkpkj = 0.
    z2_bkpkj = 0.
    z3_bkpkj = 0.
    z4_bkpkj = 0.
    z5_bkpkj = 0.
    z6_bkpkj = 0.
    z7_bkpkj = 0.
    z8_bkpkj = 0.
    z9_bkpkj = 0.
    z10_bkpkj = 0.
    z11_bkpkj = 0.
    z12_bkpkj = 0.
    z13_bkpkj = 0.
    z14_bkpkj = 0.
    z15_bkpkj = 0.
    z16_bkpkj = 0.
    z17_bkpkj = 0.
    z18_bkpkj = 0.
    z19_bkpkj = 0.
    z20_bkpkj = 0.
    z21_bkpkj = 0.
    z22_bkpkj = 0.
    z23_bkpkj = 0.
    set_seed(seed)
    : channel noise - end
}
