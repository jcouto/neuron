TITLE Calcium ion accumulation and diffusion with pump and buffer

: Ca2+ lateral diffussion, pump and buffering
: this version 02-2005
: by Klaus M. Stiefel, CNL, The Salk Institute
:
: Thanks to Michael Hines
: Written to model Cerebellar Purkinje cell Ca2+-dynamics
: as described by Schmidt et al., J.Physiology, 2004


NEURON {
	SUFFIX ca_dynamics
	
	USEION ca READ  ica, cai WRITE cai
	USEION mg READ  mgi WRITE mgi VALENCE 2
	RANGE pmp_0, Dca, cainull
	RANGE Dog, OGnull
	RANGE Dcbd1, CBnull
	RANGE Dpar, PVnull	
}

UNITS {
	(mV)	= 	(millivolt)
	(molar) = 	(1/liter)
	(mM)	= 	(millimolar)
	(um)	= 	(micron)
	(mA)	= 	(milliamp)
	(mol)	= 	(1)
	FARADAY = 	(faraday)	 (coulomb)
	PI	= 	(pi)		(1)
	R 	= 	(k-mole)	(joule/degC)
		
}

ASSIGNED {
	ica		(mA/cm2)
	eca		(mV)

	diam		(um)
	area		(um2)

	correct		(/um)

	steadystate	(mM)

	o1		(um2/ms mM)
	o2		(um2/ms)
	
	nf1		(um2/ms mM)
	nf2		(um2/ms)
	ns1		(um2/ms mM)
	ns2		(um2/ms)


	m1		(um2/ms mM)
	m2		(um2/ms)
	p1		(um2/ms mM)
	p2		(um2/ms)
	
	kpmp1		(um2/ms mM)
	kpmp2		(um2/ms)
	kpmp3		(um2/ms)
}

PARAMETER {

	celsius	=37	(degC)
	cainull	=4.5e-5	(mM)
	mginull	=.59	(mM)

	Dca	= .233	(um2/ms)
	Dog	= .014	(um2/ms)
	Dcbd1	= .028	(um2/ms)	
	Dcbd2	= 0	(um2/ms)	
	Dpar	= .043	(um2/ms)

:	values for Oregon Green BAPTA-1
	
	OGnull 	=.16		(mM)
	o1_0	=4.3e2		(/ms mM)
	o2_0	=.14		(/ms)

:	values for Calbindin
:	2 high and 2 low affinity binding sites

	CBnull=.16		(mM)
	nf1_0	=43.5		(/ms mM)
	nf2_0	=3.58e-2	(/ms)
	ns1_0	=5.5		(/ms mM)
	ns2_0	=0.26e-2	(/ms)


:	values for Parvalbumin
	
	PVnull 	= .08		(mM)
	m1_0	= 1.07e2	(/ms mM)
	m2_0	= 9.5e-4		(/ms)
	p1_0	= 0.8		(/ms mM)
	p2_0	= 2.5e-2		(/ms)

:	values for the pump
:	pump is surface based
:	parameters from Schmidt et al., 2003

	pmp_0	=	1e-6	(mM)
	kpmp1_0 =	3e6	(/ms mM)
	kpmp2_0 =	1.75e3	(/ms)
	kpmp3_0 = 	7.255e3	(/ms)

}

CONSTANT {
	cao	= 2.4	(mM) : modified from 2.0
	cor	= 1	(um)

}

STATE {
	cai		(mM)
	mgi		(mM)

:	Oregon green BAPTA-1	
	OG		(mM)
	OG_ca		(mM)	

:	Calbindin. There is a fast and a slow binding site
	CB		(mM)
	CB_f_ca		(mM)
	CB_ca_s		(mM)
	CB_ca_ca	(mM)
:	immobile fraction		
	iCB		(mM)
	iCB_f_ca	(mM)
	iCB_ca_s	(mM)
	iCB_ca_ca	(mM)	

:	Parvalbumin. Binds ca and mg
	PV		(mM)
	PV_ca		(mM)
	PV_mg		(mM)
	
:	Pump
	pmp		(mM)
	capmp		(mM)		
}

INITIAL {
	:	params()
	:	call it from hoc
	
	: since cai is a state it is set to 0 be default
	: thus make sure it is set properly to the external ion value
	VERBATIM
	cai = _ion_cai;
	mgi = _ion_mgi;
	ENDVERBATIM	
		
}

BREAKPOINT {
	
	SOLVE ca_changes METHOD sparse
}
	
KINETIC ca_changes {
	
	
	COMPARTMENT (PI*(diam/2)^2) {cai cao mgi OG OG_ca CB CB_f_ca CB_ca_s CB_ca_ca iCB iCB_f_ca iCB_ca_s iCB_ca_ca PV PV_ca PV_mg pmp capmp} 
:	LONGITUDINAL_DIFFUSION (PI*(diam/2)^2)*Dca {cai mgi}
:	LONGITUDINAL_DIFFUSION (PI*(diam/2)^2)*Dog {OG OG_ca}
:	LONGITUDINAL_DIFFUSION (PI*(diam/2)^2)*Dcbd1 {CB CB_f_ca CB_ca_s CB_ca_ca}
:	LONGITUDINAL_DIFFUSION (PI*(diam/2)^2)*Dcbd2 {iCB iCB_f_ca iCB_ca_s iCB_ca_ca}
:	LONGITUDINAL_DIFFUSION (PI*(diam/2)^2)*Dpar {PV PV_ca PV_mg}

	~cai << (-ica*PI*diam*(1e4)/(2*FARADAY))

:	Pump
	~ cai + pmp <-> capmp (kpmp1, kpmp2)
	~ capmp <-> pmp (kpmp3, 0)
	
:	OG BAPTA-1
	~ cai + OG <-> OG_ca (o1, o2)

:	CB
: 	a mobile and an immobile proportion
:	a fast and a slow binding site
:	lots of equations
	~ cai + CB <-> CB_ca_s (nf1, nf2)
	~ cai + CB <-> CB_f_ca (ns1, ns2)
	~ cai + CB_f_ca <-> CB_ca_ca (nf1, nf2)
	~ cai + CB_ca_s <-> CB_ca_ca (ns1, ns2)

	~ cai + iCB <-> iCB_ca_s (nf1, nf2)
	~ cai + iCB <-> iCB_f_ca (ns1, ns2)
	~ cai + iCB_f_ca <-> iCB_ca_ca (nf1, nf2)
	~ cai + iCB_ca_s <-> iCB_ca_ca (ns1, ns2)


:	PV
	~ cai + PV <-> PV_ca (m1, m2)
	~ mgi + PV <-> PV_mg (p1, p2)		
}

PROCEDURE params() {	

:	if the pumps are per surface	
:	*(diam/2)^2*PI to correct for COMPARTMENT, /(diam/2) because cai is per vol, pump per surface
:	net correction *(diam/2)^2*PI/(diam/2)=*diam*PI/2
: 	
:	if they are /surface and /volume, choose exponent 1 < exp < 2.

	kpmp1 = kpmp1_0*PI*diam/2
	kpmp2 = kpmp2_0*PI*diam/2
	kpmp3 = kpmp3_0*PI*diam/2
	
:	buffers per volume
:	*(diam/2)^2*PI to correct for COMPARTMENT

	o1 = o1_0*(PI*(diam/2)^2)
	o2 = o2_0*(PI*(diam/2)^2)

	nf1 = nf1_0*(PI*(diam/2)^2)
	nf2 = nf2_0*(PI*(diam/2)^2)
	ns1 = ns1_0*(PI*(diam/2)^2)
	ns2 = ns2_0*(PI*(diam/2)^2)


	m1 = m1_0*(PI*(diam/2)^2)
	m2 = m2_0*(PI*(diam/2)^2)
	p1 = p1_0*(PI*(diam/2)^2)
	p2 = p2_0*(PI*(diam/2)^2)

}

FUNCTION ssOG() {
	LOCAL alpha
	alpha=o2/(o1*cainull)
	 steadystate=(OGnull*alpha)/(1+alpha)

	VERBATIM
	return steadystate;
	ENDVERBATIM
}	


FUNCTION ssCBfast() {
	LOCAL alpha
	alpha=nf2/(nf1*cainull)
	 steadystate=(CBnull*alpha)/(1+alpha)

	VERBATIM
	return steadystate;
	ENDVERBATIM
}	


FUNCTION ssCBslow (){
	LOCAL alpha
	alpha=ns2/(ns1*cainull)
	 steadystate=(CBnull*alpha)/(1+alpha)

	VERBATIM
	return steadystate;
	ENDVERBATIM
}	

FUNCTION ssPV() {
	LOCAL alpha, beta
	alpha=(cainull*m1)/m2
	beta=(mginull*p1)/p2
	steadystate=PVnull/(1+alpha+beta)

	VERBATIM
	return steadystate;
	ENDVERBATIM
}

FUNCTION ssPVca() {
	LOCAL alpha, beta
	alpha=(cainull*m1)/m2
	steadystate=alpha

	VERBATIM
	return steadystate;
	ENDVERBATIM
}
