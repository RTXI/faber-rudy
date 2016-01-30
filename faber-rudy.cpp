/*
 * Copyright (C) 2011 Weill Medical College of Cornell University
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <iostream>
#include <math.h>
#include "faber-rudy.h"

using namespace std;
#define fastEXP RTmath->fastEXP //define shortcut for fastEXP function

// Required by RTXI
extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new Faber_Rudy_2000();
}

// These variables will automatically populate the Default GUI
// { "GUI LABEL", "GUI TOOLTIP", "VARIABLE TYPE",},
static DefaultGUIModel::variable_t vars[] = {    
	{ "Current Input (A)", "Current Input (A)", DefaultGUIModel::INPUT, }, 
	{ "Voltage Output (V)", "Voltage Output (V)", DefaultGUIModel::OUTPUT, }, 
	{ "Model Rate (Hz)", "Model Rate (Hz)", 
	   DefaultGUIModel::PARAMETER | DefaultGUIModel::INTEGER,},
	{ "gNa", "Sodium Conductance", 
	   DefaultGUIModel::PARAMETER | DefaultGUIModel::DOUBLE,},
	{ "INa", "Sodium Current", DefaultGUIModel::STATE ,} 
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

/*** Model Constructor ***/
Faber_Rudy_2000::Faber_Rudy_2000(void) : DefaultGUIModel("Faber-Rudy 2000", ::vars, ::num_vars) {
	// Function call to create graphical interface *can be overloaded for custom GUI*
	DefaultGUIModel::createGUI(vars, num_vars);
	initialize(); // Custom function for code organization, initializes parameters
	show(); // Show GUI
	resizeMe();
	refresh(); // Refresh GUI
}

/*** Model Destructor Function ***/
Faber_Rudy_2000::~Faber_Rudy_2000(void) {
	// Delete lookup tables
	delete[] lkup;
	delete[] coslkup;
	delete[] acoslkup;
	delete RTmath;
} 

/*** Realtime Execute Function ***/
void Faber_Rudy_2000::execute() {

	// Input in A for more compatibility with traditional current clamp input, 
	// Factor of 1e10 to replicate a whole cell capacitance in pF range
	I_Inject = input(0) * -1 * 1e10;

	// Run solver required number of steps, determined by RTXI thread and model rate
	for(i = 0; i < steps; i++) {
		solve();
	}

	// Output in V for more compatability with traditional amplifier output
	output(0) = V / 1000; 
}

/*** Non-Realtime Update Function ***/
void Faber_Rudy_2000::update(DefaultGUIModel::update_flags_t flag) {

	switch (flag) {
		case UNPAUSE: // Called when pause button is untoggled
			DT = (1.0/modelRate)*1000; // Calculates DT (ms)
			// Grabs RTXI thread period and converts to ms (from ns)
			RTperiod = RT::System::getInstance()->getPeriod()*1e-6; 

			// Determines number of iterations model solver will run based 
			// on RTXI thread rate
			steps = static_cast<int>(ceil(RTperiod/DT)); 
			break;
			
		case PAUSE: // Called when pause button is toggled
			// Model is reset to initial conditions when the modify button is 
			// clicked published initial conditions
			/* 
			   GUINEA_PIG[0]= 0.000003; 
			   GUINEA_PIG[1]= 0.999745;
			   GUINEA_PIG[2] = 0.004503; 
			   GUINEA_PIG[3]= 0.004503;
			   GUINEA_PIG[4]= 0.000129; 
			   GUINEA_PIG[5]= 0.000994; 
			   GUINEA_PIG[6]= 0.994041;
			   GUINEA_PIG[7]= 0.000838; 
			   GUINEA_PIG[8]= 0.993336; 
			   GUINEA_PIG[9]= 0.995484; 
			   GUINEA_PIG[10]= 12.236437; //10;
			   GUINEA_PIG[11]= 136.89149; //145;
			   GUINEA_PIG[12]= 0.000079; //0.12E-3;
			   GUINEA_PIG[13]=1.179991;
			   GUINEA_PIG[14]=1.179991;
			*/

			//Initial conditions after 300 WT HH beats at BCL = 1000ms
			/* 
			   GUINEA_PIG[0]= 0.000004;
			   GUINEA_PIG[1]= 0.999723;
			   GUINEA_PIG[2]= 0.005301;
			   GUINEA_PIG[3]= 0.033794;
			   GUINEA_PIG[4]= 0.000134;
			   GUINEA_PIG[5]= 0.001020;
			   GUINEA_PIG[6]= 0.993579;
			   GUINEA_PIG[7]= 0.000878;
			   GUINEA_PIG[8]= 0.992906;
			   GUINEA_PIG[9]= 0.995226;
			   GUINEA_PIG[10]= 13.454714;
			   GUINEA_PIG[11]= 135.525797;
			   GUINEA_PIG[12]= 0.000122;
			   GUINEA_PIG[13]= 1.531388;
			   GUINEA_PIG[14]=1.942522;
			*/

			// Initial conditions after 1100 beats at 1000ms BCL with new concentrations
			V = -83.4366;
			GUINEA_PIG[0] = 7.71003e-06;
			GUINEA_PIG[1] = 0.999098;
			GUINEA_PIG[2] = 0.0184959;
			GUINEA_PIG[3] = 0.082228;
			GUINEA_PIG[4] = 0.000261177;
			GUINEA_PIG[5] = 0.00160623;
			GUINEA_PIG[6] = 0.980753;
			GUINEA_PIG[7] = 0.00199501;
			GUINEA_PIG[8] = 0.978777;
			GUINEA_PIG[9] = 0.986735;
			GUINEA_PIG[10] = 11.753;
			GUINEA_PIG[11] = 133.987;
			GUINEA_PIG[12] = 0.000147466;
			GUINEA_PIG[13] = 1.46651;
			GUINEA_PIG[14] = 2.31135;
			break;

		case MODIFY: // Called when modify button is hit                
			// Update Paramters
			modelRate = getParameter("Model Rate (Hz)").toInt();
			gna = getParameter("gNa").toDouble();

			DT = (1.0/modelRate)*1000; // Calculates DT (ms)

			// Grabs RTXI thread period and converts to ms (from ns)
			RTperiod = RT::System::getInstance()->getPeriod()*1e-6; 

			// Determines number of iterations model solver will run based on RTXI thread rate
			steps = static_cast<int>(ceil(RTperiod/DT)); 
			break;
	}
}


/*** Parameter Initialization ***/
void Faber_Rudy_2000::initialize() {
	RTmath = new RealTimeMath(); // RealTimeMath Library

	// Initial conditions and Constant initialization
	V = -89.00;
	I_Inject = 0;
	R = 8314.0;
	T = 273 + 37;
	F = 96485.0;
	z_Na = 1;						
	z_Ca = 2;
	z_K = 1;

	//Cell Geometry
	pi = 3.141592;						
	L = 0.01; //Length 100um (cm)
	r = 0.0011; //Radius 11um (cm)
	R_cg = 2; //Ratio of capacitive membrane area to gemometric membrane area
	V_cell = 1000*pi*r*r*L; //Volume of cell 38E-6 uL
	A_geo = 2*pi*r*r + 2*pi*r*L; //Geometric membrane area; 0.767E-4 cm^2
	A_cap = R_cg*A_geo; //Capacitive membrane area
	V_myo = V_cell * 0.68; //Myoplasm volume set to 68% of volume of cell
	V_mito = V_cell * 0.26; //Mitochondria volume set to 26% of volume of cell
	V_sr = V_cell * 0.06; //Sarcoplasmic reticulum volume set to 6%
	V_nsr = V_sr * 0.92; //Network SR is 92% of SR
	V_jsr = V_sr * 0.08; //Junctional SR is 8% of SR
	V_cleft = (V_cell/0.88)*0.12; //Cleft volume
	flag = 0;
	Na_out = 137; // Was 140, changed to match patch solutions
	Ca_out = 2.0; // Was 1.8, changed to match patch solutions
	K_out = 5.4; // Was 4.5, changed to match patch solutions
	dCa_in = 0.0;
	dCa_jsr = 0.0;
	dCa_nsr = 0.0;
	dNa_in = 0.0;
	dK_in = 0.0;
	I_Total = 0.0;
	xs2_ss= 0.0; 
	xs1_ss= 0.0; 
	tau_xs1= 0.0;  
	tau_xs2= 0.0; 
	tau_on=0.5; // constant of activation of Ca release from JSR
	tau_off=0.5; // time constant of deactivation of Ca release from JSR
	tcicr=1000;
	tjsrol=1000;
	TRPN=0.0143923;
	CMDN=0.00257849;
	Ca_CSQN=6.97978;
	GP_err=1E-5;
	gna = 16;
	Km_Ca = 0.6E-3;
	P_Ca = 5.4E-4;
	g_Ca_in = 1;
	g_Ca_out = 0.341;
	P_K = 1.93E-7;
	g_K_in = 0.75;
	g_K_out = 0.75;
	P_Na = 6.75E-7;
	g_Na_in = 0.75;
	g_Na_out = 0.75;
	G_Ca_T = 0.05;
	PR_NaK = 0.01833;
	G_Kp = 0.00552;
	c1 = .00025; // Scaling factor for inaca (uA/uF)
	c2 = 0.0001; // Half-saturation concentration of NaCa exhanger (mM)
	gammas = 0.15; // Position of energy barrier controlling voltage dependance of inaca
	I_Na_K_bar = 2.25;
	Km_Na_i = 10;
	Km_Ko = 1.5;
	I_p_Ca_bar = 1.15;
	Km_p_Ca = 0.0005;
	G_Ca_b_bar = 0.003016;
	G_Na_b_bar = 0.004;//0.00141;
	tau_tr = 180;
	I_up_bar = 0.00875;
	Ca_nsr_bar = 15;
	K_leak = I_up_bar/ Ca_nsr_bar;
	Km_up = 0.00092;
	G_rel_bar = 150; // Max. constant of Ca release from JSR due to overload (ms^-1)
	G_rel_bar_OL = 0;
	CSQN_bar = 10; // Max. [Ca] buffered in CSQN (mM)
	Km_CSQN = 0.8; // Equilibrium constant of buffering for CSQN (mM)
	TRPN_bar = 0.070;						
	CMDN_bar = 0.050;						
	Km_TRPN = 0.0005;
	Km_CMDN = 0.00238;
	sigma = (exp(Na_out/67.3)-1)/7;	

	// Set GUINEA_PIG array to corresponding initial values
	//published initial conditions
	/* 
	   GUINEA_PIG[0]= 0.000003; 
	   GUINEA_PIG[1]= 0.999745;
	   GUINEA_PIG[2] = 0.004503; 
	   GUINEA_PIG[3]= 0.004503;
	   GUINEA_PIG[4]= 0.000129; 
	   GUINEA_PIG[5]= 0.000994; 
	   GUINEA_PIG[6]= 0.994041;
	   GUINEA_PIG[7]= 0.000838; 
	   GUINEA_PIG[8]= 0.993336; 
	   GUINEA_PIG[9]= 0.995484; 
	   GUINEA_PIG[10]= 12.236437; //10;
	   GUINEA_PIG[11]= 136.89149; //145;
	   GUINEA_PIG[12]= 0.000079; //0.12E-3;
	   GUINEA_PIG[13]=1.179991;
	   GUINEA_PIG[14]=1.179991;
	*/

	//Initial conditions after 300 WT HH beats at BCL = 1000ms
	/* 
	   GUINEA_PIG[0]= 0.000004;
	   GUINEA_PIG[1]= 0.999723;
	   GUINEA_PIG[2]= 0.005301;
	   GUINEA_PIG[3]= 0.033794;
	   GUINEA_PIG[4]= 0.000134;
	   GUINEA_PIG[5]= 0.001020;
	   GUINEA_PIG[6]= 0.993579;
	   GUINEA_PIG[7]= 0.000878;
	   GUINEA_PIG[8]= 0.992906;
	   GUINEA_PIG[9]= 0.995226;
	   GUINEA_PIG[10]= 13.454714;
	   GUINEA_PIG[11]= 135.525797;
	   GUINEA_PIG[12]= 0.000122;
	   GUINEA_PIG[13]= 1.531388;
	   GUINEA_PIG[14]=1.942522;
	*/

	V = -83.4366;
	GUINEA_PIG[0] = 7.71003e-06;
	GUINEA_PIG[1] = 0.999098;
	GUINEA_PIG[2] = 0.0184959;
	GUINEA_PIG[3] = 0.082228;
	GUINEA_PIG[4] = 0.000261177;
	GUINEA_PIG[5] = 0.00160623;
	GUINEA_PIG[6] = 0.980753;
	GUINEA_PIG[7] = 0.00199501;
	GUINEA_PIG[8] = 0.978777;
	GUINEA_PIG[9] = 0.986735;
	GUINEA_PIG[10] = 11.753;
	GUINEA_PIG[11] = 133.987;
	GUINEA_PIG[12] = 0.000147466;
	GUINEA_PIG[13] = 1.46651;
	GUINEA_PIG[14] = 2.31135;

	// Lookup tables to speed up code for voltage dependent calculations (Range: -1000 to 1000 mV)
	lkup= new double[20000][22];
	double Vx = 0;
	V_min = -1000; 
	int z = 0;
	for (z=0; z<20000; z++) {
		Vx = V_min+0.1*z;
		lkup[z][0] = V_min+0.1*z; //V
		lkup[z][1] = 0.32*(Vx+47.13)/(1-exp(-0.1*(Vx+47.13)));//am
		lkup[z][2] = 0.08*exp(-Vx/11); //bm

		if (Vx < -40) {
			lkup[z][3] = 0.135*exp((80+Vx)/-6.8);//ah
			lkup[z][4] = 3.56*exp(0.079*Vx)+310000*exp(0.35*Vx);//bh
			lkup[z][5] = (-127140*exp(0.2444*Vx)-0.00003474*exp(-0.04391*Vx))*((Vx+37.78)/(1+exp(0.311*(Vx+79.23))));//aj
			lkup[z][6] = (0.1212*exp(-0.01052*Vx))/(1+exp(-0.1378*(Vx+40.14)));//bj
		}
		else { 
			lkup[z][3] = 0; //ah
			lkup[z][4] = 1/(0.13*(1+exp((Vx+10.66)/-11.1))); //bh
			lkup[z][5] = 0; //aj
			lkup[z][6] = (0.3*exp(-0.0000002535*Vx))/(1+exp(-0.1*(Vx+32))); //bj
		}
		lkup[z][7] = 1/(1+exp(-(Vx+10)/6.24)); //d_inf

		if (Vx<-9.99 && Vx>-10.01) lkup[z][8] = 2.289;
		else lkup[z][8] = lkup[z][7] * (1-exp(-(Vx+10)/6.24))/(0.035*(Vx+10)); //tau_d
		lkup[z][9] = 1/(1+exp((Vx+32)/8)) + 0.6/(1+exp((50-Vx)/20)); //f_inf
		lkup[z][10] = 1/(0.0197*exp(-pow((0.0337 * (Vx+10)),2)) + 0.02); //tau_f
		lkup[z][11] = 3.7+6.1/(1+exp((Vx+25.0)/4.5)); //tau_b
		if (Vx<=0) lkup[z][12] = -0.875*Vx + 12.0; //tau_g
		else lkup[z][12] = 12.0; //tau_g
		lkup[z][13] = 1/(1+exp(-(Vx+14.0)/10.8)); //b_ss
		lkup[z][14] = 1/(1+exp((Vx+60.0)/5.6)); //g_ss
		lkup[z][15] = 1/(1+exp(-(Vx+21.5)/7.5)); //xr_ss

		if (Vx<-38.89&&Vx>-39.91) lkup[z][16] = 169;
		else lkup[z][16] = 1/(0.00138*(Vx+14.2)/(1-exp(-0.123*(Vx+14.2)))+0.00061*(Vx+38.9)/(exp(0.145*(Vx+38.9))-1)); //tau_xr

		lkup[z][17] = 1/(1+exp((Vx+9)/22.4)); //r_Kr
		lkup[z][18] = 1/(1+exp(-(Vx-1.5)/16.7)); //xs1_ss

		if(Vx<-29.99&&Vx>-30.01) lkup[z][19] = 417.9;
		else lkup[z][19] = 1/(0.0000719*(Vx+30)/(1-exp(-0.148*(Vx+30)))+0.000131*(Vx+30)/(exp(0.0687*(Vx+30))-1)); //tau_xs1
		lkup[z][20] = 1/(1+exp((7.488-Vx)/5.98));	//Kp
		lkup[z][21] = 1/(1+0.1245*exp((-0.1*Vx*F)/(R*T))+0.0365*sigma*exp((-Vx*F)/(R*T))); //f_Na_K
	}

	// Cosin and Acosine lookup tables
	coslkup = new double[3141592][2];

	double num = 0;
	cos_num_min = 0.000001;
	int z1 = 0;

	for(z1 =0;z1<3141592;z1++) {
		num = cos_num_min+.000001*z1;
		coslkup[z1][0] = cos_num_min+.000001*z1;
		coslkup[z1][1] = cos(num);
	}

	acoslkup = new double[2000000][2];

	num = 0;
	acos_num_min = -1.000000;
	z1 = 0;

	for(z1 =0;z1<2000000;z1++) {
		num = acos_num_min+.000001*z1;
		acoslkup[z1][0] = acos_num_min+.000001*z1;
		acoslkup[z1][1] = acos(num);
	}

	double countX = 0;
	double px;

	for(double p = 10; p < 20; p += .1) {
		px = p/3.141592;
		px = p - (static_cast<int>(floor(px))*3.141592);
		ilow = fabs(((px)-cos_num_min)/0.000001);
		linext = -(-(px)+coslkup[ilow+1][0])/0.000001;
		cosHolder = (coslkup[ilow+1][1]-coslkup[ilow][1])*linext+coslkup[ilow+1][1];
		countX++;
	}

	modelRate = 100000; // Default model rate is 100khz

	// Set States and Parameters, connects GUI with actual variables
	setState("INa", I_Na);
	setParameter("Model Rate (Hz)", modelRate);
	setParameter("gNa", gna);

	update(MODIFY); // Update user input parameters and rate dependent variables
}

/*** Current and gating solver ***/
void Faber_Rudy_2000::solve() {
	d = GUINEA_PIG[0]; 
	f = GUINEA_PIG[1];
	xs1 = GUINEA_PIG[2]; 
	xs2 = GUINEA_PIG[3];
	xr = GUINEA_PIG[4]; 
	b = GUINEA_PIG[5]; 
	g = GUINEA_PIG[6];
	m = GUINEA_PIG[7]; 
	hh = GUINEA_PIG[8]; 
	j = GUINEA_PIG[9]; 
	Na_in = GUINEA_PIG[10];	//10;
	K_in = GUINEA_PIG[11]; //145;
	Ca_in = GUINEA_PIG[12];	//0.12E-3;
	Ca_jsr= GUINEA_PIG[13];
	Ca_nsr= GUINEA_PIG[14];

	//obtain values from lookup tables
	ilow = fabs((V-V_min)/0.1);
	linext = (-V+lkup[ilow][0])/0.1;
	am = (lkup[ilow+1][1]-lkup[ilow][1])*linext+lkup[ilow+1][1];
	bm = (lkup[ilow+1][2]-lkup[ilow][2])*linext+lkup[ilow+1][2];
	ah = (lkup[ilow+1][3]-lkup[ilow][3])*linext+lkup[ilow+1][3];  
	bh = (lkup[ilow+1][4]-lkup[ilow][4])*linext+lkup[ilow+1][4];
	aj = (lkup[ilow+1][5]-lkup[ilow][5])*linext+lkup[ilow+1][5];
	bj = (lkup[ilow+1][6]-lkup[ilow][6])*linext+lkup[ilow+1][6];
	d_inf = (lkup[ilow+1][7]-lkup[ilow][7])*linext+lkup[ilow+1][7];	
	tau_d = (lkup[ilow+1][8]-lkup[ilow][8])*linext+lkup[ilow+1][8];
	f_inf = (lkup[ilow+1][9]-lkup[ilow][9])*linext+lkup[ilow+1][9];
	tau_f = (lkup[ilow+1][10]-lkup[ilow][10])*linext+lkup[ilow+1][10];
	tau_b = (lkup[ilow+1][11]-lkup[ilow][11])*linext+lkup[ilow+1][11];
	tau_g = (lkup[ilow+1][12]-lkup[ilow][12])*linext+lkup[ilow+1][12];
	b_ss = (lkup[ilow+1][13]-lkup[ilow][13])*linext+lkup[ilow+1][13];
	g_ss = (lkup[ilow+1][14]-lkup[ilow][14])*linext+lkup[ilow+1][14];
	xr_ss = (lkup[ilow+1][15]-lkup[ilow][15])*linext+lkup[ilow+1][15];
	tau_xr = (lkup[ilow+1][16]-lkup[ilow][16])*linext+lkup[ilow+1][16];
	r_Kr = (lkup[ilow+1][17]-lkup[ilow][17])*linext+lkup[ilow+1][17];
	xs1_ss = (lkup[ilow+1][18]-lkup[ilow][18])*linext+lkup[ilow+1][18];
	tau_xs1 = (lkup[ilow+1][19]-lkup[ilow][19])*linext+lkup[ilow+1][19];
	Kp = (lkup[ilow+1][20]-lkup[ilow][20])*linext+lkup[ilow+1][20];
	f_Na_K = (lkup[ilow+1][21]-lkup[ilow][21])*linext+lkup[ilow+1][21];

	//INa
	E_Na = (((R*T)/F)*log(Na_out/Na_in));
	mtau = 1/(am+bm);
	htau = 1/(ah+bh);
	jtau = 1/(aj+bj);
	mss = am*mtau;
	hss = ah*htau;
	jss = aj*jtau;
	m = mss-(mss-m)*fastEXP(-DT/mtau);
	hh = hss-(hss-hh)*fastEXP(-DT/htau);
	j = jss-(jss-j)*fastEXP(-DT/jtau);
	I_Na = gna*m*m*m*hh*j*(V-E_Na);

	//I_Ca_L
	ad = d_inf/ tau_d;
	bd = (1-d_inf)/ tau_d;
	d = d + ((ad*(1-d))-(bd*d)) * DT;

	af = f_inf/ tau_f;
	bf = (1-f_inf)/ tau_f;
	f = f + ((af*(1-f))-(bf*f)) * DT;
	f_Ca = 1/(1+ (Ca_in/Km_Ca));
	//f_Ca = 1/(1+ pow((Ca_in/Km_Ca),1));
	I_Ca_bar = P_Ca * (z_Ca*z_Ca) * ((V*F*F)/(R*T)) * (g_Ca_in * Ca_in * fastEXP((z_Ca*V*F)/(R*T)) - g_Ca_out * Ca_out)/(fastEXP((z_Ca*V*F)/(R*T))-1);
	I_Ca_L = d * f * f_Ca * I_Ca_bar; 

	//I_Ca_K
	I_Ca_K_bar = P_K * (z_K*z_K) * ((V*F*F)/(R*T)) * (g_K_in * K_in * fastEXP((z_K*V*F)/(R*T)) - g_K_out * K_out)/(fastEXP((z_K*V*F)/(R*T))-1);
	I_Ca_K = d * f * f_Ca * I_Ca_K_bar;

	//I_Ca_Na
	I_Ca_Na_bar = P_Na * (z_Na*z_Na) * ((V*F*F)/(R*T)) * (g_Na_in * Na_in * fastEXP((z_Na*V*F)/(R*T)) - g_Na_out * Na_out)/(fastEXP((z_Na*V*F)/(R*T))-1);
	I_Ca_Na = d * f * f_Ca * I_Ca_Na_bar;

	//I_Ca_T
	E_Ca_T = (R*T/(2*F))*log(Ca_out/Ca_in);
	b = b_ss-(b_ss-b)*fastEXP(-DT/tau_b); 
	g = g_ss-(g_ss-g)*fastEXP(-DT/tau_g); 
	I_Ca_T = G_Ca_T*b*b*g*(V-E_Ca_T);

	//I_Kr
	G_Kr = 0.02614*sqrt(K_out/5.4);
	E_Kr = ((R*T)/F)*log(K_out/K_in);
	xr = xr_ss - (xr_ss - xr)*fastEXP(-DT/tau_xr);
	I_Kr = G_Kr*xr*r_Kr*(V-E_Kr);

	//I_Ks
	E_Ks = ((R*T)/F)*log((K_out+PR_NaK*Na_out)/(K_in+PR_NaK*Na_in));
	G_Ks = 0.125*(1+0.6/(1+fastEXP(1.4*log(.000038/Ca_in))));    
	//G_Ks = 0.125*(1+0.6/(1+pow((0.000038/Ca_in),1.4))); //Scaling constant for epicardial is 0.433, midmyocardial = 0.125, endocardial = 0.289
	xs2_ss = xs1_ss;
	tau_xs2 = 4*tau_xs1;
	xs1 = xs1_ss-(xs1_ss-xs1)*fastEXP(-DT/tau_xs1);
	xs2 = xs2_ss-(xs2_ss-xs2)*fastEXP(-DT/tau_xs2);
	I_Ks = G_Ks*xs1*xs2*(V-E_Ks);

	//I_K1
	G_K1=0.75*sqrt(K_out/5.4);
	E_K1=(R*T/F)*log(K_out/K_in);
	aK1=1.02/(1+fastEXP(0.2385*(V-E_K1-59.215)));
	bK1=(0.49124*fastEXP(0.08032*(V-E_K1+5.476))+fastEXP(0.06175*(V-E_K1-594.31)))/(1+fastEXP(-0.5143*(V-E_K1+4.753)));
	K1_s=aK1/(aK1+bK1);
	I_K1=G_K1*K1_s*(V-E_K1);

	//At I_K1 only fastEXP in, APD increased to 211 ms

	//I_Kp
	E_Kp=(R*T/F)*log(K_out/K_in);	//Note: E_Kp = E_K1	
	I_Kp=G_Kp*Kp*(V-E_Kp);		

	//I_Na_Ca
	I_Na_Ca = c1*fastEXP((gammas-1)*V*F/(R*T))*((fastEXP(V*F/(R*T))*Na_in*Na_in*Na_in*Ca_out-Na_out*Na_out*Na_out*Ca_in) /
												(1+c2*fastEXP((gammas-1)*V*F/(R*T))*(fastEXP(V*F/(R*T))*Na_in*Na_in*Na_in*Ca_out+Na_out*Na_out*Na_out*Ca_in)));
	//c1*exp((gammas-1)*v*F/(R*T))*((exp(v*F/(R*T))*pow(na_in, 3)*ca_out-pow(na_out, 3)*ca_in)/(1+c2*exp((gammas-1)*v*F/(R*T))*(exp(v*F/(R*T))*pow(na_in, 3)*ca_out+pow(na_out, 3)*ca_in)));

	//I_Na_K
	I_Na_K=I_Na_K_bar*f_Na_K*(1/(1+((Km_Na_i/Na_in)*(Km_Na_i/Na_in))))*(K_out/(K_out+Km_Ko));
	//I_Na_K=I_Na_K_bar*f_Na_K*(1/(1+pow((Km_Na_i/Na_in),2)))*(K_out/(K_out+Km_Ko));

	//I_p_Ca
	I_p_Ca=I_p_Ca_bar*(Ca_in/(Km_p_Ca + Ca_in));

	//I_Ca_b
	E_Ca_N = ((R*T)/(2*F))*log(Ca_out/Ca_in);
	I_Ca_b = G_Ca_b_bar*(V-E_Ca_N);

	//I_Na_b
	E_Na_N=(R*T/F)*log(Na_out/Na_in);			
	I_Na_b = G_Na_b_bar*(V-E_Na_N);

	//Total current through membrane for each ion	
	I_Na_ion_total = I_Na + I_Na_b + 3*I_Na_K + 3*I_Na_Ca +I_Ca_Na;	
	I_K_ion_total = I_Kr + I_Ks + I_K1 + I_Kp - 2*I_Na_K + I_Ca_K;		
	I_Ca_ion_total = I_Ca_L +I_Ca_T + I_Ca_b + I_p_Ca - 2*I_Na_Ca;						
	I_Total = I_Na_ion_total + I_K_ion_total + I_Ca_ion_total + I_Inject;

	if ((-1*I_Total) > 10 && tcicr > 10 && flag == 1) {flag = 0;}

	//Na_in
	dNa_in = -1*DT*(I_Na_ion_total*A_cap)/(V_myo*z_Na*F);
	Na_in = Na_in + dNa_in;

	//K_in
	dK_in = -DT*((I_K_ion_total)*A_cap)/(V_myo*z_K*F);
	K_in = K_in + dK_in;

	//Current related to the CICR
	//I_tr
	I_tr = (Ca_nsr - Ca_jsr)/tau_tr;

	//I_leak
	I_leak = K_leak*Ca_nsr; 

	//I_up
	I_up = I_up_bar*(Ca_in/(Ca_in+Km_up));

	//I_rel
	dCa_ion_total_new = (I_Ca_ion_total - Ca_in_total_old)/DT; 
	if (V>-35 && dCa_ion_total_new > dCa_ion_total && flag==0) {flag = 1; tcicr = 0;} 

	on=1/(1+fastEXP((-tcicr+4)/tau_on));
	off=(1-1/(1+fastEXP((-tcicr+4)/tau_off)));
	magrel=1/(1+fastEXP(((I_Ca_ion_total)+5)/0.9));   //I_Ca_L -2*I_Na_Ca + I_p_Ca + I_Ca_T + I_Ca_b
	I_rel = G_rel_bar*on*off*magrel*(Ca_jsr-Ca_in); 
	tcicr=tcicr+DT;		

	//I_rel_OL
	G_rel_OL = G_rel_bar_OL*(1-fastEXP(-tjsrol/tau_on))*fastEXP(-tjsrol/tau_off);
	I_rel_OL = G_rel_OL*(Ca_jsr-Ca_in);
	tjsrol=tjsrol+DT; 

	//Ca_jsr
	Ca_CSQN=CSQN_bar*(Ca_jsr/(Ca_jsr+Km_CSQN));
	dCa_jsr=DT*(I_tr-I_rel-I_rel_OL);
	bjsr=CSQN_bar - Ca_CSQN - dCa_jsr - Ca_jsr+Km_CSQN;
	cjsr=Km_CSQN*(Ca_CSQN + dCa_jsr+Ca_jsr);
	Ca_jsr=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2; 

	//Ca_nsr
	dCa_nsr = DT*(I_up - I_leak - I_tr*V_jsr/V_nsr);
	Ca_nsr=Ca_nsr+dCa_nsr; 

	//Ca_in
	dCa_in = -DT*(((I_Ca_ion_total*A_cap)/(V_myo*z_Ca*F))+((I_up - I_leak)*V_nsr/V_myo)-(I_rel*V_jsr/V_myo)-(I_rel_OL*V_jsr/V_myo));
	TRPN = TRPN_bar*(Ca_in/(Ca_in + Km_TRPN));
	CMDN = CMDN_bar*( Ca_in/(Ca_in + Km_CMDN));
	Ca_total = TRPN + CMDN + dCa_in + Ca_in;
	b_myo = CMDN_bar + TRPN_bar - Ca_total + Km_TRPN + Km_CMDN;
	c_myo = (Km_CMDN*Km_TRPN)-(Ca_total*(Km_TRPN + Km_CMDN))+(TRPN_bar*Km_CMDN)+(CMDN_bar*Km_TRPN);
	d_myo = -Km_TRPN*Km_CMDN*Ca_total;
	g_myo = sqrt(b_myo*b_myo-3*c_myo);

	//Ca_in = (2*g_myo/3)*cos(acos(9*b_myo*c_myo-2*b_myo*b_myo*b_myo-27*d_myo)/(2*((b_myo*b_myo-3*c_myo)*sqrt(b_myo*b_myo-3*c_myo)))/3)-(b_myo/3);  //Original


	//Ca_in = (2*g_myo/3)*cos(acos((9*b_myo*c_myo-2*b_myo*b_myo*b_myo-27*d_myo)/(2*((b_myo*b_myo-3*c_myo)*sqrt(b_myo*b_myo-3*c_myo))))/3)-(b_myo/3);  //Corina change


	// Optimized Ca_in calculation
	Ca_in_x = (9*b_myo*c_myo-2*b_myo*b_myo*b_myo-27*d_myo)/(2*((b_myo*b_myo-3*c_myo)*sqrt(b_myo*b_myo-3*c_myo)));

	ilow = fabs((Ca_in_x-acos_num_min)/0.000001);
	linext = -(-Ca_in_x+acoslkup[ilow+1][0])/0.000001;
	acosHolder = (acoslkup[ilow+1][1]-acoslkup[ilow][1])*linext+acoslkup[ilow+1][1];
	acosHolder = acosHolder / 3;
	ilow = fabs(((acosHolder)-cos_num_min)/0.000001);
	linext = -(-(acosHolder)+coslkup[ilow+1][0])/0.000001;
	cosHolder = (coslkup[ilow+1][1]-coslkup[ilow][1])*linext+coslkup[ilow+1][1];
	cosHolder2 = cosHolder / 3.141592;
	cosHolder = cosHolder - (static_cast<int>(floor(cosHolder2))*3.141592);
	ilow = fabs(((acosHolder)-cos_num_min)/0.000001);
	Ca_in = (2*g_myo/3) * cosHolder - (b_myo/3);


	Ca_in_total_old = I_Ca_ion_total;
	dCa_ion_total = dCa_ion_total_new; 

	GUINEA_PIG[0] = d; 
	GUINEA_PIG[1]= f;
	GUINEA_PIG[2]= xs1; 
	GUINEA_PIG[3]= xs2;
	GUINEA_PIG[4]= xr; 
	GUINEA_PIG[5]= b; 
	GUINEA_PIG[6]= g;
	GUINEA_PIG[7]= m; 
	GUINEA_PIG[8]= hh; 
	GUINEA_PIG[9]= j; 
	GUINEA_PIG[10]= Na_in;	
	GUINEA_PIG[11]= K_in;
	GUINEA_PIG[12]= Ca_in;	
	GUINEA_PIG[13]= Ca_jsr;
	GUINEA_PIG[14]= Ca_nsr;
	V = V - I_Total*DT;
} 
