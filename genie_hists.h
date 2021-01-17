//initialisation of histograms
h1_el_Mott_crosssec = new TH1F("h1_el_Mott_crosssec","",200,0.,0.01);
h1_Wvar = new TH1F("h1_Wvar","",400,0,3);
h1_xbjk = new TH1F("h1_xbjk","",400,0,3);
h1_Q2 = new TH1F("h1_Q2","",400,0,6);
h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
h1_Nprot=new TH1F("h1_Nprot","",10,-0.5,4.5);
h1_Nprot_NonZeroProt=new TH1F("h1_Nprot_NonZeroProt","",8,0.5,4.5);
h1_Nphot=new TH1F("h1_Nphot","",10,-0.5,4.5);
h1_Npiphot=new TH1F("h1_Npiphot","",10,-0.5,4.5);
h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,-0.5,4.5);
h1_Npi=new TH1F("h1_Npi","",10,-0.5,4.5);
h1_Npi_NonZeroProt=new TH1F("h1_Npi_NonZeroProt","",10,-0.5,4.5);
h1_Npipl=new TH1F("h1_Npipl","",10,-0.5,4.5);
h1_Npimi=new TH1F("h1_Npimi","",10,-0.5,4.5);
h1_MissMomentum = new TH1F("MissMomentum","",80,0.,1.);
h1_el_mom = new TH1F("h1_el_mom","",100,0.2,6);
h1_el_mom_corr = new TH1F("h1_el_mom_corr","",100,0.,5.);
h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
h1_prot_mom = new TH1F("h1_prot_mom","",300,0,3);
h1_prot_mom_ratio = new TH1F("h1_prot_mom_ratio","",50,0.97,1.2);
h1_Wvar_weight = new TH1F("h1_Wvar_weight","",400,0,3);
h1_xbjk_weight = new TH1F("h1_xbjk_weight","",400,0,3);
h1_Q2_weight = new TH1F("h1_Q2_weight","",400,0,6);
h1_nu_weight = new TH1F("h1_nu_weight","",400,0,4);
h1_WvarCal_weight = new TH1F("h1_WvarCal_weight","",400,0,3);
h1_xbjkCal_weight = new TH1F("h1_xbjkCal_weight","",400,0,3);
h1_Q2Cal_weight = new TH1F("h1_Q2Cal_weight","",400,0,6);
h1_nuCal_weight = new TH1F("h1_nuCal_weight","",400,0,3);


	//Binning for fractional feed down energy histograms
	int N_qe = -1;
	double *x_qe;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
		N_qe=109;
		x_qe=new double[N_qe+1];
		for (int i=0;i<=64;i++)	x_qe[i]=-1+i*0.015;
		for (int i=0;i<=44;i++)	x_qe[i+65]=-0.04+(i+1)*0.01;
	}

	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
		N_qe=109;
		x_qe=new double[N_qe+1];
		for (int i=0;i<=64;i++)	x_qe[i]=-1+i*0.015;
		for (int i=0;i<=44;i++)	x_qe[i+65]=-0.04+(i+1)*0.01;
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
		N_qe=82;
		x_qe=new double[N_qe+1];
		for (int i=0;i<=29;i++)	x_qe[i]=-1+i*0.03;
		for (int i=0;i<=52;i++)	x_qe[i+30]=-0.13+(i+1)*0.01;
	}

	//Definitions of further Histograms
	CosDeltaThetaElectronPhotonAboveThreshold=new TH1F("CosDeltaThetaElectronPhotonAboveThreshold","",100,-1.,1.);
	CosDeltaPhiElectronPhotonAboveThreshold=new TH1F("CosDeltaPhiElectronPhotonAboveThreshold","",100,-1.,1.);

	RadCosThetaGammaEgamma = new TH2F("RadCosThetaGammaEgamma","",100,-1.,1.,600,0.,6.);
	RadCosDeltaThetaGammaEgamma = new TH2F("RadCosDeltaThetaGammaEgamma","",100,-1.,1.,600,0.,6.);
	NonRadThetaVsPhiGamma = new TH2F("NonRadThetaVsPhiGamma","",360,0.,360.,180,0.,180.);

	h1_E_rec_1pi_weight_frac_feed=new TH1F("h1_E_rec_1pi_weight_frac_feed","",N_qe,x_qe);
	h1_E_rec_2pi_weight_frac_feed=new TH1F("h1_E_rec_2pi_weight_frac_feed","",N_qe,x_qe);
	h1_E_rec_3pi_weight_frac_feed=new TH1F("h1_E_rec_3pi_weight_frac_feed","",N_qe,x_qe);
	h1_E_rec_4pi_weight_frac_feed=new TH1F("h1_E_rec_4pi_weight_frac_feed","",N_qe,x_qe);
	h1_E_rec_0pi_frac_feed=new TH1F("h1_E_rec_0pi_frac_feed","",N_qe,x_qe);
	h1_E_tot_cut2_fracfeed = new TH1F("h1_E_tot_cut2_fracfeed","",N_qe,x_qe);
	h1_E_rec_cut2_new_fracfeed = new TH1F("h1_E_rec_cut2_new_fracfeed","",N_qe,x_qe);
	h1_E_tot_p_bkgd_fracfeed = new TH1F("h1_E_tot_p_bkgd_fracfeed","",N_qe,x_qe);
	h1_E_rec_p_bkgd_fracfeed = new TH1F("h1_E_rec_p_bkgd_fracfeed","",N_qe,x_qe);
	h1_E_tot_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_2p1pi_2p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_2p0pi_fracfeed","",N_qe,x_qe);
	h1_E_tot_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_2p1pi_1p1pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p1pi_fracfeed","",N_qe,x_qe);
	h1_E_tot_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_tot_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_2p1pi_1p0pi_fracfeed = new TH1F("h1_E_rec_2p1pi_1p0pi_fracfeed","",N_qe,x_qe);
	h1_E_tot_3pto2p_fracfeed = new TH1F("h1_E_tot_3pto2p_fracfeed","",N_qe,x_qe);
	h1_E_rec_3pto2p_fracfeed = new TH1F("h1_E_rec_3pto2p_fracfeed","",N_qe,x_qe);
	h1_E_tot_3pto1p_fracfeed = new TH1F("h1_E_tot_3pto1p_fracfeed","",N_qe,x_qe);
	h1_E_rec_3pto1p_fracfeed = new TH1F("h1_E_rec_3pto1p_fracfeed","",N_qe,x_qe);
	h1_E_tot_4pto3p_fracfeed = new TH1F("h1_E_tot_4pto3p_fracfeed","",N_qe,x_qe);
	h1_E_rec_4pto3p_fracfeed = new TH1F("h1_E_rec_4pto3p_fracfeed","",N_qe,x_qe);
	h1_E_tot_43pto1p_fracfeed = new TH1F("h1_E_tot_43pto1p_fracfeed","",N_qe,x_qe);
	h1_E_rec_43pto1p_fracfeed = new TH1F("h1_E_rec_43pto1p_fracfeed","",N_qe,x_qe);
	h1_E_tot_4pto2p_fracfeed = new TH1F("h1_E_tot_4pto2p_fracfeed","",N_qe,x_qe);
	h1_E_rec_4pto2p_fracfeed = new TH1F("h1_E_rec_4pto2p_fracfeed","",N_qe,x_qe);
	h1_E_tot_4pto1p_fracfeed = new TH1F("h1_E_tot_4pto1p_fracfeed","",N_qe,x_qe);
	h1_E_rec_4pto1p_fracfeed = new TH1F("h1_E_rec_4pto1p_fracfeed","",N_qe,x_qe);
	h1_E_tot_1p2pi_fracfeed = new TH1F("h1_E_tot_1p2pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_1p2pi_fracfeed = new TH1F("h1_E_rec_1p2pi_fracfeed","",N_qe,x_qe);
	h1_E_tot_1p3pi_fracfeed = new TH1F("h1_E_tot_1p3pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_1p3pi_fracfeed = new TH1F("h1_E_rec_1p3pi_fracfeed","",N_qe,x_qe);
	h1_E_tot_2p2pi_fracfeed = new TH1F("h1_E_tot_2p2pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_2p2pi_fracfeed = new TH1F("h1_E_rec_2p2pi_fracfeed","",N_qe,x_qe);
	h1_E_tot_3p1pi_fracfeed = new TH1F("h1_E_tot_3p1pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_3p1pi_fracfeed = new TH1F("h1_E_rec_3p1pi_fracfeed","",N_qe,x_qe);
	h1_E_tot_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_tot_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_1p2pi_1p0pi_fracfeed = new TH1F("h1_E_rec_1p2pi_1p0pi_fracfeed","",N_qe,x_qe);
	h1_E_rec_undetfactor_fracfeed = new TH1F("h1_E_rec_undetfactor_fracfeed","",N_qe,x_qe);
	h1_E_tot_undetfactor_fracfeed = new TH1F("h1_E_tot_undetfactor_fracfeed","",N_qe,x_qe);

	h1_theta0=new TH1F("h1_theta0","",300,0,180);
	h2_Ecal_Eqe=new TH2F("h2_Ecal_Eqe","",800,0,8.,800,0,8.);
	h2_EqeEcalratio_Eqe=new TH2F("h2_EqeEcalratio_Eqe","",600,0,5,300,0,2);
	h2_EqeEcaldiff_Eqe=new TH2F("h2_EqeEcaldiff_Eqe","",600,0,5,300,-3,3);
	h2_N_prot_pi=new TH2F("h2_N_prot_pi","",10,0,5,10,0,5);
	h2_N_prot_pi_phot=new TH2F("h2_N_prot_pi_phot","",10,0,5,10,0,5);
	h2_N_prot_pi_phot_nonrad=new TH2F("h2_N_prot_pi_phot_nonrad","",10,0,5,10,0,5);
//	h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,0,180);
	h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,10,60);
	h2_el_mom_diff = new TH2F("h2_el_mom_diff","",500,0.,1.,500,-0.1,0.1);

	int NBinsNu = 300, NBinsQ2 = 300;
	double MinNu = 0., MaxNu = 4.; double MinQ2 = 0., MaxQ2 = 6.;
	h2_Q2_nu = new TH2F("h2_Q2_nu","",NBinsNu,MinNu,MaxNu,NBinsQ2,MinQ2,MaxQ2);
	h2_Q2_nu_weight = new TH2F("h2_Q2_nu_weight","",NBinsNu,MinNu,MaxNu,NBinsQ2,MinQ2,MaxQ2);
	h2_Q2_nu_weight_FirstSector = new TH2F("h2_Q2_nu_weight_FirstSector","",0.7*NBinsNu,MinNu,MaxNu,0.7*NBinsQ2,MinQ2,MaxQ2);

	h2_Q2_xbjk_weight = new TH2F("h2_Q2_xbjk_weight","",200,0,3,200,0,5);
	h2_Q2_W=new TH2F("h2_Q2_W","",200,0,3,200,0,5);
	h2_xB_W=new TH2F("h2_xB_W","",200,0,3,200,0,3);
	h2_Q2_W_weight=new TH2F("h2_Q2_W_weight","",200,0,3,200,0,5);
	h2_el_pcorr_puncorr = new TH2F("h2_el_pcorr_puncorr","",100,0,1,100,0,3);
	h2_Erec_pperp = new TH2F("h2_Erec_pperp","",400,0,1,400,0,6.);
	h2_Erec_pperp_newcut2 = new TH2F("h2_Erec_pperp_newcut2","",400,0,1,400,0,6.);
	h2_Erec_pperp_cut3 = new TH2F("h2_Erec_pperp_cut3","",400,0,1,400,0,6.);
	h2_Erec_pperp_2p = new TH2F("h2_Erec_pperp_2p","",400,0,1,400,0,6.);
	h2_Erec_pperp_321p = new TH2F("h2_Erec_pperp_321p","",400,0,1,400,0,6.);
	h2_Erec_pperp_31p = new TH2F("h2_Erec_pperp_31p","",400,0,1,400,0,6.);
	h2_Erec_pperp_4321p = new TH2F("h2_Erec_pperp_4321p","",400,0,1,400,0,6.);
	h2_Erec_pperp_431p = new TH2F("h2_Erec_pperp_431p","",400,0,1,400,0,6.);
	h2_Erec_pperp_421p = new TH2F("h2_Erec_pperp_421p","",400,0,1,400,0,6.);
	h2_Erec_pperp_41p = new TH2F("h2_Erec_pperp_41p","",400,0,1,400,0,6.);
	h2_Erec_pperp_1p1pi = new TH2F("h2_Erec_pperp_1p1pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_1p2pi_1p0pi = new TH2F("h2_Erec_pperp_1p2pi_1p0pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_1p2pi_1p1pi = new TH2F("h2_Erec_pperp_1p2pi_1p1pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_2p1pi_2p0pi = new TH2F("h2_Erec_pperp_2p1pi_2p0pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_2p1pi_1p1pi = new TH2F("h2_Erec_pperp_2p1pi_1p1pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_2p1pi_1p0pi = new TH2F("h2_Erec_pperp_2p1pi_1p0pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_1p3pi = new TH2F("h2_Erec_pperp_1p3pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_2p2pi = new TH2F("h2_Erec_pperp_2p2pi","",400,0,1,400,0,6.);
	h2_Erec_pperp_3p1pi = new TH2F("h2_Erec_pperp_3p1pi","",400,0,1,400,0,6.);
	h2_pperp_W=new TH2F("h2_pperp_W","",200,0,3,200,0,2);
	h2_Etot_pperp = new TH2F("h2_Etot_pperp","",400,0,1,400,0,6.);

	h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);

	h2_QVector_theta_phi = new TH2F("h2_QVector_theta_phi","",200,0,360,200,0,80);

	//Binning for energy reconstruction histograms
	int n_bins;
	double *x_values;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
		n_bins=38;
		x_values=new double[n_bins+1];
		for (int i=0;i<=17;i++) x_values[i]=0.4+i*0.04;
		for (int i=0;i<=20;i++) x_values[i+18]=1.08+(i+1)*0.02;
	}

	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
		n_bins=54;
		x_values=new double[n_bins+1];
		for (int i=0;i<=23;i++) x_values[i]=i*0.09;
		for (int i=0;i<=30;i++) x_values[i+24]=2.07+(i+1)*0.03;
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
		n_bins=38;
		x_values=new double[n_bins+1];
		for (int i=0;i<=21;i++)	x_values[i]=i*0.2;
		for (int i=0;i<=16;i++)	x_values[i+22]=4.2+(i+1)*0.05;
	}

	//Definitions of further Histograms
	h1_E_rec_2p_det = new TH1F("h1_E_rec_2p_det","",n_bins,x_values);
	h1_E_tot_2p_det = new TH1F("h1_E_tot_2p_det","",n_bins,x_values);
	h1_E_tot_p_bkgd = new TH1F("h1_E_tot_p_bkgd","",n_bins,x_values);
	h1_E_rec_p_bkgd = new TH1F("h1_E_rec_p_bkgd","",n_bins,x_values);
	h1_E_tot_3pto1p = new TH1F("h1_E_tot_3pto1p","",n_bins,x_values);
	h1_E_rec_3pto1p = new TH1F("h1_E_rec_3pto1p","",n_bins,x_values);
	h1_E_tot_43pto1p = new TH1F("h1_E_tot_43pto1p","",n_bins,x_values);
	h1_E_rec_43pto1p = new TH1F("h1_E_rec_43pto1p","",n_bins,x_values);
	h1_E_tot_3pto2p =new TH1F("h1_E_tot_3pto2p","",n_bins,x_values);
	h1_E_rec_3pto2p = new TH1F("h1_E_rec_3pto2p","",n_bins,x_values);
	h1_E_tot_4pto1p =new TH1F("h1_E_tot_4pto1p","",n_bins,x_values);
	h1_E_rec_4pto1p = new TH1F("h1_E_rec_4pto1p","",n_bins,x_values);
	h1_E_tot_4pto3p =new TH1F("h1_E_tot_4pto3p","",n_bins,x_values);
	h1_E_rec_4pto3p = new TH1F("h1_E_rec_4pto3p","",n_bins,x_values);
	h1_E_tot_4pto2p =new TH1F("h1_E_tot_4pto2p","",n_bins,x_values);
	h1_E_rec_4pto2p = new TH1F("h1_E_rec_4pto2p","",n_bins,x_values);
	h1_E_rec = new TH1F("h1_E_rec","",n_bins,x_values);
	h1_E_rec_0pi = new TH1F("h1_E_rec_0pi","",n_bins,x_values);
	h1_E_rec_1pi = new TH1F("h1_E_rec_1pi","",n_bins,x_values);
	h1_E_rec_1pi_weight = new TH1F("h1_E_rec_1pi_weight","",n_bins,x_values);
	h1_E_rec_2pi_weight = new TH1F("h1_E_rec_2pi_weight","",n_bins,x_values);
	h1_E_rec_3pi_weight = new TH1F("h1_E_rec_3pi_weight","",n_bins,x_values);
	h1_E_rec_4pi_weight = new TH1F("h1_E_rec_4pi_weight","",n_bins,x_values);
	h1_E_rec_20pi = new TH1F("h1_E_rec_20pi","",n_bins,x_values);
	h1_E_rec_21pi = new TH1F("h1_E_rec_21pi","",n_bins,x_values);
	h1_E_rec_30pi = new TH1F("h1_E_rec_30pi","",n_bins,x_values);
	h1_E_rec_310pi = new TH1F("h1_E_rec_310pi","",n_bins,x_values);
	h1_E_rec_320pi = new TH1F("h1_E_rec_320pi","",n_bins,x_values);
	h1_E_rec_3210pi = new TH1F("h1_E_rec_3210pi","",n_bins,x_values);
	h1_E_rec_40pi = new TH1F("h1_E_rec_40pi","",n_bins,x_values);
	h1_E_rec_410pi = new TH1F("h1_E_rec_410pi","",n_bins,x_values);
	h1_E_rec_420pi = new TH1F("h1_E_rec_420pi","",n_bins,x_values);
	h1_E_rec_4210pi = new TH1F("h1_E_rec_4210pi","",n_bins,x_values);
	h1_E_rec_430pi = new TH1F("h1_E_rec_430pi","",n_bins,x_values);
	h1_E_rec_4310pi = new TH1F("h1_E_rec_4310pi","",n_bins,x_values);
	h1_E_rec_4320pi = new TH1F("h1_E_rec_4320pi","",n_bins,x_values);
	h1_E_rec_43210pi = new TH1F("h1_E_rec_43210pi","",n_bins,x_values);
	h1_E_rec_1prot	= new TH1F("h1_E_rec_1prot","",n_bins,x_values);
	h1_E_tot_1prot	= new TH1F("h1_E_tot_1prot","",n_bins,x_values);
	h1_E_rec_cutpi1_piplpimi = new TH1F("h1_E_rec_cutpi1_piplpimi","",n_bins,x_values);
	h1_E_tot_cutpi1_piplpimi = new TH1F("h1_E_tot_cutpi1_piplpimi","",n_bins,x_values);
	h1_E_tot = new TH1F("h1_E_tot","",n_bins,x_values);
	h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
	h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
	h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
	h1_E_rec_undetfactor	= new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
	h1_E_tot_undetfactor	= new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
	h1_E_tot_1p2pi	= new TH1F("h1_E_tot_1p2pi","",n_bins,x_values);
	h1_E_rec_1p2pi	= new TH1F("h1_E_rec_1p2pi","",n_bins,x_values);
	h1_E_tot_1p3pi	= new TH1F("h1_E_tot_1p3pi","",n_bins,x_values);
	h1_E_rec_1p3pi	= new TH1F("h1_E_rec_1p3pi","",n_bins,x_values);
	h1_E_tot_2p2pi	= new TH1F("h1_E_tot_2p2pi","",n_bins,x_values);
	h1_E_rec_2p2pi	= new TH1F("h1_E_rec_2p2pi","",n_bins,x_values);
	h1_E_tot_3p1pi	= new TH1F("h1_E_tot_3p1pi","",n_bins,x_values);
	h1_E_rec_3p1pi	= new TH1F("h1_E_rec_3p1pi","",n_bins,x_values);
	h1_E_tot_1p2pi_1p0pi	= new TH1F("h1_E_tot_1p2pi_1p0pi","",n_bins,x_values);
	h1_E_rec_1p2pi_1p0pi	= new TH1F("h1_E_rec_1p2pi_1p0pi","",n_bins,x_values);
	h1_E_tot_2p1pi_2p0pi	= new TH1F("h1_E_tot_2p1pi_2p0pi","",n_bins,x_values);
	h1_E_rec_2p1pi_2p0pi	= new TH1F("h1_E_rec_2p1pi_2p0pi","",n_bins,x_values);
	h1_E_tot_2p1pi_1p1pi	= new TH1F("h1_E_tot_2p1pi_1p1pi","",n_bins,x_values);
	h1_E_rec_2p1pi_1p1pi	= new TH1F("h1_E_rec_2p1pi_1p1pi","",n_bins,x_values);
	h1_E_tot_2p1pi_1p0pi	= new TH1F("h1_E_tot_2p1pi_1p0pi","",n_bins,x_values);
	h1_E_rec_2p1pi_1p0pi	= new TH1F("h1_E_rec_2p1pi_1p0pi","",n_bins,x_values);

	// Unweighted plots for the number of events

	h1_MissMomentum_NoWeight = new TH1F("MissMomentum_NoWeight","",100,0.,1.);

	h1_ECal_Slice0_NoWeight = new TH1F("epRecoEnergy_slice_0_NoWeight","",n_bins,x_values);
	h1_ECal_Slice1_NoWeight = new TH1F("epRecoEnergy_slice_1_NoWeight","",n_bins,x_values);
	h1_ECal_Slice2_NoWeight = new TH1F("epRecoEnergy_slice_2_NoWeight","",n_bins,x_values);
	h1_ECal_Slice3_NoWeight = new TH1F("epRecoEnergy_slice_3_NoWeight","",n_bins,x_values);

	h1_EQE_Slice0_NoWeight = new TH1F("eRecoEnergy_slice_0_NoWeight","",n_bins,x_values);
	h1_EQE_Slice1_NoWeight = new TH1F("eRecoEnergy_slice_1_NoWeight","",n_bins,x_values);
	h1_EQE_Slice2_NoWeight = new TH1F("eRecoEnergy_slice_2_NoWeight","",n_bins,x_values);
	h1_EQE_Slice3_NoWeight = new TH1F("eRecoEnergy_slice_3_NoWeight","",n_bins,x_values);

	//Defintions of Histogram for each slice
	for(int h = 0; h < n_slice; h++){
		h1_Erec_p_bkgd_slice[h]= new TH1F(Form("h1_Erec_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice[h]= new TH1F(Form("h1_Etot_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_3pto1p_slice[h]= new TH1F(Form("h1_Erec_3pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_3pto1p_slice[h]= new TH1F(Form("h1_Etot_3pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_3pto2p_slice[h]= new TH1F(Form("h1_Erec_3pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_3pto2p_slice[h]= new TH1F(Form("h1_Etot_3pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_3p1pi_slice[h]= new TH1F(Form("h1_Erec_3p1pi_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_3p1pi_slice[h]= new TH1F(Form("h1_Etot_3p1pi_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_43pto1p_slice[h]= new TH1F(Form("h1_Etot_43pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_43pto1p_slice[h]= new TH1F(Form("h1_Erec_43pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_4pto3p_slice[h]= new TH1F(Form("h1_Erec_4pto3p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_4pto3p_slice[h]= new TH1F(Form("h1_Etot_4pto3p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_4pto2p_slice[h]= new TH1F(Form("h1_Erec_4pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_4pto2p_slice[h]= new TH1F(Form("h1_Etot_4pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_4pto1p_slice[h]= new TH1F(Form("h1_Erec_4pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_4pto1p_slice[h]= new TH1F(Form("h1_Etot_4pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_Npi0[h] = new TH1F(Form("h1_Etot_Npi0_%d",h+1),"",n_bins,x_values);
		h1_Erec_Npi0[h] = new TH1F(Form("h1_Erec_Npi0_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_pipl_pimi_fact[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_pipl_pimi_fact_pipl[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pipl_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_pipl_pimi_fact_pimi[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pimi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_pipl_pimi_new_fact[h]= new TH1F(Form("h1_Erec_bkgd_pipl_pimi_new_fact_%d",h+1),"",n_bins,x_values);
		h1_Etot_Npi1[h] = new TH1F(Form("h1_Etot_Npi1_%d",h+1),"",n_bins,x_values);
		h1_Erec_Npi1[h] = new TH1F(Form("h1_Erec_Npi1_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_1p2pi[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_1p2pi[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to2p0pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to2p0pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p2pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p2pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p0pi_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p0pi_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p2pi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p2pi_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_1p2pi_1p0pi[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_1p0pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_1p2pi_1p0pi[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_1p0pi_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_1p3pi[h] = new TH1F(Form("h1_Etot_bkgd_1p3pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_1p3pi[h] = new TH1F(Form("h1_Erec_bkgd_1p3pi_%d",h+1),"",n_bins,x_values);
	}

	for(int h=0;h<20;h++){
		h2_N_pi_phot[h]=new TH2F(Form("h2_N_pi_phot_%d",h),"",10,0,5,10,0,5);
	}





//1p 1pi histos
h1_en_recon1=new TH1F("calorimetric","", n_bins, x_values);
  h1_en_recon1->Sumw2();
h1_en_recon1_pimi=new TH1F("cal_pimi", "", n_bins, x_values);
  h1_en_recon1_pimi->Sumw2();
h1_en_recon1_pipl=new TH1F("cal_pipl", "",n_bins, x_values);
  h1_en_recon1_pipl->Sumw2();
h1_en_recon2=new TH1F("kin_e","",n_bins, x_values);
  h1_en_recon2->Sumw2();
h1_en_recon2_pimi=new TH1F("kin_e_pimi","",n_bins, x_values);
  h1_en_recon2_pimi->Sumw2();
h1_en_recon2_pipl=new TH1F("kin_e_pipl","",n_bins, x_values);
  h1_en_recon2_pipl->Sumw2();
h1_en_recon3=new TH1F("kin_e_pi","",n_bins, x_values);
  h1_en_recon3->Sumw2();
h1_en_recon3_pimi=new TH1F("kin_e_pi_pimi","",n_bins, x_values);
  h1_en_recon3_pimi->Sumw2();
h1_en_recon3_pipl=new TH1F("kin_e_pi_pipl","",n_bins, x_values);
  h1_en_recon3_pipl->Sumw2();

h1_en_recon1_Q2_1=new TH1F("calorimetric_q2bin","", n_bins, x_values);
  h1_en_recon1_Q2_1->Sumw2();
h1_en_recon1_pimi_Q2_1=new TH1F("cal_pimi_q2bin", "", n_bins, x_values);
  h1_en_recon1_pimi_Q2_1->Sumw2();
h1_en_recon1_pipl_Q2_1=new TH1F("cal_pipl_q2bin", "",n_bins, x_values);
  h1_en_recon1_pipl_Q2_1->Sumw2();
h1_en_recon2_Q2_1=new TH1F("kin_e_q2bin","",n_bins, x_values);
  h1_en_recon2_Q2_1->Sumw2();
h1_en_recon2_pimi_Q2_1=new TH1F("kin_e_pimi_q2bin","",n_bins, x_values);
  h1_en_recon2_pimi_Q2_1->Sumw2();
h1_en_recon2_pipl_Q2_1=new TH1F("kin_e_pipl_q2bin","",n_bins, x_values);
  h1_en_recon2_pipl_Q2_1->Sumw2();
h1_en_recon3_Q2_1=new TH1F("kin_e_pi_q2bin","",n_bins, x_values);
  h1_en_recon3_Q2_1->Sumw2();
h1_en_recon3_pimi_Q2_1=new TH1F("kin_e_pi_pimi_q2bin","",n_bins, x_values);
  h1_en_recon3_pimi_Q2_1->Sumw2();
h1_en_recon3_pipl_Q2_1=new TH1F("kin_e_pi_pipl_q2bin","",n_bins, x_values);
  h1_en_recon3_pipl_Q2_1->Sumw2();


h1_rot1_2pi_1p=new TH1F("2pi_1p_rot_cal","",n_bins, x_values);
  h1_rot1_2pi_1p->Sumw2();
h1_rot1_2pi_1p_pimi=new TH1F("2pi_1p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_2pi_1p_pimi->Sumw2();
h1_rot1_2pi_1p_pipl=new TH1F("2pi_1p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_2pi_1p_pipl->Sumw2();
h1_rot2_2pi_1p=new TH1F("2pi_1p_rot_kin_e","",n_bins, x_values);
  h1_rot2_2pi_1p->Sumw2();
h1_rot2_2pi_1p_pimi=new TH1F("2pi_1p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_2pi_1p_pimi->Sumw2();
h1_rot2_2pi_1p_pipl=new TH1F("2pi_1p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_2pi_1p_pipl->Sumw2();
h1_rot3_2pi_1p=new TH1F("2pi_1p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_2pi_1p->Sumw2();
h1_rot3_2pi_1p_pimi=new TH1F("2pi_1p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_2pi_1p_pimi->Sumw2();
h1_rot3_2pi_1p_pipl=new TH1F("2pi_1p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_2pi_1p_pipl->Sumw2();
h1_rot1_1pi_2p=new TH1F("1pi_2p_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_2p->Sumw2();
h1_rot1_1pi_2p_pimi=new TH1F("1pi_2p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_2p_pimi->Sumw2();
h1_rot1_1pi_2p_pipl=new TH1F("1pi_2p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_2p_pipl->Sumw2();
h1_rot2_1pi_2p=new TH1F("1pi_2p_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_2p->Sumw2();
h1_rot2_1pi_2p_pimi=new TH1F("1pi_2p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_2p_pimi->Sumw2();
h1_rot2_1pi_2p_pipl=new TH1F("1pi_2p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_2p_pipl->Sumw2();
h1_rot3_1pi_2p=new TH1F("1pi_2p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_2p->Sumw2();
h1_rot3_1pi_2p_pimi=new TH1F("1pi_2p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_2p_pimi->Sumw2();
h1_rot3_1pi_2p_pipl=new TH1F("1pi_2p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_2p_pipl->Sumw2();
h1_rot1_2pi_2p=new TH1F("2pi_2p_rot_cal","",n_bins, x_values);
  h1_rot1_2pi_2p->Sumw2();
h1_rot1_2pi_2p_pimi=new TH1F("2pi_2p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_2pi_2p_pimi->Sumw2();
h1_rot1_2pi_2p_pipl=new TH1F("2pi_2p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_2pi_2p_pipl->Sumw2();
h1_rot2_2pi_2p=new TH1F("2pi_2p_rot_kin_e","",n_bins, x_values);
  h1_rot2_2pi_2p->Sumw2();
h1_rot2_2pi_2p_pimi=new TH1F("2pi_2p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_2pi_2p_pimi->Sumw2();
h1_rot2_2pi_2p_pipl=new TH1F("2pi_2p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_2pi_2p_pipl->Sumw2();
h1_rot3_2pi_2p=new TH1F("2pi_2p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_2pi_2p->Sumw2();
h1_rot3_2pi_2p_pimi=new TH1F("2pi_2p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_2pi_2p_pimi->Sumw2();
h1_rot3_2pi_2p_pipl=new TH1F("2pi_2p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_2pi_2p_pipl->Sumw2();
h1_rot1_1pi_3p=new TH1F("1pi_3p_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_3p->Sumw2();
h1_rot1_1pi_3p_pimi=new TH1F("1pi_3p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_3p_pimi->Sumw2();
h1_rot1_1pi_3p_pipl=new TH1F("1pi_3p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_3p_pipl->Sumw2();
h1_rot2_1pi_3p=new TH1F("1pi_3p_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_3p->Sumw2();
h1_rot2_1pi_3p_pimi=new TH1F("1pi_3p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_3p_pimi->Sumw2();
h1_rot2_1pi_3p_pipl=new TH1F("1pi_3p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_3p_pipl->Sumw2();
h1_rot3_1pi_3p=new TH1F("1pi_3p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_3p->Sumw2();
h1_rot3_1pi_3p_pimi=new TH1F("1pi_3p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_3p_pimi->Sumw2();
h1_rot3_1pi_3p_pipl=new TH1F("1pi_3p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_3p_pipl->Sumw2();
h1_rot1_3pi_1p=new TH1F("3pi_1p_rot_cal","",n_bins, x_values);
  h1_rot1_3pi_1p->Sumw2();
h1_rot1_3pi_1p_pimi=new TH1F("3pi_1p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_3pi_1p_pimi->Sumw2();
h1_rot1_3pi_1p_pipl=new TH1F("3pi_1p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_3pi_1p_pipl->Sumw2();
h1_rot2_3pi_1p=new TH1F("3pi_1p_rot_kin_e","",n_bins, x_values);
  h1_rot2_3pi_1p->Sumw2();
h1_rot2_3pi_1p_pimi=new TH1F("3pi_1p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_3pi_1p_pimi->Sumw2();
h1_rot2_3pi_1p_pipl=new TH1F("3pi_1p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_3pi_1p_pipl->Sumw2();
h1_rot3_3pi_1p=new TH1F("3pi_1p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_3pi_1p->Sumw2();
h1_rot3_3pi_1p_pimi=new TH1F("3pi_1p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_3pi_1p_pimi->Sumw2();
h1_rot3_3pi_1p_pipl=new TH1F("3pi_1p_rot_kin_e_pi_pipl","",n_bins, x_values);
h1_rot3_3pi_1p_pipl->Sumw2();
h1_rot1_1pi_1p_1phot=new TH1F("1pi_1p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_1p_1phot->Sumw2();
h1_rot1_1pi_1p_1phot_pimi=new TH1F("1pi_1p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_1p_1phot_pimi->Sumw2();
h1_rot1_1pi_1p_1phot_pipl=new TH1F("1pi_1p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_1p_1phot_pipl->Sumw2();
h1_rot2_1pi_1p_1phot=new TH1F("1pi_1p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_1p_1phot->Sumw2();
h1_rot2_1pi_1p_1phot_pimi=new TH1F("1pi_1p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_1p_1phot_pimi->Sumw2();
h1_rot2_1pi_1p_1phot_pipl=new TH1F("1pi_1p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_1p_1phot_pipl->Sumw2();
h1_rot3_1pi_1p_1phot=new TH1F("1pi_1p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_1p_1phot->Sumw2();
h1_rot3_1pi_1p_1phot_pimi=new TH1F("1pi_1p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_1p_1phot_pimi->Sumw2();
h1_rot3_1pi_1p_1phot_pipl=new TH1F("1pi_1p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_1p_1phot_pipl->Sumw2();
h1_rot1_1pi_2p_1phot=new TH1F("1pi_2p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_2p_1phot->Sumw2();
h1_rot1_1pi_2p_1phot_pimi=new TH1F("1pi_2p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_2p_1phot_pimi->Sumw2();
h1_rot1_1pi_2p_1phot_pipl=new TH1F("1pi_2p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_2p_1phot_pipl->Sumw2();
h1_rot2_1pi_2p_1phot=new TH1F("1pi_2p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_2p_1phot->Sumw2();
h1_rot2_1pi_2p_1phot_pimi=new TH1F("1pi_2p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_2p_1phot_pimi->Sumw2();
h1_rot2_1pi_2p_1phot_pipl=new TH1F("1pi_2p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_2p_1phot_pipl->Sumw2();
h1_rot3_1pi_2p_1phot=new TH1F("1pi_2p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_2p_1phot->Sumw2();
h1_rot3_1pi_2p_1phot_pimi=new TH1F("1pi_2p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_2p_1phot_pimi->Sumw2();
h1_rot3_1pi_2p_1phot_pipl=new TH1F("1pi_2p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_2p_1phot_pipl->Sumw2();
h1_rot1_1pi_1p_2phot=new TH1F("1pi_1p_2phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_1p_2phot->Sumw2();
h1_rot1_1pi_1p_2phot_pimi=new TH1F("1pi_1p_2phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_1p_2phot_pimi->Sumw2();
h1_rot1_1pi_1p_2phot_pipl=new TH1F("1pi_1p_2phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_1p_2phot_pipl->Sumw2();
h1_rot2_1pi_1p_2phot=new TH1F("1pi_1p_2phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_1p_2phot->Sumw2();
h1_rot2_1pi_1p_2phot_pimi=new TH1F("1pi_1p_2phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_1p_2phot_pimi->Sumw2();
h1_rot2_1pi_1p_2phot_pipl=new TH1F("1pi_1p_2phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_1p_2phot_pipl->Sumw2();
h1_rot3_1pi_1p_2phot=new TH1F("1pi_1p_2phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_1p_2phot->Sumw2();
h1_rot3_1pi_1p_2phot_pimi=new TH1F("1pi_1p_2phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_1p_2phot_pimi->Sumw2();
h1_rot3_1pi_1p_2phot_pipl=new TH1F("1pi_1p_2phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_1p_2phot_pipl->Sumw2();
h1_rot1_2pi_1p_1phot=new TH1F("2pi_1p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_2pi_1p_1phot->Sumw2();
h1_rot1_2pi_1p_1phot_pimi=new TH1F("2pi_1p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_2pi_1p_1phot_pimi->Sumw2();
h1_rot1_2pi_1p_1phot_pipl=new TH1F("2pi_1p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_2pi_1p_1phot_pipl->Sumw2();
h1_rot2_2pi_1p_1phot=new TH1F("2pi_1p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_2pi_1p_1phot->Sumw2();
h1_rot2_2pi_1p_1phot_pimi=new TH1F("2pi_1p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_2pi_1p_1phot_pimi->Sumw2();
h1_rot2_2pi_1p_1phot_pipl=new TH1F("2pi_1p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_2pi_1p_1phot_pipl->Sumw2();
h1_rot3_2pi_1p_1phot=new TH1F("2pi_1p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_2pi_1p_1phot->Sumw2();
h1_rot3_2pi_1p_1phot_pimi=new TH1F("2pi_1p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_2pi_1p_1phot_pimi->Sumw2();
h1_rot3_2pi_1p_1phot_pipl=new TH1F("2pi_1p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_2pi_1p_1phot_pipl->Sumw2();
h1_rot1_1pi_3p_1phot=new TH1F("1pi_3p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_3p_1phot->Sumw2();
h1_rot1_1pi_3p_1phot_pimi=new TH1F("1pi_3p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_3p_1phot_pimi->Sumw2();
   h1_rot1_1pi_3p_1phot_pipl=new TH1F("1pi_3p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_3p_1phot_pipl->Sumw2();
   h1_rot2_1pi_3p_1phot=new TH1F("1pi_3p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_3p_1phot->Sumw2();
   h1_rot2_1pi_3p_1phot_pimi=new TH1F("1pi_3p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_3p_1phot_pimi->Sumw2();
   h1_rot2_1pi_3p_1phot_pipl=new TH1F("1pi_3p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_3p_1phot_pipl->Sumw2();
   h1_rot3_1pi_3p_1phot=new TH1F("1pi_3p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_3p_1phot->Sumw2();
   h1_rot3_1pi_3p_1phot_pimi=new TH1F("1pi_3p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_3p_1phot_pimi->Sumw2();
   h1_rot3_1pi_3p_1phot_pipl=new TH1F("1pi_3p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_3p_1phot_pipl->Sumw2();
/*   TH1F *h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01); */
/*   TH1F *h1_p_vert_corr=new TH1F("h1_p_vert_corr", "", 300, -10, 10); */
/*   TH1F *h1_pimi_vert_corr=new TH1F("h1_pimi_vert_corr", "", 300, -10, 10); */
/*   TH1F *h1_pipl_vert_corr=new TH1F("h1_pipl_vert_corr", "", 300, -10, 10); */
/*   TH1F *h1_p_vert_corr_cut=new TH1F("h1_p_vert_corr_cut", "", 300, -10, 10); */
/*   TH1F *h1_pimi_vert_corr_cut=new TH1F("h1_pimi_vert_corr_cut", "", 300, -10, 10); */
/*   TH1F *h1_pipl_vert_corr_cut=new TH1F("h1_pipl_vert_corr_cut", "", 300, -10, 10); */
/*   TH1F *h1_el_vertuncorr=new TH1F("h1_el_vertuncorr","",200,-10,10); */
/*   TH1F *h1_el_vertcorr=new TH1F("h1_el_vertcorr","",200,-10,10); */
/*   TH1F *h1_el_vertcorr_cut=new TH1F("h1_el_vertcorr_cut","",200,-10,10); */
/*   TH1F *h1_ec_beta_corr=new TH1F("h1_ec_beta_corr","",300,0,2); */
/*   TH1F *h1_ec_beta_corr_cut=new TH1F("h1_ec_beta_corr_cut","",300,0,2); */
/*   TH1F *h1_Q2 = new TH1F("h1_Q2","",200,0,6); */
/*   h1_Q2->Sumw2(); */
/*   TH1F *h1_omega = new TH1F("h1_omega","",200,0,5); */
/*   h1_omega->Sumw2(); */
/*   TH1F *h1_Wvar = new TH1F("h1_Wvar","",200,0,3); */
/*   h1_Wvar->Sumw2(); */
   h1_Q2_sub = new TH1F("h1_Q2_sub","",200,0,6); 
   h1_Q2_sub->Sumw2(); 
   h1_omega_sub = new TH1F("h1_omega_sub","",200,0,5);
   h1_omega_sub->Sumw2();
   h1_Wvar_sub = new TH1F("h1_Wvar_sub","",200,0,3);
   h1_Wvar_sub->Sumw2();
/*   TH1F *h1_feed_down_cal = new TH1F("h1_feed_down_cal","",200,0,1); */
/*   h1_feed_down_cal->Sumw2(); */
/*   TH1F *h1_feed_down_kin = new TH1F("h1_feed_down_kin","",200,0,1); */
/*   h1_feed_down_kin->Sumw2(); */
/*   TH2F *h2_e_ec_xy = new TH2F("h2_e_ec_xy","",100,-600,600,100,-600,600); */
/*   TH2F *h2_e_ec_xy_fidcut = new TH2F("h2_e_ec_xy_fidcut","",100,-600,600,100,-600,600); */
   h2_neut_costheta_phi=new TH2F("h2_neut_costheta_phi","",200,0,360,200,0,1.1);
   h2_neut_costheta_phi_cut=new TH2F("h2_neut_costheta_phi_cut","",200,0,360,200,0,1.1);
   h2_e_phi_theta=new TH2F("h2_e_phi_theta", "", 200,0,360,200,0,180);
   h2_e_phi_theta_cut=new TH2F("h2_e_phi_theta_cut", "", 200,0,360,200,0,180);
   h2_p_phi_theta=new TH2F("h2_p_phi_theta", "", 200,0,360,200,0,180);
   h2_p_phi_theta_cut=new TH2F("h2_p_phi_theta_cut", "", 200,0,360,200,0,180);
   h2_pimi_phi_theta=new TH2F("h2_pimi_phi_theta", "", 200,0,360,200,0,180);
   h2_pimi_phi_theta_cut=new TH2F("h2_pimi_phi_theta_cut", "", 200,0,360,200,0,180);
   h2_pipl_phi_theta=new TH2F("h2_pipl_phi_theta", "", 200,0,360,200,0,180);
   h2_pipl_phi_theta_cut=new TH2F("h2_pipl_phi_theta_cut", "", 200,0,360,200,0,180);
   h2_Np_Npi=new TH2F("h2_Np_Npi","",11,-0.5,4.5,11,-0.5, 4.5);
   h2_phot_pi_1p=new TH2F("h2_phot_pi_1p", "",11,-0.5,4.5,11,-0.5,4.5);
   h2_phot_pi_2p=new TH2F("h2_phot_pi_2p", "",11,-0.5,4.5,11,-0.5,4.5);
   h2_phot_pi_3p=new TH2F("h2_phot_pi_3p", "",11,-0.5,4.5,11,-0.5,4.5);
   prot_Deltat_p= new TH2F("prot_Deltat_p","",300,0.,2.7,300,-15,15);
   prot_Deltat_p_cut= new TH2F("prot_Deltat_p_cut","",300,0.,2.7,300,-15,15);
   pimi_delt_p= new TH2F("pimi_delt_p","",300,0.,2.5,300,-15,15);
   pimi_delt_p_cut= new TH2F("pimi_delt_p_cut","",300,0.,2.5,300,-15,15);
   pipl_delt_p= new TH2F("pipl_delt_p","",300,0.,2.5,300,-15,15);
   pipl_delt_p_cut= new TH2F("pipl_delt_p_cut","",300,0.,2.5,300,-15,15);
//   h2_el_E_p_ratio = new TH2F("h2_el_E_p_ratio","",200,0,4.5,200,0,0.5);
//   h2_el_E_p_ratio_cut = new TH2F("h2_el_E_p_ratio_cut","",200,0,4.5,200,0,0.5);
   h2_Wvar_Q2 = new TH2F("h2_Wvar_Q2","", 200,0,3,200,0,5);
  h2_Wvar_Q2->Sumw2();
   h2_Q2_omega = new TH2F("h2_Q2_omega","",200,0,3.5,200,0,5);
  h2_Q2_omega->Sumw2();
   h2_Wvar_Q2_sub = new TH2F("h2_Wvar_Q2_sub","", 200,0,3,200,0,5);
  h2_Wvar_Q2_sub->Sumw2();
   h2_Q2_omega_sub = new TH2F("h2_Q2_omega_sub","",200,0,3.5,200,0,5);
  h2_Q2_omega_sub->Sumw2();
   h2_kin_e_Wvar = new TH2F("h2_kin_e_Wvar","",200,0,5,200,0,5);
  h2_kin_e_Wvar->Sumw2();
   h2_kin_e_pi_Wvar = new TH2F("h2_kin_e_pi_Wvar","",200,0,5,200,0,5);
  h2_kin_e_pi_Wvar->Sumw2();
   h2_cal_Wvar = new TH2F("h2_cal_Wvar","",200,0,5,200,0,5);
  h2_cal_Wvar->Sumw2();
 h3_Npi_Np_Nphoton=new TH3F("h3_Npi_Np_Nphoton", "", 11, -0.5, 4.5, 11, -0.5, 4.5, 11, -0.5, 4.5);
