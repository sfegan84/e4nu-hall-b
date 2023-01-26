 //Definition and initialization of Histograms
h1_el_vertuncorr=new TH1F("h1_el_vertuncorr","",200,-10,10);
h1_el_vertcorr=new TH1F("h1_el_vertcorr","",200,-10,10);
h1_el_Mott_crosssec = new TH1F("h1_el_Mott_crosssec","",200,0.,0.01);
h1_el_Etot = new TH1F("h1_el_Etot","",600,0,4);
h1_el_Ein = new TH1F("h1_el_Ein","",800,0,1.5);
h1_el_Etot_cut = new TH1F("h1_el_Etot_cut","",600,0,4);
h1_el_Ein_cut = new TH1F("h1_el_Ein_cut","",800,0,1.5);
h1_el_cc_nphe = new TH1F("h1_el_cc_nphe","",300,0,30);
h1_el_cc_nphe_cut = new TH1F("h1_el_cc_nphe_cut","",300,0,30);
h1_el_cc_nphe_cut2 = new TH1F("h1_el_cc_nphe_cut2","",300,0,30);
h1_el_cc_chi2 = new TH1F("h1_el_cc_chi2","",200,0,0.7);
h1_Wvar = new TH1F("h1_Wvar","",400,0,3);
h1_Wvar_weight = new TH1F("h1_Wvar_weight","",400,0,3);
h1_Wepp = new TH1F("h1_Wepp","",400,0,3);
h1_Wepp_uncorr = new TH1F("h1_Wepp_uncorr","",400,0,3);
h1_xbjk = new TH1F("h1_xbjk","",400,0,3);
h1_xbjk_weight = new TH1F("h1_xbjk_weight","",400,0,3);
h1_Q2 = new TH1F("h1_Q2","",400,0,6);
h1_Q2_weight = new TH1F("h1_Q2_weight","",400,0,6);
h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
h1_Nprot=new TH1F("h1_Nprot","",10,0,5);
h1_Nprot_NonZeroProt=new TH1F("h1_Nprot_NonZeroProt","",8,1,5);
h1_Nphot=new TH1F("h1_Nphot","",10,0,5);
h1_Npiphot=new TH1F("h1_Npiphot","",10,0,5);
h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,0,5);
h1_photon_E=new TH1F("h1_photon_E","",200,0,4.5);
h1_photon_EC_E=new TH1F("h1_photon_EC_E","",200,0,4.5);
h1_phot_e_angle= new TH1F("h1_phot_e_angle","",300,0,180);
h1_time_ec=new TH1F("h1_time_ec","",200,-20,20);
h1_Npi=new TH1F("h1_Npi","",10,0,5);
h1_Npi_NonZeroProt=new TH1F("h1_Npi_NonZeroProt","",10,0,5);
h1_Npipl=new TH1F("h1_Npipl","",10,0,5);
h1_Npimi=new TH1F("h1_Npimi","",10,0,5);
h1_el_mom = new TH1F("h1_el_mom","",100,0.0,5.0);
h1_el_mom_corr = new TH1F("h1_el_mom_corr","",100,0.0,5.0);
h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
h1_el_prot_vertdiff_all= new TH1F("h1_el_prot_vertdiff_all","",300,-10,10);
//TH1F *h1_prot_vertdiff= new TH1F("h1_prot_vertdiff","",300,-10,10);
h1_el_prot_vertdiff= new TH1F("h1_el_prot_vertdiff","",300,-10,10);
h1_el_prot_vertdiff1= new TH1F("h1_el_prot_vertdiff1","",300,-10,10);
h1_el_prot_vertdiff2= new TH1F("h1_el_prot_vertdiff2","",300,-10,10);
h1_el_3prot_vertdiff1= new TH1F("h1_el_3prot_vertdiff1","",300,-10,10);
h1_el_3prot_vertdiff2= new TH1F("h1_el_3prot_vertdiff2","",300,-10,10);
h1_el_3prot_vertdiff3= new TH1F("h1_el_3prot_vertdiff3","",300,-10,10);
h1_el_4prot_vertdiff1= new TH1F("h1_el_4prot_vertdiff1","",300,-10,10);
h1_el_4prot_vertdiff2= new TH1F("h1_el_4prot_vertdiff2","",300,-10,10);
h1_el_4prot_vertdiff3= new TH1F("h1_el_4prot_vertdiff3","",300,-10,10);
h1_el_4prot_vertdiff4= new TH1F("h1_el_4prot_vertdiff4","",300,-10,10);
h1_pipl_prot_vertdiff= new TH1F("h1_pipl_prot_vertdiff","",300,-10,10);
h1_pimi_prot_vertdiff= new TH1F("h1_pimi_prot_vertdiff","",300,-10,10);
h1_prot_mom = new TH1F("h1_prot_mom","",300,0,3);
h1_prot_mom_ratio = new TH1F("h1_prot_mom_ratio","",50,0.97,1.2);
h1_prot_mom_pimi = new TH1F("h1_prot_mom_pimi","",300,0,3);
h1_prot_mom_ratio_pimi = new TH1F("h1_prot_mom_ratio_pimi","",50,0.97,1.2);
h1_prot_mom_pipl = new TH1F("h1_prot_mom_pipl","",300,0,3);
h1_prot_mom_ratio_pipl = new TH1F("h1_prot_mom_ratio_pipl","",50,0.97,1.2);
h1_pimi_mom = new TH1F("h1_pimi_mom","",300,0,3);
h1_pimi_mom_ratio = new TH1F("h1_pimi_mom_ratio","",50,0.97,1.2);
h1_pipl_mom = new TH1F("h1_pipl_mom","",300,0,3);
h1_pipl_mom_ratio = new TH1F("h1_pipl_mom_ratio","",50,0.97,1.2);
h1_pi0_mom = new TH1F("h1_pi0_mom","",300,0,3);
h1_pi0_gamma1_mom = new TH1F("h1_pi0_gamma1_mom","",300,0,3);
h1_pi0_gamma2_mom = new TH1F("h1_pi0_gamma2_mom","",300,0,3);
h1_neg_m=new TH1F("h1_neg_m","",300,0,2);
h1_pos_m=new TH1F("h1_pos_m","",300,0,2);
h1_kp_m=new TH1F("h1_kp_m","",300,0,2);
h1_km_m=new TH1F("h1_km_m","",300,0,2);


  //Binning for fractional feed down energy histograms
  int N_qe = -1;
  double *x_qe;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
    N_qe=109;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=64;i++)  x_qe[i]=-1+i*0.015;
    for (int i=0;i<=44;i++)  x_qe[i+65]=-0.04+(i+1)*0.01;
  }

  else if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
    N_qe=109;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=64;i++)  x_qe[i]=-1+i*0.015;
    for (int i=0;i<=44;i++)  x_qe[i+65]=-0.04+(i+1)*0.01;
  }

  else if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
    N_qe=82;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=29;i++)  x_qe[i]=-1+i*0.03;
    for (int i=0;i<=52;i++)  x_qe[i+30]=-0.13+(i+1)*0.01;
  }

  else{
    std::cout << "Beam Energy not determined. Aborting" << std::endl;
    x_qe = NULL;
    //continue;
  }

 //Definitions of further Histograms
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

h1_beta_ec = new TH1F("h1_beta_ec","",300,0,2);
h1_beta_ec_corr = new TH1F("h1_beta_ec_corr","",300,0,2);
h1_beta_ec_corr_cut = new TH1F("h1_beta_ec_corr_cut","",300,0,2);
h1_el_ec_sc_timediff = new TH1F("h1_el_ec_sc_timediff","",600,-30,30);
h1_el_ec_sc_timediff_corr = new TH1F("h1_el_ec_sc_timediff_corr","",600,-30,30);
h1_el_ec_sc_timediff_allSCpd = new TH1F("h1_el_ec_sc_timediff_allSCpd","",600,-30,30);
h1_el_ec_sc_timediff_corr_allSCpd = new TH1F("h1_el_ec_sc_timediff_corr_allSCpd","",600,-30,30);
h1_theta0=new TH1F("h1_theta0","",300,0,180);

h2_Ecal_Eqe=new TH2F("h2_Ecal_Eqe","",600,0,5,600,0,5);
h2_Ecal_Ekin=new TH2F("h2_Ecal_Ekin","",600,0,5,600,0,5);
h2_Ecal_Ekin_pipl=new TH2F("h2_Ecal_Ekin_pipl","",600,0,5,600,0,5);
h2_Ecal_Ekin_pimi=new TH2F("h2_Ecal_Ekin_pimi","",600,0,5,600,0,5);
h2_EqeEcalratio_Eqe=new TH2F("h2_EqeEcalratio_Eqe","",600,0,5,300,0,2);
h2_EqeEcaldiff_Eqe=new TH2F("h2_EqeEcaldiff_Eqe","",600,0,5,300,-3,3);
h2_N_prot_pi=new TH2F("h2_N_prot_pi","",10,0,5,10,0,5);
h2_N_prot_pi_phot=new TH2F("h2_N_prot_pi_phot","",10,0,5,10,0,5);
h2_N_prot_pi_phot_nonrad=new TH2F("h2_N_prot_pi_phot_nonrad","",10,0,5,10,0,5);
h2_el_E_p_ratio_cut = new TH2F("h2_el_E_p_ratio_cut","",200,0,4.5,200,0,0.5);
h2_el_E_p_ratio = new TH2F("h2_el_E_p_ratio","",200,0,4.5,200,0,0.5);
h2_el_E_p_ratio_withoutCC = new TH2F("h2_el_E_p_ratio_withoutCC","",200,0,4.5,200,0,0.5);
h2_el_E_p_ratio_withCC = new TH2F("h2_el_E_p_ratio_withCC","",200,0,4.5,200,0,0.5);

h2_el_Ein_Eout=new TH2D("h2_el_Ein_Eout","",400,0,2,400,0,1.5);
h2_el_Einout_Etot=new TH2D("h2_el_Einout_Etot","",400,0,4,400,0,4);

h2_el_ec_xy = new TH2F("h2_el_ec_xy","",100,-600,600,100,-600,600);
h2_el_ec_xy_fidcut = new TH2F("h2_el_ec_xy_fidcut","",100,-600,600,100,-600,600);
h2_el_phi_vert = new TH2F("h2_el_phi_vert","",500,-6,7,120,0,380);
h2_el_phi_vert_uncorr = new TH2F("h2_el_phi_vert_uncorr","",500,-6,7,120,0,380);
h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,0,180);
h2_neutral_costheta_phi_EC_all = new TH2F("h2_neutral_costheta_phi_EC_all","",200,0,360,200,0,1.1);
h2_neutral_theta_phi_EC_all = new TH2F("h2_neutral_theta_phi_EC_all","",200,0,360,200,0,180);
h2_neutral_theta_phi_EC_all_fidcut = new TH2F("h2_neutral_theta_phi_EC_all_fidcut","",200,0,360,200,0,180);
h2_pimi_theta_phi = new TH2F("h2_pimi_theta_phi","",200,0,360,200,0,180);
h2_pipl_theta_phi = new TH2F("h2_pipl_theta_phi","",200,0,360,200,0,180);
h2_kpl_theta_phi = new TH2F("h2_kpl_theta_phi","",200,0,360,200,0,180);
h2_pimi_theta_phi_beffid = new TH2F("h2_pimi_theta_phi_beffid","",200,0,360,200,0,180);
h2_pipl_theta_phi_beffid = new TH2F("h2_pipl_theta_phi_beffid","",200,0,360,200,0,180);
h2_kpl_theta_phi_beffid = new TH2F("h2_kpl_theta_phi_beffid","",200,0,360,200,0,180);
h2_prot_theta_phi = new TH2F("h2_prot_theta_phi","",200,0,360,200,0,180);
h2_prot_px_py_p = new TH2F("h2_prot_px_py_p","",100,-1,1,100,-1,1);
h2_prot_px_py_p_fidcut = new TH2F("h2_prot_px_py_p_fidcut","",100,-1,1,100,-1,1);
h2_pipl_theta_phi_p = new TH2F("h2_pipl_theta_phi_p","",200,0,360,200,0,180);
h2_pipl_theta_phi_fidcut = new TH2F("h2_pipl_theta_phi_fidcut","",200,0,360,200,0,180);
h2_kpl_theta_phi_p = new TH2F("h2_kpl_theta_phi_p","",200,0,360,200,0,180);
h2_kpl_theta_phi_fidcut = new TH2F("h2_kpl_theta_phi_fidcut","",200,0,360,200,0,180);
h2_pimi_theta_phi_p = new TH2F("h2_pimi_theta_phi_p","",200,0,360,200,0,180);
h2_pimi_theta_phi_fidcut = new TH2F("h2_pimi_theta_phi_fidcut","",200,0,360,200,0,180);
h2_el_mom_diff = new TH2F("h2_el_mom_diff","",500,0.,1.,500,-0.1,0.1);
h2_Q2_nu = new TH2F("h2_Q2_nu","",200,0,3.5,200,0,5);
h2_Q2_nu_weight = new TH2F("h2_Q2_nu_weight","",200,0,3.5,200,0,5);
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
h2_pipl_delt_p= new TH2F("h2_pipl_delt_p","",300,0.,2.5,300,-15,15);
h2_pimi_delt_p= new TH2F("h2_pimi_delt_p","",300,0.,2.5,300,-15,15);
h2_pos_delt_p= new TH2F("h2_pos_delt_p","",300,0.,2.5,300,-15,15);
h2_kp_delt_p= new TH2F("h2_kp_delt_p","",300,0.,2.5,300,-15,15);
h2_neg_delt_p= new TH2F("h2_neg_delt_p","",300,0.,2.5,300,-15,15);
h2_pipl_beta_p = new TH2F("h2_pipl_beta_p","",200,0,4,200,0,1.2);
h2_kp_beta_p = new TH2F("h2_kp_beta_p","",200,0,4,200,0,1.2);
h2_pimi_beta_p = new TH2F("h2_pimi_beta_p","",200,0,4,200,0,1.2);
h2_neg_beta_p = new TH2F("h2_neg_beta_p","",200,0,4,200,0,1.2);
h2_pos_beta_p = new TH2F("h2_pos_beta_p","",200,0,4,200,0,1.2);
h2_prot_beta_p = new TH2F("h2_prot_beta_p","",200,0,4,200,0,1.2);
h2_neg_E_p = new TH2F("h2_neg_E_p","",200,0,3,200,0,200);
h2_pos_E_p = new TH2F("h2_pos_E_p","",200,0,3,200,0,200);
h2_kp_E_p = new TH2F("h2_kp_E_p","",200,0,3,200,0,200);
h2_pimi_E_p = new TH2F("h2_pimi_E_p","",200,0,3,200,0,200);
h2_pipl_E_p = new TH2F("h2_pipl_E_p","",200,0,3,200,0,200);
h2_prot_E_p = new TH2F("h2_prot_E_p","",200,0,3,200,0,200);
h2_prot_Deltat_p= new TH2F("h2_prot_Deltat_p","",300,0.,2.7,300,-15,15);
h2_el_vertcorr_runN= new TH2F("h2_el_vertcorr_runN","",30,18370,18436,300,-7,7);
h2_phot_e_angle_vsphotE= new TH2F("h2_phot_e_angle_vsphotE","",300,0,3,300,0,180);
h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);
h2_Wepp_ephi= new TH2F("h2_Wepp_ephi","",720,0,360,200,0.85,1.05);
h2_Wepp_ephi_corr= new TH2F("h2_Wepp_ephi_corr","",720,0,360,200,0.85,1.05);
h2_Wepp_ephi_uncorrprot= new TH2F("h2_Wepp_ephi_uncorrprot","",720,0,360,200,0.85,1.05);
h2_Wepp_ephi_corr_uncorrprot= new TH2F("h2_Wepp_ephi_corr_uncorrprot","",720,0,360,200,0.85,1.05);

  //Binning for energy reconstruction histograms
  int n_bins = -1;
  double *x_values;

  if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
    n_bins=38;
    x_values=new double[n_bins+1];
    for (int i=0;i<=17;i++) x_values[i]=0.4+i*0.04;
    for (int i=0;i<=20;i++) x_values[i+18]=1.08+(i+1)*0.02;
  }

  else if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
    n_bins=54;
    x_values=new double[n_bins+1];
    for (int i=0;i<=23;i++) x_values[i]=i*0.09;
    for (int i=0;i<=30;i++) x_values[i+24]=2.07+(i+1)*0.03;
  }

  else if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
    n_bins=38;
    x_values=new double[n_bins+1];
    for (int i=0;i<=21;i++)  x_values[i]=i*0.2;
    for (int i=0;i<=16;i++)  x_values[i+22]=4.2+(i+1)*0.05;
  }

  else{
    std::cout << "Beam Energy not determined. Aborting" << std::endl;
    x_values = NULL;
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
h1_E_rec_1prot  = new TH1F("h1_E_rec_1prot","",n_bins,x_values);
h1_E_tot_1prot  = new TH1F("h1_E_tot_1prot","",n_bins,x_values);
h1_E_rec_cutpi1_piplpimi = new TH1F("h1_E_rec_cutpi1_piplpimi","",n_bins,x_values);
h1_E_tot_cutpi1_piplpimi = new TH1F("h1_E_tot_cutpi1_piplpimi","",n_bins,x_values);
h1_Etot = new TH1F("h1_Etot","",n_bins,x_values);
h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
h1_E_rec_undetfactor  = new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
h1_E_tot_undetfactor  = new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
h1_E_tot_undetfactor_pipl  = new TH1F("h1_E_tot_undetfactor_pipl","",n_bins,x_values);
h1_E_tot_undetfactor_pimi  = new TH1F("h1_E_tot_undetfactor_pimi","",n_bins,x_values);
h1_E_tot_1p2pi  = new TH1F("h1_E_tot_1p2pi","",n_bins,x_values);
h1_E_rec_1p2pi  = new TH1F("h1_E_rec_1p2pi","",n_bins,x_values);
h1_E_tot_1p3pi  = new TH1F("h1_E_tot_1p3pi","",n_bins,x_values);
h1_E_rec_1p3pi  = new TH1F("h1_E_rec_1p3pi","",n_bins,x_values);
h1_E_tot_2p2pi  = new TH1F("h1_E_tot_2p2pi","",n_bins,x_values);
h1_E_rec_2p2pi  = new TH1F("h1_E_rec_2p2pi","",n_bins,x_values);
h1_E_tot_3p1pi  = new TH1F("h1_E_tot_3p1pi","",n_bins,x_values);
h1_E_rec_3p1pi  = new TH1F("h1_E_rec_3p1pi","",n_bins,x_values);
h1_E_tot_1p2pi_1p0pi  = new TH1F("h1_E_tot_1p2pi_1p0pi","",n_bins,x_values);
h1_E_rec_1p2pi_1p0pi  = new TH1F("h1_E_rec_1p2pi_1p0pi","",n_bins,x_values);
h1_E_tot_2p1pi_2p0pi  = new TH1F("h1_E_tot_2p1pi_2p0pi","",n_bins,x_values);
h1_E_rec_2p1pi_2p0pi  = new TH1F("h1_E_rec_2p1pi_2p0pi","",n_bins,x_values);
h1_E_tot_2p1pi_1p1pi  = new TH1F("h1_E_tot_2p1pi_1p1pi","",n_bins,x_values);
h1_E_rec_2p1pi_1p1pi  = new TH1F("h1_E_rec_2p1pi_1p1pi","",n_bins,x_values);
h1_E_tot_2p1pi_1p0pi  = new TH1F("h1_E_tot_2p1pi_1p0pi","",n_bins,x_values);
h1_E_rec_2p1pi_1p0pi  = new TH1F("h1_E_rec_2p1pi_1p0pi","",n_bins,x_values);


  //Defintions of sector-specific histograms
  for(int h=0;h<nsect;h++) {
   h1_el_SCpdfidcut[h] = new TH1F(Form("h1_el_SCpdfidcut_%d",h+1),"",31,-0.5,30.5);
   h1_el_SCpd[h] = new TH1F(Form("h1_el_SCpd_%d",h+1),"",31,-0.5,30.5);
   h1_el_cc_deltat[h] = new TH1F(Form("h1_el_cc_deltat_%d",h+1),"",500,-100,100);
   h1_el_cc_deltat_cut[h] = new TH1F(Form("h1_el_cc_deltat_cut_%d",h+1),"",500,-100,100);
   h1_el_ec_sc_timediff_sect[h] = new TH1F(Form("h1_el_ec_sc_timediff_sect_%d",h+1),"",600,-30,30);
   h1_el_ec_sc_timediff_sect_corr[h] = new TH1F(Form("h1_el_ec_sc_timediff_sect_corr_%d",h+1),"",600,-30,30);
   h2_el_ec_sc_timediff_ecu[h] = new TH2F(Form("h2_el_ec_sc_timediff_ecu_%d",h+1),"",800,0,400,600,-30,30);
   h2_el_ec_sc_timediff_ecv[h] = new TH2F(Form("h2_el_ec_sc_timediff_ecv_%d",h+1),"",800,0,400,600,-30,30);
   h2_el_ec_sc_timediff_ecw[h] = new TH2F(Form("h2_el_ec_sc_timediff_ecw_%d",h+1),"",800,0,400,600,-30,30);
   h2_el_ec_sc_timediff_SCpd[h] = new TH2F(Form("h2_el_ec_sc_timediff_SCpd_%d",h+1),"",75,0,25,600,-30,30);
   h1_beta_ec_corr_sect[h] = new TH1F(Form("h1_beta_ec_corr_sect_%d",h+1),"",300,0,1.3);
   h1_e_mom_corrfuct[h] = new TH1F(Form("h1_e_mom_corrfuct_%d",h+1),"",400,0.7,1.3);
   h2_el_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_el_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_el_theta_phi_p_fidcut[h] = new TH2F(Form("h2_el_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_el_theta_phi_p_beffidcut2[h] = new TH2F(Form("h2_el_theta_phi_p_beffidcut2_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_el_theta_phi_p_fidcut2[h] = new TH2F(Form("h2_el_theta_phi_p_fidcut2_%d",h+1),"",200,h*60,(h+1)*60,200,0,70);
   h2_prot_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_prot_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_prot_theta_phi_p_fidcut[h] = new TH2F(Form("h2_prot_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pipl_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_pipl_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pipl_theta_phi_p_fidcut[h] = new TH2F(Form("h2_pipl_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_kpl_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_kpl_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_kpl_theta_phi_p_fidcut[h] = new TH2F(Form("h2_kpl_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pimi_theta_phi_p_beffidcut[h] = new TH2F(Form("h2_pimi_theta_phi_p_beffidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_pimi_theta_phi_p_fidcut[h] = new TH2F(Form("h2_pimi_theta_phi_p_fidcut_%d",h+1),"",200,h*60,(h+1)*60,200,0,180);
   h2_el_theta_p[h] = new TH2F(Form("h2_el_theta_p_%d",h+1),"",400,0,4,360,0,180);
   h2_el_theta_p_cut[h] = new TH2F(Form("h2_el_theta_p_cut_%d",h+1),"",400,0,4,360,0,180);
   h2_prot_theta_p[h] = new TH2F(Form("h2_prot_theta_p_%d",h+1),"",400,0,4,360,0,180);
   h2_prot_theta_p_cut[h] = new TH2F(Form("h2_prot_theta_p_cut_%d",h+1),"",400,0,4,360,0,180);
   h2_pimi_theta_p[h] = new TH2F(Form("h2_pimi_theta_p_%d",h+1),"",800,0,4,360,0,180);
   h2_pimi_theta_p_cut[h] = new TH2F(Form("h2_pimi_theta_p_cut_%d",h+1),"",400,0,4,360,0,180);
   h2_pipl_theta_p[h] = new TH2F(Form("h2_pipl_theta_p_%d",h+1),"",800,0,4,360,0,180);
   h2_pipl_theta_p_cut[h] = new TH2F(Form("h2_pipl_theta_p_cut_%d",h+1),"",400,0,4,360,0,180);
   h2_kpl_theta_p[h] = new TH2F(Form("h2_kpl_theta_p_%d",h+1),"",800,0,4,360,0,180);
   h2_kpl_theta_p_cut[h] = new TH2F(Form("h2_kpl_theta_p_cut_%d",h+1),"",400,0,4,360,0,180);
  }


  for(int h=0;h<n_slice;h++){
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
   h1_Erec_Npi0_new[h] = new TH1F(Form("h1_Erec_Npi0_new_%d",h+1),"",n_bins,x_values);
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


h1_E_tot_undetfactor09  = new TH1F("h1_E_tot_undetfactor09","",1,0,6);
h1_E_tot_cut2_09 = new TH1F("h1_E_tot_cut2_09","",1,0,6);
h1_E_tot_p_bkgd09  = new TH1F("h1_E_tot_p_bkgd09","",1,0,6);
h1_Etot_p321_bkgd09  = new TH1F("h1_Etot_p321_bkgd09","",1,0,6);
h1_Etot_p31_bkgd09  = new TH1F("h1_Etot_p31_bkgd09","",1,0,6);
h1_Etot_p4321_bkgd09  = new TH1F("h1_Etot_p4321_bkgd09","",1,0,6);
h1_Etot_p431_bkgd09  = new TH1F("h1_Etot_p431_bkgd09","",1,0,6);
h1_Etot_p421_bkgd09  = new TH1F("h1_Etot_p421_bkgd09","",1,0,6);
h1_Etot_p41_bkgd09  = new TH1F("h1_Etot_p41_bkgd09","",1,0,6);
h1_Etot_bkgd09_2p1pi_2p0pi  = new TH1F("h1_Etot_bkgd09_2p1pi_2p0pi","",1,0,6);
h1_Etot_bkgd09_2p1pi_1p1pi  = new TH1F("h1_Etot_bkgd09_2p1pi_1p1pi","",1,0,6);
h1_Etot_bkgd09_2p1pi_1p0pi  = new TH1F("h1_Etot_bkgd09_2p1pi_1p0pi","",1,0,6);
h1_Etot_bkgd09_1p2pi_1p0pi  = new TH1F("h1_Etot_bkgd09_1p2pi_1p0pi","",1,0,6);
h1_Etot_bkgd09_1p2pi_1p1pi  = new TH1F("h1_Etot_bkgd09_1p2pi_1p1pi","",1,0,6);
h1_Etot_bkgd09_1p3pi  = new TH1F("h1_Etot_bkgd09_1p3pi","",1,0,6);
h1_Etot_bkgd09_2p2pi  = new TH1F("h1_Etot_bkgd09_2p2pi","",1,0,6);
h1_Etot_bkgd09_3p1pi  = new TH1F("h1_Etot_bkgd09_3p1pi","",1,0,6);


  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){

      h1_Etot_Npi0_Ecalcut[i][j]=  new TH1F(Form("h1_Etot_Npi0_Ecalcut_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j]=  new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_Ecalcut_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut41[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut41_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut421[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut421_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut431[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut431_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut4321[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut4321_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut31[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut31_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut321[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut321_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut3_2p1pito2p0pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_1p3pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_2p2pi_%d_%d",i+1,j+1),"",1,0,6);
      h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[i][j]=  new TH1F(Form("h1_Etot_p_bkgd_slice_Ecalcut_3p1pi_%d_%d",i+1,j+1),"",1,0,6);
    }
  }



//1p 1pi histos
h1_en_recon1=new TH1F("calorimetric","", n_bins, x_values);
  h1_en_recon1->Sumw2();
h1_en_recon1_pimi=new TH1F("cal_pimi", "", n_bins, x_values);
  h1_en_recon1_pimi->Sumw2();
h1_en_recon1_pipl=new TH1F("cal_pipl", "",n_bins, x_values);
  h1_en_recon1_pipl->Sumw2();
h1_en_recon1_pi0=new TH1F("cal_pi0", "",n_bins, x_values);
  h1_en_recon1_pi0->Sumw2();
h1_en_recon1_eta=new TH1F("cal_eta", "",n_bins, x_values);
  h1_en_recon1_eta->Sumw2();
h1_en_recon1_klambda=new TH1F("cal_klambda", "",n_bins, x_values);
  h1_en_recon1_klambda->Sumw2();
h1_en_recon1_kpkm=new TH1F("cal_kpkm", "",n_bins, x_values);
  h1_en_recon1_kpkm->Sumw2();
h1_en_recon1_pippim=new TH1F("cal_pippim", "",n_bins, x_values);
  h1_en_recon1_pippim->Sumw2();
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

h2_cal_kin=new TH2F("cal_vs_kin","", n_bins, x_values, n_bins, x_values);
h2_cal_kin->Sumw2();
h2_cal_kin_pimi=new TH2F("cal_vs_kin_pimi", "", n_bins, x_values, n_bins, x_values);
h2_cal_kin_pimi->Sumw2();
h2_cal_kin_pipl=new TH2F("cal_vs_kin_pipl", "",n_bins, x_values, n_bins, x_values);
h2_cal_kin_pipl->Sumw2();


h2_rot_1pi_2p=new TH2F("1pi_2p_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_2p_pimi=new TH2F("1pi_2p_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_2p_pipl=new TH2F("1pi_2p_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_2pi_1p=new TH2F("2pi_1p_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_2pi_1p_pimi=new TH2F("2pi_1p_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_2pi_1p_pipl=new TH2F("2pi_1p_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_2pi_2p=new TH2F("2pi_2p_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_2pi_2p_pimi=new TH2F("2pi_2p_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_2pi_2p_pipl=new TH2F("2pi_2p_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_1pi_3p=new TH2F("1pi_3p_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_3p_pimi=new TH2F("1pi_3p_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_3p_pipl=new TH2F("1pi_3p_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_3pi_1p=new TH2F("3pi_1p_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_3pi_1p_pimi=new TH2F("3pi_1p_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_3pi_1p_pipl=new TH2F("3pi_1p_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_1pi_1p_1phot=new TH2F("1pi_1p_1phot_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_1p_1phot_pimi=new TH2F("1pi_1p_1phot_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_1p_1phot_pipl=new TH2F("1pi_1p_1phot_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_1pi_2p_1phot=new TH2F("1pi_2p_1phot_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_2p_1phot_pimi=new TH2F("1pi_2p_1phot_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_2p_1phot_pipl=new TH2F("1pi_2p_1phot_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_1pi_1p_2phot=new TH2F("1pi_1p_2phot_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_1p_2phot_pimi=new TH2F("1pi_1p_2phot_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_1p_2phot_pipl=new TH2F("1pi_1p_2phot_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_2pi_1p_1phot=new TH2F("2pi_1p_1phot_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_2pi_1p_1phot_pimi=new TH2F("2pi_1p_1phot_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_2pi_1p_1phot_pipl=new TH2F("2pi_1p_1phot_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_1pi_3p_1phot=new TH2F("1pi_3p_1phot_rot_kin_cal","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_3p_1phot_pimi=new TH2F("1pi_3p_1phot_rot_kin_cal_pimi","",n_bins, x_values,n_bins, x_values);
h2_rot_1pi_3p_1phot_pipl=new TH2F("1pi_3p_1phot_rot_kin_cal_pipl","",n_bins, x_values,n_bins, x_values);

h2_rot_1pi_2p->Sumw2();
h2_rot_2pi_1p->Sumw2();
h2_rot_2pi_2p->Sumw2();
h2_rot_1pi_3p->Sumw2();
h2_rot_3pi_1p->Sumw2();
h2_rot_1pi_1p_1phot->Sumw2();
h2_rot_1pi_2p_1phot->Sumw2();
h2_rot_1pi_1p_2phot->Sumw2();
h2_rot_2pi_1p_1phot->Sumw2();
h2_rot_1pi_3p_1phot->Sumw2();

h2_rot_1pi_2p_pimi->Sumw2();
h2_rot_2pi_1p_pimi->Sumw2();
h2_rot_2pi_2p_pimi->Sumw2();
h2_rot_1pi_3p_pimi->Sumw2();
h2_rot_3pi_1p_pimi->Sumw2();
h2_rot_1pi_1p_1phot_pimi->Sumw2();
h2_rot_1pi_2p_1phot_pimi->Sumw2();
h2_rot_1pi_1p_2phot_pimi->Sumw2();
h2_rot_2pi_1p_1phot_pimi->Sumw2();
h2_rot_1pi_3p_1phot_pimi->Sumw2();

h2_rot_1pi_2p_pipl->Sumw2();
h2_rot_2pi_1p_pipl->Sumw2();
h2_rot_2pi_2p_pipl->Sumw2();
h2_rot_1pi_3p_pipl->Sumw2();
h2_rot_3pi_1p_pipl->Sumw2();
h2_rot_1pi_1p_1phot_pipl->Sumw2();
h2_rot_1pi_2p_1phot_pipl->Sumw2();
h2_rot_1pi_1p_2phot_pipl->Sumw2();
h2_rot_2pi_1p_1phot_pipl->Sumw2();
h2_rot_1pi_3p_1phot_pipl->Sumw2();

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
   h2_W_Ecal_pi0 = new TH2F("h2_W_Ecal_pi0","",200,0,5,200,0,5);
  h2_W_Ecal_pi0->Sumw2();
   h2_W_Ecal_eta = new TH2F("h2_W_Ecal_eta","",200,0,5,200,0,5);
  h2_W_Ecal_eta->Sumw2();
   h2_W_Ecal_pi0_2 = new TH2F("h2_W_Ecal_pi0_2","",200,0,5,200,0,5);
  h2_W_Ecal_pi0_2->Sumw2();
   h2_W_Ecal_eta2 = new TH2F("h2_W_Ecal_eta2","",200,0,5,200,0,5);
  h2_W_Ecal_eta2->Sumw2();
   h2_Ecal_pi0_p_perp = new TH2F("h2_Ecal_pi0_p_perp","",200,0,5,200,0,5);
  h2_Ecal_pi0_p_perp->Sumw2();
   h2_Ecal_pip_p_perp = new TH2F("h2_Ecal_pip_p_perp","",200,0,5,200,0,5);
  h2_Ecal_pip_p_perp->Sumw2();
   h2_Ecal_pim_p_perp = new TH2F("h2_Ecal_pim_p_perp","",200,0,5,200,0,5);
  h2_Ecal_pim_p_perp->Sumw2();
   h2_Ecal_eta_p_perp = new TH2F("h2_Ecal_eta_p_perp","",200,0,5,200,0,5);
  h2_Ecal_eta_p_perp->Sumw2();

 h3_Npi_Np_Nphoton=new TH3F("h3_Npi_Np_Nphoton", "", 11, -0.5, 4.5, 11, -0.5, 4.5, 11, -0.5, 4.5);





  //initialization of correct error arrays for the histograms
  h1_E_tot_cut2_09->Sumw2();
  h1_E_tot_undetfactor09->Sumw2();
  h1_E_tot_p_bkgd09->Sumw2();
  h1_Etot_p321_bkgd09->Sumw2();
  h1_Etot_p31_bkgd09->Sumw2();
  h1_Etot_p4321_bkgd09->Sumw2();
  h1_Etot_p431_bkgd09->Sumw2();
  h1_Etot_p421_bkgd09->Sumw2();
  h1_Etot_p41_bkgd09->Sumw2();
  h1_Etot_bkgd09_2p1pi_2p0pi->Sumw2();
  h1_Etot_bkgd09_2p1pi_1p1pi->Sumw2();
  h1_Etot_bkgd09_2p1pi_1p0pi->Sumw2();
  h1_Etot_bkgd09_1p2pi_1p0pi->Sumw2();
  h1_Etot_bkgd09_1p2pi_1p1pi->Sumw2();
  h1_Etot_bkgd09_1p3pi->Sumw2();
  h1_Etot_bkgd09_2p2pi->Sumw2();
  h1_Etot_bkgd09_3p1pi->Sumw2();

  for (int i=0;i<N_pperp;i++){
    for(int j=0;j<N_Ecal;j++){
      h1_Etot_Npi0_Ecalcut[i][j]->Sumw2();
      h1_Etot_bkgd_pipl_pimi_fact_Ecalcut[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut41[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut421[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut431[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut4321[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut31[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut321[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito2p0pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p1pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p1pito1p0pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p0pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_1p2pito1p1pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_1p3pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_2p2pi[i][j]->Sumw2();
      h1_Etot_p_bkgd_slice_Ecalcut_3p1pi[i][j]->Sumw2();
    }
  }

  for(int h=0;h<n_slice;h++){
   h1_Erec_p_bkgd_slice[h]->Sumw2();
   h1_Etot_p_bkgd_slice[h]->Sumw2();
   h1_Erec_3pto1p_slice[h]->Sumw2();
   h1_Etot_3pto1p_slice[h]->Sumw2();
   h1_Erec_3pto2p_slice[h]->Sumw2();
   h1_Etot_3pto2p_slice[h]->Sumw2();
   h1_Erec_3p1pi_slice[h]->Sumw2();
   h1_Etot_3p1pi_slice[h]->Sumw2();
   h1_Etot_43pto1p_slice[h]->Sumw2();
   h1_Erec_43pto1p_slice[h]->Sumw2();
   h1_Erec_4pto3p_slice[h]->Sumw2();
   h1_Etot_4pto3p_slice[h]->Sumw2();
   h1_Erec_4pto2p_slice[h]->Sumw2();
   h1_Etot_4pto2p_slice[h]->Sumw2();
   h1_Erec_4pto1p_slice[h]->Sumw2();
   h1_Etot_4pto1p_slice[h]->Sumw2();
   h1_Etot_Npi0[h]->Sumw2();
   h1_Erec_Npi0_new[h]->Sumw2();
   h1_Etot_bkgd_pipl_pimi_fact[h]->Sumw2();
   h1_Etot_bkgd_pipl_pimi_fact_pipl[h]->Sumw2();
   h1_Etot_bkgd_pipl_pimi_fact_pimi[h]->Sumw2();
   h1_Erec_bkgd_pipl_pimi_new_fact[h]->Sumw2();
   h1_Etot_Npi1[h]->Sumw2();
   h1_Erec_Npi1[h]->Sumw2();
   h1_Etot_bkgd_1p2pi[h]->Sumw2();
   h1_Erec_bkgd_1p2pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p1pi_to1p1pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p1pi_to2p0pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p1pi_to1p1pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p1pi_to2p0pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p2pi[h]->Sumw2();
   h1_Erec_p_bkgd_slice_2p1pi_to1p0pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p1pi_to1p0pi[h]->Sumw2();
   h1_Etot_p_bkgd_slice_2p2pi[h]->Sumw2();
   h1_Etot_bkgd_1p2pi_1p0pi[h]->Sumw2();
   h1_Erec_bkgd_1p2pi_1p0pi[h]->Sumw2();
   h1_Erec_bkgd_1p3pi[h]->Sumw2();
   h1_Etot_bkgd_1p3pi[h]->Sumw2();
  }

  h1_E_rec_1pi_weight_frac_feed->Sumw2();
  h1_E_rec_2pi_weight_frac_feed->Sumw2();
  h1_E_rec_3pi_weight_frac_feed->Sumw2();
  h1_E_rec_4pi_weight_frac_feed->Sumw2();
  h1_E_rec_0pi_frac_feed->Sumw2();
  h1_E_tot_cut2_fracfeed ->Sumw2();
  h1_E_rec_cut2_new_fracfeed ->Sumw2();
  h1_E_tot_p_bkgd_fracfeed ->Sumw2();
  h1_E_rec_p_bkgd_fracfeed ->Sumw2();
  h1_E_tot_2p1pi_2p0pi_fracfeed ->Sumw2();
  h1_E_rec_2p1pi_2p0pi_fracfeed ->Sumw2();
  h1_E_tot_2p1pi_1p1pi_fracfeed ->Sumw2();
  h1_E_rec_2p1pi_1p1pi_fracfeed ->Sumw2();
  h1_E_tot_2p1pi_1p0pi_fracfeed ->Sumw2();
  h1_E_rec_2p1pi_1p0pi_fracfeed ->Sumw2();
  h1_E_tot_1p3pi_fracfeed ->Sumw2();
  h1_E_rec_1p3pi_fracfeed ->Sumw2();
  h1_E_tot_2p2pi_fracfeed ->Sumw2();
  h1_E_rec_2p2pi_fracfeed ->Sumw2();
  h1_E_tot_3p1pi_fracfeed ->Sumw2();
  h1_E_rec_3p1pi_fracfeed ->Sumw2();
  h1_E_tot_3pto2p_fracfeed->Sumw2();
  h1_E_rec_3pto2p_fracfeed->Sumw2();
  h1_E_tot_3pto1p_fracfeed->Sumw2();
  h1_E_rec_3pto1p_fracfeed->Sumw2();
  h1_E_tot_4pto3p_fracfeed->Sumw2();
  h1_E_rec_4pto3p_fracfeed->Sumw2();
  h1_E_tot_43pto1p_fracfeed->Sumw2();
  h1_E_rec_43pto1p_fracfeed->Sumw2();
  h1_E_tot_4pto2p_fracfeed->Sumw2();
  h1_E_rec_4pto2p_fracfeed->Sumw2();
  h1_E_tot_4pto1p_fracfeed->Sumw2();
  h1_E_rec_4pto1p_fracfeed->Sumw2();
  h1_E_tot_1p2pi_fracfeed->Sumw2();
  h1_E_rec_1p2pi_fracfeed->Sumw2();
  h1_E_tot_1p2pi_1p0pi_fracfeed->Sumw2();
  h1_E_rec_1p2pi_1p0pi_fracfeed->Sumw2();
  h1_E_rec_undetfactor_fracfeed->Sumw2();
  h1_E_tot_undetfactor_fracfeed->Sumw2();
  h1_E_rec_2p_det->Sumw2();
  h1_E_tot_2p_det->Sumw2();
  h1_E_tot_p_bkgd->Sumw2();
  h1_E_rec_p_bkgd->Sumw2();
  h1_E_tot_3pto1p->Sumw2();
  h1_E_rec_3pto1p->Sumw2();
  h1_E_tot_43pto1p->Sumw2();
  h1_E_rec_43pto1p->Sumw2();
  h1_E_tot_3pto2p->Sumw2();
  h1_E_rec_3pto2p->Sumw2();
  h1_E_tot_4pto1p->Sumw2();
  h1_E_rec_4pto1p->Sumw2();
  h1_E_tot_4pto3p->Sumw2();
  h1_E_rec_4pto3p->Sumw2();
  h1_E_tot_4pto2p->Sumw2();
  h1_E_rec_4pto2p->Sumw2();
  h1_E_rec->Sumw2();
  h1_E_rec_0pi->Sumw2();
  h1_E_rec_1pi->Sumw2();
  h1_E_rec_4pi_weight->Sumw2();
  h1_E_rec_40pi->Sumw2();
  h1_E_rec_410pi->Sumw2();
  h1_E_rec_420pi->Sumw2();
  h1_E_rec_4210pi->Sumw2();
  h1_E_rec_430pi->Sumw2();
  h1_E_rec_4310pi->Sumw2();
  h1_E_rec_4320pi->Sumw2();
  h1_E_rec_43210pi->Sumw2();
  h1_E_rec_1pi_weight->Sumw2();
  h1_E_rec_2pi_weight->Sumw2();
  h1_E_rec_3pi_weight->Sumw2();
  h1_E_rec_20pi->Sumw2();
  h1_E_rec_21pi->Sumw2();
  h1_E_rec_30pi->Sumw2();
  h1_E_rec_310pi->Sumw2();
  h1_E_rec_320pi->Sumw2();
  h1_E_rec_3210pi->Sumw2();
  h1_E_rec_1prot->Sumw2();
  h1_E_tot_1prot->Sumw2();
  h1_E_rec_cutpi1_piplpimi->Sumw2();
  h1_E_tot_cutpi1_piplpimi->Sumw2();
  h1_Etot->Sumw2();
  h1_E_rec_cut2_new->Sumw2();
  h1_E_tot_cut2->Sumw2();
  h1_E_rec_cut005_newcut3->Sumw2();
  h1_E_rec_undetfactor->Sumw2();
  h1_E_tot_undetfactor->Sumw2();
  h1_E_tot_undetfactor_pipl->Sumw2();
  h1_E_tot_undetfactor_pimi->Sumw2();
  h1_E_tot_1p2pi->Sumw2();
  h1_E_rec_1p2pi->Sumw2();
  h1_E_tot_1p3pi->Sumw2();
  h1_E_rec_1p3pi->Sumw2();
  h1_E_tot_2p2pi->Sumw2();
  h1_E_rec_2p2pi->Sumw2();
  h1_E_tot_3p1pi->Sumw2();
  h1_E_rec_3p1pi->Sumw2();
  h1_E_tot_1p2pi_1p0pi->Sumw2();
  h1_E_rec_1p2pi_1p0pi->Sumw2();
  h1_E_tot_2p1pi_2p0pi->Sumw2();
  h1_E_rec_2p1pi_2p0pi->Sumw2();
  h1_E_tot_2p1pi_1p1pi->Sumw2();
  h1_E_rec_2p1pi_1p1pi->Sumw2();
  h1_E_tot_2p1pi_1p0pi->Sumw2();
  h1_E_rec_2p1pi_1p0pi->Sumw2();
