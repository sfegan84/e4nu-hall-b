prot_mom_lim=0.95;
min_good_mom=0.4;
max_mom=1.1;
pipl_maxmom=0.65;
pimi_maxmom=0.6;


for (int i=0;i<N_Ecal;i++){
  Ecal_lowlim[i]=0.45+i*0.18;
  Ecal_uplim[i]=0.63+i*0.18;
 }

Ecal_lowlim[5]=0.0;
Ecal_uplim[5]=1.35;

for (int i=0;i<N_Ecal;i++){
  std::cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<std::endl;
 }

vert_min["3He"]=-3.05;
vert_min["C12"]=4.95;
vert_min["CH2"]=4.85;

vert_max["3He"]=-0.18;
vert_max["C12"]=5.76;
vert_max["CH2"]=5.62;


EC_photon_beta["3He"]=0.89;
EC_photon_beta["C12"]=0.89;
EC_photon_beta["CH2"]=0.91;

LEC_photon_beta["3He"]=0.97;
LEC_photon_beta["C12"]=0.97;
LEC_photon_beta["CH2"]=0.97;

EC_time_offset[std::make_pair("3He",1)]=-0.73;  EC_time_offset[std::make_pair("3He",2)]=-0.81; EC_time_offset[std::make_pair("3He",3)]=-0.91;
EC_time_offset[std::make_pair("3He",4)]=-0.94;  EC_time_offset[std::make_pair("3He",5)]=-0.92; EC_time_offset[std::make_pair("3He",6)]=-0.81;

EC_time_offset[std::make_pair("C12",1)]=-0.71;  EC_time_offset[std::make_pair("C12",2)]=-0.77; EC_time_offset[std::make_pair("C12",3)]=-0.87;
EC_time_offset[std::make_pair("C12",4)]=-0.91;  EC_time_offset[std::make_pair("C12",5)]=-0.89; EC_time_offset[std::make_pair("C12",6)]=-0.79;

EC_time_offset[std::make_pair("CH2",1)]=-0.70;  EC_time_offset[std::make_pair("CH2",2)]=-0.80; EC_time_offset[std::make_pair("CH2",3)]=-0.91;
EC_time_offset[std::make_pair("CH2",4)]=-0.92;  EC_time_offset[std::make_pair("CH2",5)]=-0.91; EC_time_offset[std::make_pair("CH2",6)]=-0.80;

elmom_corr_fact[0]=1.007;
elmom_corr_fact[1]=0.988;
elmom_corr_fact[2]=1.008;
elmom_corr_fact[3]=1.011;
elmom_corr_fact[4]=1.014;
elmom_corr_fact[5]=1.013;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
