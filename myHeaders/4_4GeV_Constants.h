prot_mom_lim=2.7;

if(ftarget=="3He"){//the EC threshold was different for He3
  min_good_mom=1.3;
 }
 else{
   min_good_mom=1.1;
 }

max_mom=3.7;
pipl_maxmom=1.9;
pimi_maxmom=1.6;


for (int i=0;i<N_Ecal;i++){
  Ecal_lowlim[i]=1.5+i*0.5;
  Ecal_uplim[i]=2.+i*0.5;
 }

Ecal_lowlim[5]=0.;
Ecal_uplim[5]=4.;

for (int i=0;i<N_Ecal;i++){
  std::cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<std::endl;;
 }

vert_min["3He"]=-3.27;
vert_min["4He"]=-2.51;
vert_min["C12"]=4.7;
vert_min["56Fe"]=4.6;

vert_max["3He"]=0.07;
vert_max["4He"]=1.71;
vert_max["C12"]=5.3;
vert_max["56Fe"]=5.4;


EC_photon_beta["3He"]=0.92;
EC_photon_beta["4He"]=0.91;
EC_photon_beta["C12"]=0.92;
EC_photon_beta["56Fe"]=0.91;


LEC_photon_beta["3He"]=0.97;
LEC_photon_beta["4He"]=0.97;
LEC_photon_beta["C12"]=0.95;
LEC_photon_beta["56Fe"]=0.96;

EC_time_offset[std::make_pair("3He",1)]=-0.15;  EC_time_offset[std::make_pair("3He",2)]=-0.26; EC_time_offset[std::make_pair("3He",3)]=-0.41;
EC_time_offset[std::make_pair("3He",4)]=-0.29;  EC_time_offset[std::make_pair("3He",5)]=-0.25; EC_time_offset[std::make_pair("3He",6)]=-0.23;

EC_time_offset[std::make_pair("4He",1)]=-0.01;  EC_time_offset[std::make_pair("4He",2)]=-0.11; EC_time_offset[std::make_pair("4He",3)]=-0.23;
EC_time_offset[std::make_pair("4He",4)]=-0.26;  EC_time_offset[std::make_pair("4He",5)]=-0.21; EC_time_offset[std::make_pair("4He",6)]=-0.09;

EC_time_offset[std::make_pair("C12",1)]=-0.01;  EC_time_offset[std::make_pair("C12",2)]=-0.11; EC_time_offset[std::make_pair("C12",3)]=-0.23;
EC_time_offset[std::make_pair("C12",4)]=-0.27;  EC_time_offset[std::make_pair("C12",5)]=-0.21; EC_time_offset[std::make_pair("C12",6)]=-0.08;

EC_time_offset[std::make_pair("56Fe",1)]=-0.49;  EC_time_offset[std::make_pair("56Fe",2)]=-0.14; EC_time_offset[std::make_pair("56Fe",3)]=-0.32;
EC_time_offset[std::make_pair("56Fe",4)]=-0.25;  EC_time_offset[std::make_pair("56Fe",5)]=-0.17; EC_time_offset[std::make_pair("56Fe",6)]=-0.35;

elmom_corr_fact[0]=1.001;
elmom_corr_fact[1]=0.991;
elmom_corr_fact[2]=1.005;
elmom_corr_fact[3]=1.004;
elmom_corr_fact[4]=1.006;
elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
