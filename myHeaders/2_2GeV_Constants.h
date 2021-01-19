prot_mom_lim=2.15;
min_good_mom=0.55;
max_mom=2.1;
pipl_maxmom=1.4;
pimi_maxmom=1.3;


for (int i=0;i<N_Ecal;i++){
  Ecal_lowlim[i]=0.75+i*0.25;
  Ecal_uplim[i]=1.+i*0.25;
 }

Ecal_lowlim[5]=0.;
Ecal_uplim[5]=2.;

for (int i=0;i<N_Ecal;i++){
  std::cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<std::endl;
 }

vert_min["3He"]=-3.29;
vert_min["4He"]=-2.53;
vert_min["C12"]=4.8;
vert_min["56Fe"]=4.6;

vert_max["3He"]=-0.23;
vert_max["4He"]=1.73;
vert_max["C12"]=5.5;
vert_max["56Fe"]=5.3;


EC_photon_beta["3He"]=0.93;
EC_photon_beta["4He"]=0.92;
EC_photon_beta["C12"]=0.92;
EC_photon_beta["56Fe"]=0.90;

LEC_photon_beta["3He"]=0.96;
LEC_photon_beta["4He"]=0.94;
LEC_photon_beta["C12"]=0.94;
LEC_photon_beta["56Fe"]=0.95;

EC_time_offset[std::make_pair("3He",1)]=-1.37;  EC_time_offset[std::make_pair("3He",2)]=-1.42; EC_time_offset[std::make_pair("3He",3)]=-1.55;
EC_time_offset[std::make_pair("3He",4)]=-1.53;  EC_time_offset[std::make_pair("3He",5)]=-1.49; EC_time_offset[std::make_pair("3He",6)]=-1.44;

EC_time_offset[std::make_pair("4He",1)]=0.72;  EC_time_offset[std::make_pair("4He",2)]=0.27; EC_time_offset[std::make_pair("4He",3)]=0.16;
EC_time_offset[std::make_pair("4He",4)]=0.21;  EC_time_offset[std::make_pair("4He",5)]=0.22; EC_time_offset[std::make_pair("4He",6)]=0.21;

EC_time_offset[std::make_pair("C12",1)]=0.50;  EC_time_offset[std::make_pair("C12",2)]=0.39; EC_time_offset[std::make_pair("C12",3)]=0.29;
EC_time_offset[std::make_pair("C12",4)]=0.29;  EC_time_offset[std::make_pair("C12",5)]=0.32; EC_time_offset[std::make_pair("C12",6)]=0.33;

EC_time_offset[std::make_pair("56Fe",1)]=0.75;  EC_time_offset[std::make_pair("56Fe",2)]=0.49; EC_time_offset[std::make_pair("56Fe",3)]=0.37;
EC_time_offset[std::make_pair("56Fe",4)]=0.39;  EC_time_offset[std::make_pair("56Fe",5)]=0.43; EC_time_offset[std::make_pair("56Fe",6)]=0.44;

elmom_corr_fact[0]=1.001;
elmom_corr_fact[1]=0.991;
elmom_corr_fact[2]=1.005;
elmom_corr_fact[3]=1.004;
elmom_corr_fact[4]=1.006;
elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
