#include <iostream>
#include <fstream>
#include <string>

bool input_reader(std::string input_loc,std::string* neutral_pes,std::string* cation_pes,std::string* neutral_dipole,std::string* cation_dipole,std::string* neutral_nac,std::string* ionization_coupling_file,std::string* out_file,std::string* read_file,std::string* wf_out_file,std::string* spectrum_out_file,std::string* mfpad_out_file,std::string* pi_cs_file,std::string* ionization_rate_file,std::string *dist_file,int* gsize_x,int* small_gsize_x,int* n_states_neut,int* n_states_cat,int* n_angles,int* n_k,int* k_mfpad,double* kmin,double* kmax,double* xmin,double* xmax,double* mass,double* total_time,double* h,double* efield_thresh,double* pot_vec_thresh,double* pump_strength,double* pump_origin,double* pump_sigma,double* pump_energy,double* pump_CEP,double* pprobe_delay,double* probe_strength,double* probe_sigma,double* probe_energy,double* probe_CEP )
{
   std::ifstream input;

   input.open(input_loc.c_str());
   if(!input.is_open())
   {
      std::cout<<"ERROR WHEN TRYING TO OPEN INPUT FILE "<<input_loc.c_str()<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
      exit(EXIT_SUCCESS);
   }
   else
   {
      std::cout<<"READING INPUT FILE"<<std::endl;
      input>>*neutral_pes;
      input>>*cation_pes;
      input>>*neutral_dipole;
      input>>*cation_dipole;
      input>>*neutral_nac;
      input>>*ionization_coupling_file;
      input>>*out_file;
      input>>*read_file;
      input>>*wf_out_file;
      input>>*spectrum_out_file;
      input>>*mfpad_out_file;
      input>>*pi_cs_file;
      input>>*ionization_rate_file;
      input>>*dist_file;
      input>>*gsize_x;
      input>>*small_gsize_x;
      input>>*n_states_neut;
      input>>*n_states_cat;
      input>>*n_angles;
      input>>*n_k;
      input>>*k_mfpad;
      input>>*kmin;
      input>>*kmax;
      input>>*xmin;
      input>>*xmax;
      input>>*mass;
      input>>*total_time;
      input>>*h;
      input>>*efield_thresh;
      input>>*pot_vec_thresh;
      input>>*pump_strength;
      input>>*pump_origin;
      input>>*pump_sigma;
      input>>*pump_energy;
      input>>*pump_CEP;
      input>>*pprobe_delay;
      input>>*probe_strength;
      input>>*probe_sigma;
      input>>*probe_energy;
      input>>*probe_CEP;

      std::cout<<"INPUT FILE CORRECTLY READ! "<<std::endl
      <<"neutral_pes :"<<*neutral_pes<<std::endl
      <<"cation_pes :"<<*cation_pes<<std::endl
      <<"neutral_dipole :"<<*neutral_dipole<<std::endl
      <<"cation_dipole :"<<*cation_dipole<<std::endl
      <<"neutral_nac :"<<*neutral_nac<<std::endl
      <<"ionization_coupling_file :"<<*ionization_coupling_file<<std::endl
      <<"out_file :"<<*out_file<<std::endl
      <<"read_file :"<<*read_file<<std::endl
      <<"wf_out_file :"<<*wf_out_file<<std::endl
      <<"spectrum_out_file :"<<*spectrum_out_file<<std::endl
      <<"mfpad_out_file :"<<*mfpad_out_file<<std::endl
      <<"ionization_rate_file"<<*ionization_rate_file<<std::endl
      <<"pi_cs_file :"<<*pi_cs_file<<std::endl
      <<"dist_file :"<<*dist_file<<std::endl
      <<"gsize_x :"<<*gsize_x<<std::endl
      <<"small_gsize_x :"<<*small_gsize_x<<std::endl
      <<"n_states_neut :"<<*n_states_neut<<std::endl
      <<"n_states_cat :"<<*n_states_cat<<std::endl
      <<"n_angles : "<<*n_angles<<std::endl
      <<"n_k :"<<*n_k<<std::endl
      <<"k_mfpad :"<<*k_mfpad<<std::endl
      <<"kmin :"<<*kmin<<std::endl
      <<"kmax :"<<*kmax<<std::endl
      <<"xmin :"<<*xmin<<std::endl
      <<"xmax :"<<*xmax<<std::endl
      <<"mass :"<<*mass<<std::endl
      <<"total_time :"<<*total_time<<std::endl
      <<"h :"<<*h<<std::endl;
      
      return 1;
   }
}
