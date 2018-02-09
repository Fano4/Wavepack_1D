#include "wavepack_1d.hpp"

int main( int argc, char * argv [])
{
    using namespace std;

    omp_set_num_threads(32);

    //PATH LINKING OTHER FILES
    string neutral_pes("/data1/home/stephan/Wavepack_1D/wavepack_int_input/Int_pes_");
    string cation_pes("/data1/home/stephan/LiH_gridtest/LiH_cat_");
    string neutral_dipole("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_neut_");
    string neutral_nac("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_NAC_");
    string ionization_coupling_file("/data1/home/stephan/LiH_512_points_pice/LiH_PICE_");//LiH_PICE_R_i_j.txt
//    string phase_file("/data1/home/stephan/LiH_gridtest/phase_");
    string out_file="wavepack_CEP0/Output.log";
    string read_file="wavepack_CEP0/gnu-out.txt";
    string wf_out_file="wavepack_CEP0/neut_wf_state_";
    string wf1d_out_file="wavepack_CEP0/neut_wf1d_state_";
    stringstream ss_wf;
    string s_wf;
    ofstream output;
    ofstream read;
    ofstream wf_out;
    //string elec_wavepack_file="LiH_elec_wvpck.txt";
    //PARAMETERS OF THE SIMULATION
    int gsize_x(512);
    int tgsize_x(552);
    int n_states_neut(15);
    int n_states_cat(0);
    int n_angles(0);
    int n_k(0);
    double xmin(0.8/0.529);//!!! THESE VALUES ARE IN ATOMIC UNITS AND NOT IN ANGSTROM
    double xmax(21.6/0.529);
    double mass(1836*(1.007825*6.015122795/(1.007825+6.015122795)));
    double total_time(100/0.02418884);
    double h(0.005/0.02418884);
    int n_times(int(total_time/h));
    int time_index(0);
    double dipole[3];
    double efield[3];
    double efield_thresh(0);

    wavefunction* Psi= new wavefunction(gsize_x,tgsize_x, n_states_neut,n_states_cat,n_angles*n_k);
    hamilton_matrix* H=new hamilton_matrix(gsize_x,tgsize_x,n_states_neut,n_states_cat,n_k,n_angles,xmin,xmax,mass,n_times,h,efield_thresh);


    H->set_pot_neut(neutral_pes.c_str());
    H->set_pot_cat(cation_pes.c_str());
    H->set_dm_neut(neutral_dipole.c_str());
    H->set_NAC(neutral_nac);
    H->set_PICE(ionization_coupling_file.c_str());

//    H->set_phase(phase_file.c_str());
    H->print_dipole_neut();

    Psi->initialize(H);

    read.open(read_file.c_str());
    for(int i=0;i!=gsize_x;i++)
    {
       read<<H->pot_neut(0,i)<<","<<Psi->show_neut_psi(i,0).real()<<","<<Psi->show_neut_psi(i,0).imag()<<std::endl;
    }
    read.close();
    std::cout<<"##################"<<std::endl;

    output.open(out_file.c_str());
    output<<"initial energy of the system "<<setprecision(15)<<H->energy(Psi,0)<<std::endl<<"initial norm of the system ";
    output<<setprecision(15)<<Psi->norm(H)<<std::endl;
    output.close();
       for(int m=0;m!=n_states_neut;m++)
       {
          ss_wf.str("");
          ss_wf<<wf_out_file<<m<<".wvpck";
          s_wf=ss_wf.str();
          wf_out.open(s_wf.c_str());
          wf_out.close();
          ss_wf.str("");
          ss_wf<<wf1d_out_file<<m<<".wvpck";
          s_wf=ss_wf.str();
          wf_out.open(s_wf.c_str());
          wf_out.close();
       }

    //wavefunction* dPsi= new wavefunction(gsize_x, n_states_neut,n_states_cat,n_angles*n_k);

    while(time_index <= n_times)
    {
       propagate(Psi,H,&time_index,2);
       H->electric_field(time_index,efield);
       output.open(out_file.c_str(),ios_base::app);

       output<<"Time "<<time_index*h*0.02418884<<std::endl;
       output<<"Electric field (X)"<<efield[0]<<std::endl;
       output<<"Electric field (Y)"<<efield[1]<<std::endl;
       output<<"Electric field (Z)"<<efield[2]<<std::endl;
       output<<"Energy of the system "<<setprecision(15)<<H->energy(Psi,time_index)<<std::endl;
       Psi->show_dipole(dipole);
       output<<"Total dipole "<<setprecision(15)<<dipole[0]<<", "<<dipole[1]<<", "<<dipole[2]<<std::endl;
       output<<"Norm of the system "<<setprecision(15)<<Psi->norm(H)<<std::endl;

       for(int m=0;m!=n_states_neut;m++)
       {
          output<<"Population on neutral state "<<m+1<<" = "<<setprecision(15)<<Psi->state_pop(0,m)<<std::endl;
          ss_wf.str("");
          ss_wf<<wf_out_file<<m<<".wvpck";
          s_wf=ss_wf.str();
          wf_out.open(s_wf.c_str(),ios_base::app);
          for(int k=0;k!=gsize_x;k++)
          {
              wf_out<<time_index*h*0.02418884<<"   "<<xmin+k*(xmax-xmin)/gsize_x<<"   "<<real(Psi->show_neut_psi(k,m))<<"   "<<imag(Psi->show_neut_psi(k,m))<<std::endl;
          }wf_out<<std::endl;
          wf_out.close();

          ss_wf.str("");
          ss_wf<<wf1d_out_file<<m<<".wvpck";
          s_wf=ss_wf.str();
          wf_out.open(s_wf.c_str(),ios_base::app);
          for(int k=0;k!=gsize_x;k++)
          {
              wf_out<<time_index*h*0.02418884<<"   "<<xmin+k*(xmax-xmin)/gsize_x<<"   "<<H->pot_neut(m,k)<<"   "<<real(Psi->show_neut_psi(k,m))<<"   "<<imag(Psi->show_neut_psi(k,m))<<std::endl;
          }wf_out<<std::endl<<std::endl;
          wf_out.close();

       }
       for(int m=0;m!=n_states_cat;m++)
       {
          output<<"Population on cation state "<<m+1<<" = "<<setprecision(15)<<Psi->state_pop(1,m,H)<<std::endl;
       }
       output.close();
    }

    delete Psi;
   // delete dPsi;
    delete H;


   return 0;
}
