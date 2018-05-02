#include "wavepack_1d.hpp"

int main( int argc, char * argv [])
{
    using namespace std;

    omp_set_num_threads(16);

    //PATH LINKING OTHER FILES
    string neutral_pes("/data1/home/stephan/LiH_gridtest_++2df2p/Int_pes_");
    string cation_pes("/data1/home/stephan/LiH_gridtest_++2df2p/Int_pes_cat_");
    string neutral_dipole("/data1/home/stephan/LiH_gridtest_++2df2p/LiH_neut_");
    string cation_dipole("/data1/home/stephan/LiH_gridtest_++2df2p/LiH_cat_");
    string neutral_nac("/data1/home/stephan/LiH_gridtest_++2df2p/LiH_NAC_");
    string ionization_coupling_file("/data1/home/stephan/LiH_PICE_++2df2p/LiH_");//LiH_PICE_R_i_j.txt
//    string phase_file("/data1/home/stephan/LiH_gridtest/phase_");
    string out_file="wvpck_res_2df_astridpulse_CEPPI/Output.log";
    string read_file="wvpck_res_2df_astridpulse_CEPPI/PI_spectrum.txt";
    string wf_out_file="wvpck_res_2df_astridpulse_CEPPI/neut_wf_state_";
    string wf1d_out_file="wvpck_res_2df_astridpulse_CEPPI/neut_wf1d_state_";
    stringstream ss_wf;
    string s_wf;
    ofstream output;
    ofstream read;
    ofstream wf_out;
    //string elec_wavepack_file="LiH_elec_wvpck.txt";
    //PARAMETERS OF THE SIMULATION
    int gsize_x(512);
    int tgsize_x(gsize_x+10);
    int small_gsize_x(64);
    int dgsize(tgsize_x-gsize_x);
    int n_states_neut(14);//15);
    int n_states_cat(1);//1);
    int n_angles(64);//128);
    int n_k(75);//100);
    double xmin(0.8/0.529);//!!! THESE VALUES ARE IN ATOMIC UNITS AND NOT IN ANGSTROM
    double xmax(21.6/0.529);
    double mass(1836*(1.007825*6.015122795/(1.007825+6.015122795)));
    double total_time(100/0.02418884);
    double h(0.0025/0.02418884);
    int n_times(int(total_time/h));
    int time_index(0);
    double dipole[3];
    double efield[3];
    double efield_thresh(1e-5);
    double norm;
    double temp;

    wavefunction* Psi= new wavefunction(gsize_x,tgsize_x, n_states_neut,n_states_cat,n_angles*n_k);
    hamilton_matrix* H=new hamilton_matrix(gsize_x,tgsize_x,small_gsize_x,n_states_neut,n_states_cat,n_k,n_angles,xmin,xmax,mass,n_times,h,efield_thresh);


    H->set_pot_neut(neutral_pes.c_str());
    H->set_pot_cat(cation_pes.c_str());
    H->set_dm_neut(neutral_dipole.c_str());
    H->set_dm_cat(cation_dipole.c_str());
    H->set_NAC(neutral_nac);
    H->set_PICE(ionization_coupling_file.c_str());

/*    int i(10);
    int j(tgsize_x*1+12);
    int state_index_1;
    int state_index_2;
    int state_index_cont_1;
    int state_index_cont_2;
    int grid_index_1;
    int grid_index_2;
    H->show_indexes(i,j,&state_index_1,&grid_index_1,&state_index_cont_1,&state_index_2,&grid_index_2,&state_index_cont_2);
    std::cout<<"NAC "<<state_index_1<<"-"<<state_index_2<<" at "<<grid_index_1<<","<<grid_index_2<<"="<<H->show_nac(state_index_1,state_index_2,grid_index_1,grid_index_2)<<std::endl;
 
    H->show_indexes(i,j,&state_index_1,&grid_index_1,&state_index_cont_1,&state_index_2,&grid_index_2,&state_index_cont_2);
    std::cout<<state_index_1<<"-"<<state_index_2<<" at "<<grid_index_1<<","<<grid_index_2<<std::endl;
    
    std::cout<<"H element "<<i<<","<<j<<"="<<H->hamilt_element(0,i,j)<<std::endl;
    std::cout<<"H element "<<j<<","<<i<<"="<<H->hamilt_element(0,j,i)<<std::endl;
    return 0;
   */ 
    
//    for(int i=0;i!=H->n_k();i++)
//       H->dk(i);

//    exit(EXIT_SUCCESS);
//    H->set_phase(phase_file.c_str());
//    H->print_dipole_neut();

/*    ofstream test;
    test.open("test.txt");
    for(int i=0;i!=tgsize_x;i++)
    {
       for(int n=0;n!=n_states_neut;n++)
       {
          test<<H->pot_neut(n,i)<<"    ";
       }test<<std::endl;
    }
       test.close();
       exit(EXIT_SUCCESS);
       */
    Psi->initialize(H);

    read.open(read_file.c_str());
//    std::cout<<"Initial state:"<<std::endl;
    for(int i=0;i!=tgsize_x;i++)
    {
//        std::cout<<"probe "<<i<<std::endl;
       read<<H->pot_neut(0,i)<<","<<Psi->show_neut_psi(i,0).real()<<","<<Psi->show_neut_psi(i,0).imag()<<std::endl;
    }
    read.close();
    std::cout<<"##################"<<std::endl;

    output.open(out_file.c_str());
    output<<"initial energy of the system "<<setprecision(15)<<H->energy(Psi,0)<<std::endl<<"initial norm of the system = 1"<<std::endl;;
//    output<<setprecision(15)<<Psi->norm(H)<<std::endl;
    output.close();
       for(int m=0;m!=n_states_neut;m++)
       {
          ss_wf.str("");
          ss_wf<<wf_out_file<<m<<".wvpck";
          s_wf=ss_wf.str();
          wf_out.open(s_wf.c_str());
          wf_out.close();
//          ss_wf.str("");
//          ss_wf<<wf1d_out_file<<m<<".wvpck";
//          s_wf=ss_wf.str();
//          wf_out.open(s_wf.c_str());
//          wf_out.close();
       }

    //wavefunction* dPsi= new wavefunction(gsize_x, n_states_neut,n_states_cat,n_angles*n_k);

    while(time_index <= n_times)
    {
       propagate(Psi,H,&time_index,16);
       H->electric_field(time_index,efield);
       output.open(out_file.c_str(),ios_base::app);

       output<<"Time "<<time_index*h*0.02418884<<std::endl;
       output<<"Electric field (X)"<<efield[0]<<std::endl;
       output<<"Electric field (Y)"<<efield[1]<<std::endl;
       output<<"Electric field (Z)"<<efield[2]<<std::endl;
       output<<"Energy of the system "<<setprecision(15)<<H->energy(Psi,time_index)<<std::endl;
       Psi->show_dipole(dipole);
       output<<"Total dipole "<<setprecision(15)<<dipole[0]<<", "<<dipole[1]<<", "<<dipole[2]<<std::endl;

       norm=0;
       for(int m=0;m!=n_states_neut;m++)
       {
          temp=Psi->state_pop(0,m);
          norm+=temp;
          output<<"Population on neutral state "<<m+1<<" = "<<setprecision(15)<<temp<<std::endl;
          ss_wf.str("");
          ss_wf<<wf_out_file<<m<<".wvpck";
          s_wf=ss_wf.str();
          wf_out.open(s_wf.c_str(),ios_base::app);
          for(int k=dgsize;k!=tgsize_x;k++)
          {
              wf_out<<time_index*h*0.02418884<<"   "<<xmin+(k-dgsize)*(xmax-xmin)/gsize_x<<"   "<<real(Psi->show_neut_psi(k,m))<<"   "<<imag(Psi->show_neut_psi(k,m))<<std::endl;
          }wf_out<<std::endl;
          wf_out.close();

/*          ss_wf.str("");
          ss_wf<<wf1d_out_file<<m<<".wvpck";
          s_wf=ss_wf.str();
          wf_out.open(s_wf.c_str(),ios_base::app);
          for(int k=0;k!=gsize_x;k++)
          {
              wf_out<<time_index*h*0.02418884<<"   "<<xmin+k*(xmax-xmin)/gsize_x<<"   "<<H->pot_neut(m,k)<<"   "<<real(Psi->show_neut_psi(k,m))<<"   "<<imag(Psi->show_neut_psi(k,m))<<std::endl;
          }wf_out<<std::endl<<std::endl;
          wf_out.close();
*/
       }
/*       double spectrum(0);
       for(int k=0;k!=n_k;k++)
       {
          spectrum=0;
          for(int r=0;r!=tgsize_x;r++)
          {
             for(int o=0;o!=n_angles;o++)
             {
                spectrum+=std::norm(Psi->show_cat_psi(r,0,k*n_angles+o));
             }
          }
       }*/
       for(int m=0;m!=n_states_cat;m++)
       {
          temp=Psi->state_pop(1,m,H);
          norm+=temp;
          output<<"Population on cation state "<<m+1<<" = "<<setprecision(15)<<temp<<std::endl;
       }
       output<<"Norm of the system "<<setprecision(15)<<norm<<std::endl;
       output.close();
    }

    delete Psi;
   // delete dPsi;
    delete H;


   return 0;
}
