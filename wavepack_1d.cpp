#include "wavepack_1d.hpp"

int main( int argc, char * argv [])
{
    using namespace std;

    omp_set_num_threads(32);

    //PATH LINKING OTHER FILES
    string neutral_pes("/data1/home/stephan/Wavepack_1D/wavepack_int_input/Int_pes_");
    string cation_pes("/data1/home/stephan/Wavepack_1D/wavepack_int_input/Int_pes_cat_");
    string neutral_dipole("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_neut_");
    string cation_dipole("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_cat_");
    string neutral_nac("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_NAC_");
    string ionization_coupling_file("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_pice_data_newnorm.h5");//LiH_PICE_R_i_j.txt
//    string phase_file("/data1/home/stephan/LiH_gridtest/phase_");
    string out_file="/data1/home/stephan/wvpck_mfpad_probe/Output_test_probe_1.log";
    string read_file="/data1/home/stephan/wvpck_mfpad_probe/read_test_probe_1.txt";//="/data1/home/stephan/wvpck_newpice/read_03_10_18.txt";
    string wf_out_file="/data1/home/stephan/wvpck_mfpad_probe/neut_test_probe_1_wf_state_";//"/data1/home/stephan/wvpck_newpice/neut_03_10_18_wf_state_";
    string spectrum_out_file="/data1/home/stephan/wvpck_mfpad_probe/spect_test_probe_1.txt";//"/data1/home/stephan/wvpck_newpice/neut_03_10_18_wf_state_";
    string mfpad_out_file="/data1/home/stephan/wvpck_mfpad_probe/mfpad_test_probe_1.txt";//"/data1/home/stephan/wvpck_newpice/neut_03_10_18_wf_state_";
//    string wf1d_out_file="wvpck_test5/neut_wf1d_state_";
    string pi_cs_file="/data1/home/stephan/wvpck_mfpad_probe/cs_";
    stringstream ss_wf;
    string s_wf;
    ofstream output;
//    ofstream read;
    ofstream wf_out;
    ofstream spectrum;
    ofstream mfpad;
    //string elec_wavepack_file="LiH_elec_wvpck.txt";
    //PARAMETERS OF THE SIMULATION
    int gsize_x(512);
    int tgsize_x(gsize_x+10);
    int small_gsize_x(64);
    int dgsize(tgsize_x-gsize_x);
    int n_states_neut(10);
    int n_states_cat(1);
    int n_angles(45);
    int n_k(30);
    double kp(13);
    double kmin=1e-4;
    double kmax=1.0;//
    double xmin(0.8/0.529);//!!! THESE VALUES ARE IN ATOMIC UNITS AND NOT IN ANGSTROM
    double xmax(21.6/0.529);
    double mass(1836*(1.007825*6.015122795/(1.007825+6.015122795)));
    double total_time(100/0.02418884);
    double h(0.0005/0.02418884);
    int n_times(int(total_time/h));
    int time_index(0);
    double dipole[3];
    double efield[3];
    double efield_thresh(1e-4);
    double pot_vec_thresh(1.);
    double norm;
    double temp;
   double respectrum(0);
   double imspectrum(0);

    wavefunction* Psi= new wavefunction(gsize_x,tgsize_x, n_states_neut,n_states_cat,n_angles*n_k);
    hamilton_matrix* H=new hamilton_matrix(gsize_x,tgsize_x,small_gsize_x,n_states_neut,n_states_cat,n_k,n_angles,kmin,kmax,xmin,xmax,mass,n_times,h,efield_thresh,pot_vec_thresh,ionization_coupling_file.c_str());



    H->set_pot_neut(neutral_pes.c_str());
    H->set_pot_cat(cation_pes.c_str());
    H->set_dm_neut(neutral_dipole.c_str());
    H->set_dm_cat(cation_dipole.c_str());
    H->set_NAC(neutral_nac);
    H->set_PICE();
/*
    stringstream tempstr;
    string filename;
    for(int i=0;i!=n_states_neut;i++)
    {
       tempstr.str("");
       tempstr<<pi_cs_file.c_str()<<i<<"_"<<0<<".txt";
       filename=tempstr.str();
       H->plot_integrated_cross_section(filename.c_str(),i,0);
    }
    
    return EXIT_SUCCESS;
*/
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

//     read.open(read_file.c_str());
//    std::cout<<"Initial state:"<<std::endl;
//    for(int i=0;i!=tgsize_x;i++)
//    {
//        std::cout<<"probe "<<i<<std::endl;
//       read<<H->pot_neut(0,i)<<","<<Psi->show_neut_psi(i,0).real()<<","<<Psi->show_neut_psi(i,0).imag()<<std::endl;
//    }
//    read.close();
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
       propagate(Psi,H,&time_index,100);
       H->electric_field(time_index,efield);
       Psi->set_dipole(H);
       output.open(out_file.c_str(),ios_base::app);

       output<<"Time "<<time_index*h*0.02418884<<std::endl;
       output<<"Electric field (X)"<<efield[0]<<std::endl;
       output<<"Electric field (Y)"<<efield[1]<<std::endl;
       output<<"Electric field (Z)"<<efield[2]<<std::endl;
       output<<"Energy of the system "<<setprecision(15)<<H->energy(Psi,time_index)<<std::endl;
       Psi->show_dipole(dipole,0);
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
          if(!wf_out.is_open())
          {
             output<<"!!!!!!!!!!!! UNABLE TO WRITE IN "<<s_wf.c_str()<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
             exit(EXIT_FAILURE);
          }
          
          for(int k=dgsize;k!=tgsize_x;k++)
          {
              wf_out<<time_index*h*0.02418884<<"   "<<xmin+(k-dgsize)*(xmax-xmin)/gsize_x<<"   "<<real(Psi->show_neut_psi(k,m))<<"   "<<imag(Psi->show_neut_psi(k,m))<<std::endl;
          }wf_out<<std::endl;
          wf_out.close();
          
       }
//read.open(read_file.c_str(),ios_base::app);
       /*if(!read.is_open())
       {
         output<<"!!!!!!!!!!!! UNABLE TO WRITE IN "<<read_file.c_str()<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
         exit(EXIT_FAILURE);
       }*/
       //mfpad.open(mfpad_out_file.c_str(),ios_base::app);
       for(int o=0;o!=n_angles;o++)
       {
          temp=0;
          for(int r=0;r!=tgsize_x;r++)
          {
             temp+=std::norm(Psi->show_cat_psi(r,0,kp*n_angles+o));
          }
//          mfpad<<time_index*h*0.02418884<<"   "<<H->k_spher_orient(0,o)<<"   "<<H->k_spher_orient(1,o)<<"   "<<temp<<std::endl;
       }//mfpad<<std::endl;
       //mfpad.close();

       spectrum.open(spectrum_out_file.c_str(),ios_base::app);
       if(!spectrum.is_open())
       {
         output<<"!!!!!!!!!!!! UNABLE TO WRITE IN "<<read_file.c_str()<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
         exit(EXIT_FAILURE);
       }
       for(int k=0;k!=n_k;k++)
       {
          temp=0;
          for(int c=0;c!=n_states_cat;c++)
          {
             for(int o=0;o!=n_angles;o++)
             {
                for(int r=0;r!=tgsize_x;r++)
                {
      //             read<<time_index*h*0.02418884<<"   "<<H->k_mod_val(k)<<"    "<<H->k_spher_orient(0,o)<<"    "<<H->k_spher_orient(1,o)<<"    "<<Psi->show_cat_psi(r,c,k*n_angles+o)<<"    "<<Psi->show_cat_psi(r,c,k*n_angles+o)<<std::endl;
                   temp+=std::norm(Psi->show_cat_psi(r,c,k*n_angles+o));
                }
             }
          }
          spectrum<<time_index*h*0.02418884<<"   "<<H->k_mod_val(k)<<"   "<<temp<<std::endl;
       }//std::cout<<std::endl;
          spectrum<<std::endl;
    //   read.close();
       spectrum.close();
       
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
