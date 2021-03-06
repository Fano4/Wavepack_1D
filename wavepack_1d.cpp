#include "wavepack_1d.hpp"

int main( int argc, char * argv [])
{
    using namespace std;

    std::string input_file_loc;
    std::string restart_file_loc;
    int init_time_index(0);
    int nproc(0);
    bool restart_bool(0);

    if(argc<3)
    {
       std::cout<<"Error: too few arguments for calling Wavepack. I need the location of the input file and the number of processes used. EXIT"<<std::endl;
       exit(EXIT_SUCCESS);
    }
    else if (argc > 3 && argc < 5)
    {
       std::cout<<"Error : Too few arguments for restart. I need the localization of the wavepack file and the inital time index. EXIT"<<std::endl;
       exit(EXIT_SUCCESS);
    }
    else
    {
       input_file_loc=std::string(argv[1]);
       nproc=std::atoi(std::string(argv[2]).c_str());
       std::cout<<"Running wavepack using input file "<<input_file_loc.c_str()<<" with "<<nproc<<"processes !"<<std::endl;
       if(argc > 3)
       {
          restart_file_loc=std::string(argv[3]);
          init_time_index=std::atoi(std::string(argv[4]).c_str());
          std::cout<<"Restarting wave packet at time_index "<<init_time_index<<" from input file "<<restart_file_loc.c_str()<<std::endl;
          restart_bool=1;
       }
    }

    omp_set_num_threads(nproc);

    //PATH LINKING OTHER FILES
    string neutral_pes;
    string cation_pes;
    string neutral_dipole;
    string cation_dipole;
    string neutral_nac;
    string neutral_spinorb;
    string ionization_coupling_file;//LiH_PICE_R_i_j.txt
    string out_file;
    string read_file;
    string wf_out_file;
    string spectrum_out_file;
    string mfpad_out_file;
    string ionization_rate_file;
    string pi_cs_file;
    stringstream ss_wf;
    string s_wf;
    string dist_file;
    string average_mom_file;
    ofstream output;
    ofstream read;
    ofstream wf_out;
    ofstream spectrum;
    ofstream mfpad;
//    string elec_wavepack_file="/data1/home/stephan/wvpck_recollision_IR_10fs/LiH_elec_wvpck.txt";
    //PARAMETERS OF THE SIMULATION
    int gsize_x;
    int small_gsize_x;
    int n_states_neut;
    int n_states_cat;
    int n_angles;
    int n_k;
    int kp;
    double kmin;
    double kmax;//
    double xmin;//!!! THESE VALUES ARE IN ATOMIC UNITS AND NOT IN ANGSTROM
    double xmax;
    double mass;
    double total_time;
    double h;
    int time_index(init_time_index);
    double dipole[3];
    double efield[3];
    double efield_thresh;
    double pot_vec_thresh;
    double norm;
    double temp;
    double pump_strength;
    double probe_strength;
    double pump_origin;
    double pprobe_delay;
    double pump_sigma;
    double probe_sigma;
    double pump_energy;
    double probe_energy;
    double pump_CEP;
    double probe_CEP;
    time_t date_num;
    std::string date_str;
    
    time(&date_num);
    date_str=ctime(&date_num);

    std::stringstream ss_savefile;
    std::string s_savefile;
    char* tempstr=new char[6];
    mkstemp(tempstr);
    std::string temp_string(tempstr);
    ss_savefile.str("");
    ss_savefile<<wf_out_file<<"_"<<temp_string<<".svwvpck";    
    s_savefile=ss_savefile.str();

   double respectrum(0);
   double imspectrum(0);

    input_reader(input_file_loc,&neutral_pes,&cation_pes,&neutral_dipole,&cation_dipole,&neutral_nac,&neutral_spinorb,&ionization_coupling_file,&out_file,&read_file,&wf_out_file,&spectrum_out_file,&mfpad_out_file,&pi_cs_file,&ionization_rate_file,&average_mom_file,&dist_file,&gsize_x,&small_gsize_x,&n_states_neut,&n_states_cat,&n_angles,&n_k,&kp,&kmin,&kmax,&xmin,&xmax,&mass,&total_time,&h,&efield_thresh,&pot_vec_thresh,&pump_strength,&pump_origin,&pump_sigma,&pump_energy,&pump_CEP,&pprobe_delay,&probe_strength,&probe_sigma,&probe_energy,&probe_CEP);

    int tgsize_x(small_gsize_x+6);
    int n_times(int(total_time/h));
    int dgsize(tgsize_x-small_gsize_x);

    wavefunction* Psi= new wavefunction(gsize_x,tgsize_x, n_states_neut,n_states_cat,n_angles*n_k);
    //hamilton_matrix* H=new hamilton_matrix(gsize_x,tgsize_x,small_gsize_x,n_states_neut,n_states_cat,n_k,n_angles,kmin,kmax,xmin,xmax,mass,n_times,h,pump_strength,pump_origin,pump_sigma,pump_energy,pump_CEP,probe_strength,pprobe_delay,probe_sigma,probe_energy,probe_CEP,efield_thresh,pot_vec_thresh,ionization_coupling_file.c_str());
    hamilton_matrix* H=new hamilton_matrix(gsize_x,tgsize_x,small_gsize_x,n_states_neut,n_states_cat,n_k,n_angles,kmin,kmax,xmin,xmax,mass,n_times,h,pump_strength,probe_strength,pump_origin,pprobe_delay,pump_sigma,probe_sigma,pump_energy,probe_energy,pump_CEP,probe_CEP,efield_thresh,pot_vec_thresh,ionization_coupling_file.c_str());

    H->set_pice_mapping();

    H->set_pot_neut(neutral_pes.c_str());
    H->set_pot_cat(cation_pes.c_str());
    H->set_dm_neut(neutral_dipole.c_str());
    H->set_dm_cat(cation_dipole.c_str());
    H->set_NAC(neutral_nac);
    H->set_spinorb_neut(neutral_spinorb);
    

    /*
    for(int i=0;i!=tgsize_x;i++)
    {
       std::cout<<xmin+(i-dgsize)*(xmax-xmin)/small_gsize_x<<"    ";
       for(int n=0;n!=n_states_neut;n++)
       {
          std::cout<<H->pot_neut(n,i)<<"    ";
       }
       for(int n=0;n!=n_states_cat;n++)
       {
          std::cout<<H->pot_cat(n,i)<<"    ";
       }
       std::cout<<std::endl;
    }
    exit(EXIT_SUCCESS);
    */
    //SET UP ARRAYS FOR OUTPUTING IONIZATION/RECOMBINATION RATE
    double**ionization_rate=new double*[n_states_neut];
    for(int i=0;i!=n_states_neut;i++)
    {
       ionization_rate[i]=new double[n_states_cat];
    }

    if(time_index != 0)
    {
       if(n_states_cat != 0)
       {
          double *init_pot_vec=new double[3];
          if(fabs(H->pot_vec_mod()) >= H->pot_vec_thresh())
          {
             H->potential_vector(time_index,init_pot_vec);
             H->sphere_dist_read(dist_file);
             H->set_PICE();
          }
          else
          {
             H->sphere_dist_read(dist_file);
             H->set_PICE();
          }
          delete [] init_pot_vec;
       }
    }
    else
    {
       if(n_states_cat != 0 && n_angles != 0)
       {
          H->sphere_dist_read(dist_file);
          H->sphere_dist_save(dist_file);
          H->set_PICE();
       }
        output.open(out_file.c_str());
        output<<"Output from Wavepack_1D, developped by Stephan van den Wildenberg (Theoretical Physical Chemistry, University of Liege)"<<std::endl<<"File generated on "<<date_str<<std::endl;
        output<<"Reading input from "<<input_file_loc.c_str()<<std::endl;

        if(restart_bool)
           output<<"Restarting wave packet at time_index "<<init_time_index<<" from input file "<<restart_file_loc.c_str()<<std::endl;

        output<<"Pump and probe pulses have the following parameters:"<<std::endl
          <<"Strength: "<<pump_strength<<" (pump) ; "<<probe_strength<<" (probe) "<<std::endl
          <<"Sigma: "<<pump_sigma<<" (pump) ; "<<probe_sigma<<" (probe) "<<std::endl
          <<"Carrier Energy: "<<pump_energy<<" (pump) ; "<<probe_energy<<" (probe) "<<std::endl
          <<"CEP: "<<pump_CEP<<" (pump) ; "<<probe_CEP<<" (probe) "<<std::endl;
        output.close();
    }


    stringstream tempsstr;
    string filename;
    for(int i=0;i!=n_states_neut;i++)
    {
       tempsstr.str("");
       tempsstr<<pi_cs_file.c_str()<<i<<"_"<<0<<".txt";
       filename=tempsstr.str();
       H->plot_integrated_cross_section(filename.c_str(),i,0);
    }
    

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
/*
    ofstream test;
    test.open("pes.txt");
    for(int i=0;i!=tgsize_x;i++)
    {
       for(int n=0;n!=n_states_neut;n++)
       {
          test<<H->pot_neut(n,i)<<"    ";
       }
       for(int n=0;n!=n_states_cat;n++)
       {
          test<<H->pot_cat(n,i)<<"    ";
       }
       test<<std::endl;
    }
       test.close();
       exit(EXIT_SUCCESS);
 */      
    if(time_index == 0)
    {
        Psi->initialize(H);
    }
    else
    {
        Psi->initialize(H);
        Psi->load_wf(restart_file_loc.c_str());
        std::cout<<"Wave function restarted from checkpoint file!"<<std::endl;
    }


       H->diagonalize_Hamilton();
       double *eigenmat=new double[n_states_neut*tgsize_x*tgsize_x*n_states_neut];
       H->eigenstates_matrix(0,eigenmat);
       /*
       ofstream eigenout;
       eigenout.open("/data1/home/stephan/test_diago_hamilt/LiH_eigenvectors_40Ang.txt");
       for(int n=0;n!=n_states_neut*tgsize_x;n++)
       {
          for(int m=0;m!=n_states_neut*tgsize_x;m++)
          {
             eigenout<<eigenmat[m*n_states_neut*tgsize_x+n]<<",";

          }eigenout<<std::endl;
       }
       eigenout.close();
       eigenout.open("/data1/home/stephan/test_diago_hamilt/LiH_eigenvalues_40Ang.txt");
       for(int n=0;n!=n_states_neut;n++)
       {
          for(int g=0;g!=tgsize_x;g++)
          {
             eigenout<<H->eigenvalue_neut(n,g)<<std::endl;
          }
       }
       eigenout.close();
       exit(EXIT_SUCCESS);
       */
      // Psi->projection_eigenstates(1);
      
/*
       for(int n=0;n!=n_states_neut;n++)
       {
          for( int g=0;g!=tgsize_x;g++)
          {
             std::cout<<eigenval[n*tgsize_x+g]<<"  "<<real(proj_state->show_neut_psi(g,n))<<"  "<<imag(proj_state->show_neut_psi(g,n))<<"  "<<pow(abs(proj_state->show_neut_psi(g,n)),2)<<std::endl;
          }
       }*/
       /*
       double *dipole_mat=new double[n_states_neut*tgsize_x*tgsize_x*n_states_neut];
       H->change_basis_dipole(dipole_mat);
       ofstream dipole_output;
       dipole_output.open("/data1/home/stephan/test_diago_hamilt/dipole_mat_eigenbasis_256.txt");
       for(int n=0;n!=tgsize_x*n_states_neut;n++)
       {
          for(int k=0;k!=tgsize_x*n_states_neut;k++)
          {
             dipole_output<<dipole_mat[n*tgsize_x*n_states_neut+k]<<",";
          }dipole_output<<endl;
       }

       dipole_output.close();
       exit(EXIT_SUCCESS);
*/

       read.open(ionization_rate_file.c_str());
       read.close();

       read.open(average_mom_file.c_str());
       read.close();

       double average_mom_x(0);
       double average_mom_y(0);
       double average_mom_z(0);
       /*
       read.open("/data1/home/stephan/wavepack_photoelectron_020519/photoelec_Z.txt");//!!! YOU HAVE TO REPLACE THIS WITH A NON-CONSTANT USER DEFINED STRING
       read.close();
       const int ncx=2;
       const int ncy=2;
       const int ncz=100;
       const double cxmin(-9);
       const double cxmax(9);
       const double cymin(-9);
       const double cymax(9);
       const double czmin(-30);
       const double czmax(30);
       const double dx((cxmax-cxmin)/ncx);
       const double dy((cymax-cymin)/ncy);
       const double dz((czmax-czmin)/ncz);
       double z(0);
       double *cube_photoelec_dens=new double [ncx*ncy*ncz];
       double dens_sum(0);
      */
//     read.open(read_file.c_str());
//    std::cout<<"Initial state:"<<std::endl;
//    for(int i=0;i!=tgsize_x;i++)
//    {
//        std::cout<<"probe "<<i<<std::endl;
//       read<<H->pot_neut(0,i)<<","<<Psi->show_neut_psi(i,0).real()<<","<<Psi->show_neut_psi(i,0).imag()<<std::endl;
//    }
    //read.close();
    std::cout<<"##################"<<std::endl;

    if(time_index == 0)
    {
       output.open(out_file.c_str(),ios_base::app);
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
    }


    bool test(0);
    while(time_index <= n_times)
    {
       if( test )
       {
           propagate(Psi,H,&time_index,250);
       }
       test=1;

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
              wf_out<<time_index*h*0.02418884<<"   "<<xmin+(k-dgsize)*(xmax-xmin)/small_gsize_x<<"   "<<real(Psi->show_neut_psi(k,m))<<"   "<<imag(Psi->show_neut_psi(k,m))<<std::endl;
          }wf_out<<std::endl;
          wf_out.close();
          
       }
       /*read.open(read_file.c_str(),ios_base::app);
       if(!read.is_open())
       {
         output<<"!!!!!!!!!!!! UNABLE TO WRITE IN "<<read_file.c_str()<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
         exit(EXIT_FAILURE);
       }
       */
       spectrum.open(spectrum_out_file.c_str(),ios_base::app);
       if(!spectrum.is_open())
       {
         output<<"!!!!!!!!!!!! UNABLE TO WRITE IN "<<read_file.c_str()<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
         exit(EXIT_FAILURE);
       }
      average_mom_x=0;
      average_mom_y=0;
      average_mom_z=0;
      H->potential_vector(time_index,efield);
       for(int k=0;k!=n_k;k++)
       {
          temp=0;
          for(int c=0;c!=n_states_cat;c++)
          {
             for(int o=0;o!=n_angles;o++)
             {
                for(int r=0;r!=tgsize_x;r++)
                {
//                   if(Psi->state_pop(1,0,H) != 0)
                   {
                      average_mom_x+=(H->k_mod_val(k)*sin(H->k_spher_orient(0,o))*cos(H->k_spher_orient(1,o))+efield[0])*std::norm(Psi->show_cat_psi(r,c,k*n_angles+o));
                      average_mom_y+=(H->k_mod_val(k)*sin(H->k_spher_orient(0,o))*sin(H->k_spher_orient(1,o))+efield[1])*std::norm(Psi->show_cat_psi(r,c,k*n_angles+o));
                      average_mom_z+=(H->k_mod_val(k)*cos(H->k_spher_orient(0,o))+efield[2])*std::norm(Psi->show_cat_psi(r,c,k*n_angles+o));
                   }
        //           read<<time_index*h*0.02418884<<"   "<<H->k_mod_val(k)<<"    "<<H->k_spher_orient(0,o)<<"    "<<H->k_spher_orient(1,o)<<"    "<<Psi->show_cat_psi(r,c,k*n_angles+o)<<std::endl;
                   temp+=std::norm(Psi->show_cat_psi(r,c,k*n_angles+o));
                }
             }
          }
          spectrum<<time_index*h*0.02418884<<"   "<<H->k_mod_val(k)<<"   "<<temp<<std::endl;
       }//std::cout<<std::endl;
          spectrum<<std::endl;
      // read.close();
       spectrum.close();

       read.open(average_mom_file.c_str(),ios_base::app);
       read<<time_index*h*0.02418884<<"   "<<average_mom_x<<"   "<<average_mom_y<<"   "<<average_mom_z<<std::endl;
       read.close();
       
       for(int m=0;m!=n_states_cat;m++)
       {
          temp=Psi->state_pop(1,m,H);
          norm+=temp;
          output<<"Population on cation state "<<m+1<<" = "<<setprecision(15)<<temp<<std::endl;
       }
       output<<"Norm of the system "<<setprecision(15)<<norm<<std::endl;
       output.close();

       H->PI_rate(time_index,ionization_rate,Psi);
       read.open(ionization_rate_file.c_str(),ios_base::app);
       if(!read.is_open())
       {
         output<<"!!!!!!!!!!!! UNABLE TO WRITE IN "<<ionization_rate_file.c_str()<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
         exit(EXIT_FAILURE);
       }
       read<<time_index*h*0.02418884<<" ";
       for(int i=0;i!=n_states_neut;i++)
       {
          for(int j=0;j!=n_states_cat;j++)
          {
             read<<ionization_rate[i][j]<<" ";
          }
       }read<<std::endl;
       read.close();

       Psi->save_wf(s_savefile.c_str());
    }

       mfpad.open(mfpad_out_file.c_str(),ios_base::app);
       for(kp=0;kp!=n_k;kp++)
       {
          for(int o=0;o!=n_angles;o++)
          {
             temp=0;
             for(int r=0;r!=tgsize_x;r++)
             {
                temp+=std::norm(Psi->show_cat_psi(r,0,kp*n_angles+o));
             }
             mfpad<<H->k_mod_val(kp)<<"   "<<H->k_spher_orient(0,o)<<"   "<<H->k_spher_orient(1,o)<<"   "<<temp<<std::endl;
          }mfpad<<std::endl;
       }
       mfpad.close();


    delete Psi;
   // delete dPsi;
    delete H;


   return 0;
}
