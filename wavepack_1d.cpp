#include "wavepack_1d.hpp"

int main( int argc, char * argv [])
{
    using namespace std;

    std::string input_file_loc;
    std::string restart_file_loc;
    int init_time_index(0);
    int nproc(0);

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
       }
    }

    omp_set_num_threads(nproc);

    //PATH LINKING OTHER FILES
    string neutral_pes;
    string cation_pes;
    string neutral_dipole;
    string cation_dipole;
    string neutral_nac;
    string ionization_coupling_file;//LiH_PICE_R_i_j.txt
    string out_file;
    string read_file;
    string wf_out_file;
    string spectrum_out_file;
    string mfpad_out_file;
    string pi_cs_file;
    stringstream ss_wf;
    string s_wf;
    ofstream output;
//    ofstream read;
    ofstream wf_out;
    ofstream spectrum;
    ofstream mfpad;
    //string elec_wavepack_file="LiH_elec_wvpck.txt";
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

    input_reader(input_file_loc,&neutral_pes,&cation_pes,&neutral_dipole,&cation_dipole,&neutral_nac,&ionization_coupling_file,&out_file,&read_file,&wf_out_file,&spectrum_out_file,&mfpad_out_file,&pi_cs_file,&gsize_x,&small_gsize_x,&n_states_neut,&n_states_cat,&n_angles,&n_k,&kp,&kmin,&kmax,&xmin,&xmax,&mass,&total_time,&h,&efield_thresh,&pot_vec_thresh,&pump_strength,&pump_origin,&pump_sigma,&pump_energy,&pump_CEP,&pprobe_delay,&probe_strength,&probe_sigma,&probe_energy,&probe_CEP);

    int tgsize_x(gsize_x+10);
    int n_times(int(total_time/h));
    int dgsize(tgsize_x-gsize_x);

    wavefunction* Psi= new wavefunction(gsize_x,tgsize_x, n_states_neut,n_states_cat,n_angles*n_k);
    //hamilton_matrix* H=new hamilton_matrix(gsize_x,tgsize_x,small_gsize_x,n_states_neut,n_states_cat,n_k,n_angles,kmin,kmax,xmin,xmax,mass,n_times,h,pump_strength,pump_origin,pump_sigma,pump_energy,pump_CEP,probe_strength,pprobe_delay,probe_sigma,probe_energy,probe_CEP,efield_thresh,pot_vec_thresh,ionization_coupling_file.c_str());
    hamilton_matrix* H=new hamilton_matrix(gsize_x,tgsize_x,small_gsize_x,n_states_neut,n_states_cat,n_k,n_angles,kmin,kmax,xmin,xmax,mass,n_times,h,pump_strength,probe_strength,pump_origin,pprobe_delay,pump_sigma,probe_sigma,pump_energy,probe_energy,pump_CEP,probe_CEP,efield_thresh,pot_vec_thresh,ionization_coupling_file.c_str());



    H->set_pot_neut(neutral_pes.c_str());
    H->set_pot_cat(cation_pes.c_str());
    H->set_dm_neut(neutral_dipole.c_str());
    H->set_dm_cat(cation_dipole.c_str());
    H->set_NAC(neutral_nac);

    if(time_index != 0)
    {
        double *init_pot_vec=new double[3];
        if(fabs(H->pot_vec_mod()) >= H->pot_vec_thresh())
        {
           H->potential_vector(time_index,init_pot_vec);
           H->set_PICE(init_pot_vec);
        }
        else
        {
           H->set_PICE();
        }
        delete [] init_pot_vec;
    }
    else
    {
        H->set_PICE();
        output.open(out_file.c_str());
        output<<"Output from Wavepack_1D, developped by Stephan van den Wildenberg (Theoretical Physical Chemistry, University of Liege)"<<std::endl<<"File generated on "<<date_str<<std::endl;
        output<<"Pump and probe pulses have the following parameters:"<<std::endl
          <<"Strength: "<<pump_strength<<" (pump) ; "<<probe_strength<<" (probe) "<<std::endl
          <<"Sigma: "<<pump_sigma<<" (pump) ; "<<probe_sigma<<" (probe) "<<std::endl
          <<"Carrier Energy: "<<pump_energy<<" (pump) ; "<<probe_energy<<" (probe) "<<std::endl
          <<"CEP: "<<pump_CEP<<" (pump) ; "<<probe_CEP<<" (probe) "<<std::endl;
        output.close();
    }

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

//     read.open(read_file.c_str());
//    std::cout<<"Initial state:"<<std::endl;
//    for(int i=0;i!=tgsize_x;i++)
//    {
//        std::cout<<"probe "<<i<<std::endl;
//       read<<H->pot_neut(0,i)<<","<<Psi->show_neut_psi(i,0).real()<<","<<Psi->show_neut_psi(i,0).imag()<<std::endl;
//    }
//    read.close();
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


    while(time_index <= n_times)
    {
       propagate(Psi,H,&time_index,25);
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
       mfpad.open(mfpad_out_file.c_str(),ios_base::app);
       for(int o=0;o!=n_angles;o++)
       {
          temp=0;
          for(int r=0;r!=tgsize_x;r++)
          {
             temp+=std::norm(Psi->show_cat_psi(r,0,kp*n_angles+o));
          }
          mfpad<<time_index*h*0.02418884<<"   "<<H->k_spher_orient(0,o)<<"   "<<H->k_spher_orient(1,o)<<"   "<<temp<<std::endl;
       }mfpad<<std::endl;
       mfpad.close();

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

       Psi->save_wf(s_savefile.c_str());
    }

    delete Psi;
   // delete dPsi;
    delete H;


   return 0;
}
