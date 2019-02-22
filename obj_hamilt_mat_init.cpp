

// !!! THERE IS A MANUAL SIGN CORRECTION ON THE PICE !! IT NEEDS AUTOMATION TO BE USED WITH ANOTHER BASIS SET !!!


//THIS FILE CONTAINS THE ROUTINES ASSOCIATED TO THE INITIALIZATION AND SETTINGS OF THE HAMILTON MATRIX OBJECT.
//
//TABLE OF CONTENT:
//
//HAMILTON MATRIX OBJECT INITIALIZATION AND SETTINGS
//constructor=hamilton_matrix(int gsize_x,int n_states_neut,int n_states_cat,int n_k,int n_angles,double xmin,double xmax,int n_times,double h)
//destructor=~hamilton_matrix()
//set_pot_neut(std::string file_address)
//set_pot_cat(std::string file_address)
//set_dm_neut(std::string file_address)
//set_dm_cat(std::string file_address)
//set_PICE(std::string file_address)
//set_NAC(std::string file_address)
//set_phase(std::string file_address)
//void hamilton_matrix::print_dipole_neut()
//double hamilton_matrix::show_dm_neut(int state_index_1,int state_index_2,int grid_index,int component)
//void hamilton_matrix::sphere_dist_read(std::string dist_file_s)
//void hamilton_matrix::sphere_dist_save(std::string dist_file_s)
//void hamilton_matrix::sphere_dist_gen(bool randiso,int n_phi)
//double hamilton_matrix::k_mod_val(int k_index) const
//double hamilton_matrix::show_dm_cat(int state_index_1,int state_index_2,int grid_index,int component)
//double hamilton_matrix::k_spher_orient(bool component, int index) const
//void hamilton_matrix::set_pice_mapping_size();

//HAMILTON MATRIX OBJECT INITIALIZATION AND SETTINGS
//##########################################################################
//
//Constructor of the Hamilton_matrix object. This initializes all the arrays and values relative to the Hamiltonian operator
//
//##########################################################################
hamilton_matrix::hamilton_matrix(int gsize_x,int tgsize_x,int small_gsize_x,int n_states_neut,int n_states_cat,int n_k,int n_angles,double kmin,double kmax,double xmin,double xmax,double mass,int n_times,double h,double pump_strength,double probe_strength,double pump_origin,double pprobe_delay,double pump_sigma,double probe_sigma,double pump_energy,double probe_energy,double pump_CEP,double probe_CEP,double efield_thresh,double pot_vec_thresh,std::string pice_data_loc)
{
   using namespace std;
   //initialize grid parameters and time settings
   this->m_gsize_x=gsize_x;
   this->m_tgsize_x=tgsize_x;
   this->m_small_gsize_x=small_gsize_x;
   this->m_n_states_neut=n_states_neut;
   this->m_n_states_cat=n_states_cat;
   this->m_n_k=n_k;
   this->m_n_angles=n_angles;
   this->m_n_states_cont=n_k*n_angles;
   this->m_xmin=xmin;
   this->m_xmax=xmax;
   this->m_n_times=n_times;
   this->m_h=h;
   this->m_mass=mass;
   this->m_efield_thresh=efield_thresh;
   this->m_pot_vec_thresh=pot_vec_thresh;
   this->m_kmin=kmin;
   this->m_kmax=kmax;
   double delta_x((xmax-xmin)/small_gsize_x);
   std::cout<<"Size of the pixel is "<<delta_x<<std::endl;
   //
   //initialize pulses parameters

   this->m_pump_strength=pump_strength;
   this->m_probe_strength=probe_strength;
   this->m_pump_origin=pump_origin;
   this->m_pprobe_delay=pprobe_delay;
   this->m_pump_sigma=pump_sigma;
   this->m_probe_sigma=probe_sigma;
   this->m_pump_energy=pump_energy;
   this->m_probe_energy=probe_energy;
   this->m_pump_CEP=pump_CEP;
   this->m_probe_CEP=probe_CEP;

   //TEMPORARY PICE SIGN CORRECTION 
   //
   this->sign_corr=new double *[this->m_n_states_neut];
   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      this->sign_corr[i]=new double[this->m_small_gsize_x];
      for(int j=0;j!=this->m_small_gsize_x;j++)
      {
         this->sign_corr[i][j]=1;
      }
   }
   //
   //initialize potential energy surfaces arrays
   //
   std::cout<<"initializing PES arrays...";
   this->m_pot_neut=new double [tgsize_x*n_states_neut];
   std::cout<<" ...";
   this->m_pot_cat=new double [tgsize_x*n_states_cat*this->m_n_states_cont];
   std::cout<<"PES arrays initialized!"<<std::endl;
   //
   //initialize dipole moment surfaces arrays
   //of the neutral
   //
   std::cout<<"initializing dipole arrays...";
   this->m_dmx_neut=new double*[n_states_neut*n_states_neut];
   std::cout<<" ...";
   this->m_dmy_neut=new double*[n_states_neut*n_states_neut];
   std::cout<<" ...";
   this->m_dmz_neut=new double*[n_states_neut*n_states_neut];
   std::cout<<" ... .";
   for(int i=0;i!=n_states_neut*n_states_neut;i++)
   {
      this->m_dmx_neut[i]=new double[tgsize_x];
      this->m_dmy_neut[i]=new double[tgsize_x];
      this->m_dmz_neut[i]=new double[tgsize_x];
   }
   std::cout<<"dipole arrays of the neutral initialized!"<<std::endl;
   //
   //initialize dipole moment surfaces arrays
   //of the cation
   //
   this->m_dmx_cat=new double*[n_states_cat*n_states_cat];
   std::cout<<" ...";
   this->m_dmy_cat=new double*[n_states_cat*n_states_cat];
   std::cout<<" ...";
   this->m_dmz_cat=new double*[n_states_cat*n_states_cat];
   std::cout<<" ... .";
   for(int i=0;i!=n_states_cat*n_states_cat;i++)
   {
      this->m_dmx_cat[i]=new double[tgsize_x];
      this->m_dmy_cat[i]=new double[tgsize_x];
      this->m_dmz_cat[i]=new double[tgsize_x];
   }
   std::cout<<"dipole arrays of the cation initialized!"<<std::endl;
   //
   //initialize photoionization coupling elements surfaces arrays
   //
   std::cout<<"initializing PICE arrays...";

   this->pice_data=new pice_set(pice_data_loc); //Initialize a pice dataset object

   this->m_dk_vec=new double[this->m_n_states_cont];
   this->m_PICE_x=new std::complex<double> *[n_states_neut*n_states_cat*this->m_n_states_cont];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_cat*this->m_n_states_cont;i++)
   {
      this->m_PICE_x[i]=new std::complex<double>[tgsize_x];
   }
   std::cout<<" ...";
   this->m_PICE_y=new std::complex<double> *[n_states_neut*n_states_cat*this->m_n_states_cont];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_cat*this->m_n_states_cont;i++)
   {
      this->m_PICE_y[i]=new std::complex<double>[tgsize_x];
   }
   std::cout<<" ...";
   this->m_PICE_z=new std::complex<double> *[n_states_neut*n_states_cat*this->m_n_states_cont];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_cat*this->m_n_states_cont;i++)
   {
      this->m_PICE_z[i]=new std::complex<double>[tgsize_x];
   }
   std::cout<<" ...";
   std::cout<<"PICE arrays initialized!"<<std::endl;
   std::cout<<"initializing momentum vectors arrays...";
   this->k_orientation=new double*[2];
   std::cout<<" ...";
   this->k_orientation[0]=new double[n_angles];//[0] is theta
   std::cout<<" ...";
   this->k_orientation[1]=new double[n_angles];//[1] is phi
   std::cout<<" ...";
   this->k_modulus=new double[n_k];
   this->sphere_dist_gen(1);
   std::cout<<" ...";
   this->pice_time_mapping=new int[this->m_n_times];

   //
   //Initialize the coordinates of the reciprocal space for describing the continuum
   //
   for(int i=0;i!=this->m_n_k;i++)
   {
      this->k_modulus[i]=kmin+i*(kmax-kmin)/this->m_n_k;
      for(int j=0;j!=this->m_n_angles;j++)
      {
         //rho(k) * dk * dO = k*k * dk * dO
          m_dk_vec[i*this->m_n_angles+j]=4*acos(-1)*pow(this->k_modulus[i],2)*(kmax-kmin)/(this->m_n_k*this->m_n_angles);
      }
   }
   std::cout<<"momentum vectors arrays initialized!"<<std::endl;
   //
   //initialize Non-adiabatic coupling surfaces arrays
   //
   std::cout<<"initializing NAC arrays...";
   this->m_NAC=new double*[n_states_neut*n_states_neut];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_neut;i++)
   {
      this->m_NAC[i]=new double [tgsize_x];
   }
   std::cout<<" ...";
   std::cout<<" NAC initialized!"<<std::endl;
   //
   //Initialize and set up kinetic energy matrix
   //
   std::cout<<"initializing kinetic energy array...";
   this->kinetic_energy=new double[tgsize_x*tgsize_x];
   this->derivative_matrix=new double[tgsize_x*tgsize_x];
   //FINITE DIFFERENCE METHOD IMPLEMENTATION OF KINETIC ENERGY (ABRAMOWITZ EQ. 25.3.24)
   for(int i=0;i!=tgsize_x;i++)
   {
      for(int j=0;j!=tgsize_x;j++)
      {
         this->kinetic_energy[i*tgsize_x+j]=0;
         if(i==j+2)//(i,i-2)
         {
            this->kinetic_energy[i*tgsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(-1);
         }
         else if(i==j+1)//(i,i-1)
         {
            this->kinetic_energy[i*tgsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(16);
         }
         else if(i==j)//(i,i)
         {
            this->kinetic_energy[i*tgsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(-30);
         }
         else if(i==j-1)//(i,i+1)
         {
            this->kinetic_energy[i*tgsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(16);
         }
         else if(i==j-2)//(i,i+2)
         {
            this->kinetic_energy[i*tgsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(-1);
         }
         else
            this->kinetic_energy[i*tgsize_x+j]=0;
      }
   }
   std::cout<<"kinetic energy array initialized!"<<std::endl;
   //FINITE DIFFERENCE METHOD IMPLEMENTATION OF THE DERIVATIVE MATRIX
   for(int i=0;i!=tgsize_x;i++)
   {
      for(int j=0;j!=tgsize_x;j++)
      {
         this->derivative_matrix[i*tgsize_x+j]=0;
         if(i==j+2)//(i,i-2)
         {
            this->derivative_matrix[i*tgsize_x+j]=(1/(12*delta_x))*(1);
         }
         else if(i==j+1)//(i,i-1)
         {
            this->derivative_matrix[i*tgsize_x+j]=(1/(12*delta_x))*(-8);
         }
         else if(i==j-1)//(i,i+1)
         {
            this->derivative_matrix[i*tgsize_x+j]=(1/(12*delta_x))*(8);
         }
         else if(i==j-2)//(i,i+2)
         {
            this->derivative_matrix[i*tgsize_x+j]=(1/(12*delta_x))*(-1);
         }
         else
            this->derivative_matrix[i*tgsize_x+j]=0;
      }
   }
   std::cout<<"derivative matrix array initialized!"<<std::endl;

}
//##########################################################################
//
//Destructor of the Hamilton_matrix object
//
//##########################################################################
hamilton_matrix::~hamilton_matrix()
{
   delete [] this->m_pot_neut;
   delete [] this->m_pot_cat;
   delete [] this->m_dmx_neut;
   delete [] this->m_dmy_neut;
   delete [] this->m_dmz_neut;
   delete [] this->m_dmx_cat;
   delete [] this->m_dmy_cat;
   delete [] this->m_dmz_cat;
   delete [] this->m_PICE_x;
   delete [] this->m_PICE_y;
   delete [] this->m_PICE_z;
   delete [] this->m_NAC;
   delete [] this->kinetic_energy;
}
//##########################################################################
//
// This routine sets up the Potential energy curves for the neutral electronic states from an input file located at "file_address"
// The potential energy surfaces are extrapolated with a harmonic border so that the potential bordr is repulsive and that the wave packet never collides with the edges
//
//
//##########################################################################
void hamilton_matrix::set_pot_neut(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   double var(0);
   int dgsize (this->m_tgsize_x-this->m_small_gsize_x);

   ifstream input_file;

   cout<<"gathering PES of the neutral"<<endl;

   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      name_indenter.str("");
      name_indenter<<file_address<<i+1<<".input";
      filename=name_indenter.str();
      input_file.open(filename.c_str());
      
      cout<<"state "<<i+1<<"  ########################## in "<<filename.c_str()<<endl;
      if(input_file.is_open())
      {
         input_file.seekg(0);

         input_file>>var;
         for(int j=dgsize;j!=this->m_tgsize_x;j++)
         {
            this->m_pot_neut[this->m_tgsize_x*i+j]=var;
            for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
            {
               input_file>>var;
            }
         }
         input_file.close();
         for(int j=1;j<=dgsize;j++)
         {
            //Here, the PEC is extrapolated with a harmonic potential.
            this->m_pot_neut[this->m_tgsize_x*i+(dgsize-j)]=this->m_pot_neut[this->m_tgsize_x*i+dgsize]+j*(this->m_pot_neut[this->m_tgsize_x*i+dgsize]-this->m_pot_neut[this->m_tgsize_x*i+dgsize+1])+5e-3*j*j;
         }
      }
      else
      {
         cout<<"ERROR POTENTIAL FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
         exit(EXIT_FAILURE);
      }
   }
}
//##########################################################################
//
// This routine sets up the Potential energy curves for the cation electronic states from an input file located at "file_address"
// The potential energy surfaces are extrapolated with a harmonic border so that the potential bordr is repulsive and that the wave packet never collides with the edges
//
//##########################################################################
void hamilton_matrix::set_pot_cat(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   double var(0);
   int dgsize (this->m_tgsize_x-this->m_small_gsize_x);

   ifstream input_file;

   cout<<"gathering PES of the cation"<<endl;

   for(int i=0;i!=m_n_states_cat;i++)
   {
      name_indenter.str("");
      name_indenter<<file_address<<i+1<<".input";
      filename=name_indenter.str();
      input_file.open(filename.c_str());

      cout<<"state "<<i+1<<"########################## in "<<filename.c_str()<<endl;

      if(input_file.is_open())
      {
         input_file>>var;
         for(int j=dgsize;j!=this->m_tgsize_x;j++)
         {
            input_file>>this->m_pot_cat[this->m_tgsize_x*i+j];
            for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
            {
               input_file>>var;
            }
         }
         input_file.close();
         for(int j=1;j<=dgsize;j++)
         {
            //Here, the PEC is extrapolated with a harmonic potential.
            this->m_pot_cat[this->m_tgsize_x*i+(dgsize-j)]=this->m_pot_cat[this->m_tgsize_x*i+dgsize]+j*(this->m_pot_cat[this->m_tgsize_x*i+dgsize]-this->m_pot_cat[this->m_tgsize_x*i+dgsize+1])+5e-3*j*j;
         }
      }
      else
      {
         cout<<"ERROR POTENTIAL FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
         exit(EXIT_FAILURE);
      }
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_dm_neut(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   double temp;
   int dgsize(this->m_tgsize_x-this->m_small_gsize_x);

   ifstream input_file;
   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      for(int j=i;j!=this->m_n_states_neut;j++)
      {   
            name_indenter.str("");
            name_indenter<<file_address<<"DMX"<<"_"<<j+1<<"_"<<i+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
                input_file>>temp;
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    this->m_dmx_neut[i*this->m_n_states_neut+j][k]=temp;
                    this->m_dmx_neut[j*this->m_n_states_neut+i][k] = this->m_dmx_neut[i*this->m_n_states_neut+j][k];
                    for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
                    {
                       input_file>>temp;
                    }
                }
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmx_neut[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmx_neut[i*this->m_n_states_neut+j][dgsize];
                    this->m_dmx_neut[j*this->m_n_states_neut+i][dgsize-k]=this->m_dmx_neut[i*this->m_n_states_neut+j][dgsize-k];
                }
                input_file.close();
            }
      else
      {
         cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
         exit;
      }
            
            name_indenter.str("");
            name_indenter<<file_address<<"DMY"<<"_"<<j+1<<"_"<<i+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
               input_file>>temp;
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    this->m_dmy_neut[i*this->m_n_states_neut+j][k]=temp;
                    this->m_dmy_neut[j*this->m_n_states_neut+i][k]=this->m_dmy_neut[i*this->m_n_states_neut+j][k];
                    for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
                    {
                       input_file>>temp;
                    }
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmy_neut[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmy_neut[i*this->m_n_states_neut+j][dgsize];
                    this->m_dmy_neut[j*this->m_n_states_neut+i][dgsize-k]=this->m_dmy_neut[i*this->m_n_states_neut+j][dgsize-k];
                }
            }
            else
            {
               cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit;
            }
            
            name_indenter.str("");
            name_indenter<<file_address<<"DMZ"<<"_"<<j+1<<"_"<<i+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
               input_file>>temp;
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    this->m_dmz_neut[i*this->m_n_states_neut+j][k]=temp;
                    this->m_dmz_neut[j*this->m_n_states_neut+i][k]=this->m_dmz_neut[i*this->m_n_states_neut+j][k];
                    for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
                    {
                       input_file>>temp;
                    }
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmz_neut[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmz_neut[i*this->m_n_states_neut+j][dgsize];
                    this->m_dmz_neut[j*this->m_n_states_neut+i][dgsize-k]=this->m_dmz_neut[i*this->m_n_states_neut+j][dgsize-k];
                }
            }
            else
            {
               cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit;
            }
      }
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_dm_cat(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   double temp(0);
   int dgsize(this->m_tgsize_x-this->m_small_gsize_x);

   ifstream input_file;
   for(int i=0;i!=this->m_n_states_cat;i++)
   {
      for(int j=i;j!=this->m_n_states_cat;j++)
      {   
            name_indenter.str("");
            name_indenter<<file_address<<"DMX"<<"_"<<i+1<<"_"<<j+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
               input_file>>temp;
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    this->m_dmx_cat[i*this->m_n_states_cat+j][k]=temp;
                    this->m_dmx_cat[j*this->m_n_states_cat+i][k]=temp;
                    for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
                    {
                       input_file>>temp;
                    }
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmx_cat[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmx_cat[i*this->m_n_states_neut+j][dgsize];
                    this->m_dmx_cat[j*this->m_n_states_neut+i][dgsize-k]=this->m_dmx_cat[i*this->m_n_states_neut+j][dgsize-k];
                }
            }
            else
            {
               cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit(EXIT_FAILURE);
            }
            
            name_indenter.str("");
            name_indenter<<file_address<<"DMY"<<"_"<<i+1<<"_"<<j+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
               input_file>>temp;
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    this->m_dmy_cat[i*this->m_n_states_cat+j][k]=temp;
                    this->m_dmy_cat[j*this->m_n_states_cat+i][k]=temp;
                    for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
                    {
                       input_file>>temp;
                    }
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmy_cat[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmy_cat[i*this->m_n_states_neut+j][dgsize];
                    this->m_dmy_cat[j*this->m_n_states_neut+i][dgsize-k]=this->m_dmy_cat[i*this->m_n_states_neut+j][dgsize-k];
                }
            }
            else
            {
               cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit(EXIT_FAILURE);
            }
            
            name_indenter.str("");
            name_indenter<<file_address<<"DMZ"<<"_"<<i+1<<"_"<<j+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
               input_file>>temp;
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    this->m_dmz_cat[i*this->m_n_states_cat+j][k]=temp;
                    this->m_dmz_cat[j*this->m_n_states_cat+i][k]=temp;
                    for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
                    {
                       input_file>>temp;
                    }
                }
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmz_cat[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmz_cat[i*this->m_n_states_neut+j][dgsize];
                    this->m_dmz_cat[j*this->m_n_states_neut+i][dgsize-k]=this->m_dmz_cat[i*this->m_n_states_neut+j][dgsize-k];
                }

                input_file.close();
            }
            else
            {
               cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit(EXIT_FAILURE);
            }
      }
   }
}

//##########################################################################
//
//##########################################################################

void hamilton_matrix::set_PICE()
{
   const int nstneut=this->m_n_states_neut;
   const int nstcat=this->m_n_states_cat;
   const int nk=this->m_n_k;
   const int ndist=this->m_n_angles;
   std::stringstream ss_file_cs;
   std::string s_file_cs;
   int dgsize_x(this->m_tgsize_x-this->m_small_gsize_x);

   std::complex<double> pice_x;
   std::complex<double> pice_y;
   std::complex<double> pice_z;

   double *pot_vec=new double[3];
   using namespace std;
   double duration;
   clock_t begin;
   clock_t end;

   int ratio(this->m_gsize_x/this->m_small_gsize_x);

   int k(0);
   int l(0);
   int x(0);
   int i(0);
   int j(0);
   for(int t=0;t!=this->mapping_size;t++)
   {
      pot_vec[0]=this->pot_vec_reduced_mapping[t][0];
      pot_vec[1]=this->pot_vec_reduced_mapping[t][1];
      pot_vec[2]=this->pot_vec_reduced_mapping[t][2];

       for( x=0;x<this->m_small_gsize_x;x++)
       {
          std::cout<<x<<" position "<<0.529*(this->m_xmin+x*(this->m_xmax-this->m_xmin)/this->m_small_gsize_x)<<"Angstrom"<<std::endl;
          for( i=0;i<nstneut;i++)
          {
             std::cout<<"stneut"<<i<<std::endl;
             for( j=0;j<nstcat;j++)
             { 
                 std::cout<<"stcat"<<j<<std::endl;
                 begin = clock();
                 #pragma omp parallel for private(k,l) shared(x,i,j,nk,ndist,nstcat,dgsize_x,pot_vec,ratio)
                 for(k=0;k<nk;k++)
                 {
                     for(l=0;l<ndist;l++)
                     {
                        this->pice_data->fill_pice(&this->m_PICE_sto_x[t][i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x+dgsize_x],&this->m_PICE_sto_y[t][i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x+dgsize_x],&this->m_PICE_sto_z[t][i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x+dgsize_x],int(ratio*x),i,j,this->k_orientation[0][l],this->k_orientation[1][l],this->k_modulus[k],pot_vec);
                     }
                  }
                  end = clock();
                  duration=double(end-begin)/CLOCKS_PER_SEC;
                  std::cout<<duration/32<<"s"<<std::endl;

             }
          }
       }
   }
      /*
   for( x=0;x<this->m_small_gsize_x;x++)
   {
      std::cout<<x<<" position "<<0.529*(this->m_xmin+x*(this->m_xmax-this->m_xmin)/this->m_small_gsize_x)<<"Angstrom"<<std::endl;
      for( i=0;i<nstneut;i++)
      {
         std::cout<<"stneut"<<i<<std::endl;
         for( j=0;j<nstcat;j++)
         {
         std::cout<<"stcat"<<j<<std::endl;
         begin = clock();
            #pragma omp parallel for private(k,l) shared(x,i,j,nk,ndist,nstcat,dgsize_x,pot_vec,ratio)
            for(k=0;k<nk;k++)
            {
               for(l=0;l<ndist;l++)
               {
                  this->pice_data->fill_pice(&this->m_PICE_x[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x+dgsize_x],&this->m_PICE_y[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x+dgsize_x],&this->m_PICE_z[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x+dgsize_x],int(ratio*x),i,j,this->k_orientation[0][l],this->k_orientation[1][l],this->k_modulus[k],pot_vec);
               }
            }
            end = clock();
           duration=double(end-begin)/CLOCKS_PER_SEC;
           std::cout<<duration/32<<"s"<<std::endl;

         }
      }
   }
   */
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_NAC(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   double temp(0);
   int dgsize(this->m_tgsize_x-this->m_small_gsize_x);

   ifstream input_file;
   for(int i=0;i!=m_n_states_neut;i++)
   {
      for(int j=0;j!=m_n_states_neut;j++)
      {
         for(int k=0;k!=this->m_tgsize_x;k++)
         {
            this->m_NAC[j*this->m_n_states_neut+i][k]=0;
         }
      }
   }

   for(int i=0;i!=m_n_states_neut;i++)
   {
      for(int j=i;j!=m_n_states_neut;j++)
      {   
            name_indenter.str("");
            name_indenter<<file_address<<j+1<<"_"<<i+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file>>temp;
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                   this->m_NAC[i*this->m_n_states_neut+j][k]=temp;
                   for(int t=0;t!=this->m_gsize_x/this->m_small_gsize_x-1;t++)
                   {
                      input_file>>temp;
                   }


//                   this->m_NAC[i*this->m_n_states_neut+j][k]*=5e6; /// !!!!!!REMOVE THIS LINE !!!!!

//                   std::cout<<"NACME "<<i<<"-"<<j<<" at "<<k<<"="<<this->m_NAC[i*this->m_n_states_neut+j][k]<<std::endl;
                   this->m_NAC[j*this->m_n_states_neut+i][k]=-this->m_NAC[i*this->m_n_states_neut+j][k];
                }
                input_file.close();
                for (int k=0;k!=dgsize;k++)
                {
                   this->m_NAC[i*this->m_n_states_neut+j][k]=this->m_NAC[i*this->m_n_states_neut+j][dgsize];
                   this->m_NAC[j*this->m_n_states_neut+i][k]=-this->m_NAC[i*this->m_n_states_neut+j][k];
                }
            }
            else
            {
               cout<<"ERROR NAC FILE NOT FOUND:"<<filename.c_str()<<endl<<"SETTING ALL NAC TO ZERO"<<endl;
               //exit(EXIT_FAILURE);
            }
      }
   }

}
double hamilton_matrix::kinetic_energy_matrix(int i,int j)
{
   return this->kinetic_energy[this->m_tgsize_x*i+j];
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::pot_neut(int state_index,int grid_index)
{
//   std::cout<<"probe pot print "<<state_index*(this->m_tgsize_x)+grid_index<<std::endl;
   return this->m_pot_neut[state_index*(this->m_tgsize_x)+grid_index];
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::pot_cat(int state_index,int grid_index)
{
   return this->m_pot_cat[state_index*(this->m_tgsize_x)+grid_index];
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::rescale_pot(double min_pot)
{
   for(int m=0;m!=this->m_n_states_neut;m++)
   {
      for(int g=0;g!=this->m_tgsize_x;g++)
      {
         this->m_pot_neut[m*(this->m_tgsize_x)+g]-=min_pot;
      }   
   }
   for(int m=0;m!=this->m_n_states_cat;m++)
   {
      for(int g=0;g!=this->m_tgsize_x;g++)
      {
         this->m_pot_cat[m*(this->m_tgsize_x)+g]-=min_pot;
      }   
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_phase(std::string file_address)
{
   using namespace std;
   ifstream input_file1;
   ifstream input_file2;
   stringstream ss_infile1;
   stringstream ss_infile2;
   string s_infile1;
   string s_infile2;

   int factor1;
   int factor2;
   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      for(int j=i;j!=this->m_n_states_neut;j++)
      {
         ss_infile1.str("");
         ss_infile1<<file_address.c_str();
         ss_infile1<<i+1<<".input";
         s_infile1=ss_infile1.str();
         ss_infile2.str("");
         ss_infile2<<file_address.c_str();
         ss_infile2<<j+1<<".input";
         s_infile2=ss_infile2.str();

         input_file1.open(s_infile1.c_str());
         input_file2.open(s_infile2.c_str());
         if(!input_file1.is_open() || !input_file2.is_open())
         {
            cout<<"ERROR WHILE OPENING PHASE FACTOR FILE: ";
               if(!input_file1.is_open())
                  cout<<s_infile1.c_str()<<std::endl;
               else
                  cout<<s_infile2.c_str()<<std::endl;
            exit(EXIT_FAILURE);
         }

         for(int k=0;k!=this->m_gsize_x;k++)
         {
            input_file1>>factor1;
            input_file2>>factor2;
            this->m_dmx_neut[j*this->m_n_states_neut+i][k]*=factor1*factor2;
            this->m_dmy_neut[j*this->m_n_states_neut+i][k]*=factor1*factor2;
            this->m_dmz_neut[j*this->m_n_states_neut+i][k]*=factor1*factor2;
            this->m_dmx_neut[i*this->m_n_states_neut+j][k]*=factor1*factor2;
            this->m_dmy_neut[i*this->m_n_states_neut+j][k]*=factor1*factor2;
            this->m_dmz_neut[i*this->m_n_states_neut+j][k]*=factor1*factor2;
         }
      }
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::print_dipole_neut()
{
   using namespace std;
   string sdipole_x_out("dipole_neut_x.txt");
   string sdipole_y_out("dipole_neut_y.txt");
   string sdipole_z_out("dipole_neut_z.txt");
   ofstream dipole_x;
   ofstream dipole_y;
   ofstream dipole_z;
   dipole_x.open(sdipole_x_out.c_str());
   dipole_y.open(sdipole_y_out.c_str());
   dipole_z.open(sdipole_z_out.c_str());

   for(int k=0;k!=this->m_gsize_x;k++)
   {
      for(int i=0;i!=this->m_n_states_neut;i++)
      {
         for(int j=i;j!=this->m_n_states_neut;j++)
         {
            dipole_x<<setw(12)<<setprecision(8)<<this->m_dmx_neut[i*this->m_n_states_neut+j][k];
            dipole_y<<setw(12)<<setprecision(8)<<this->m_dmy_neut[i*this->m_n_states_neut+j][k];
            dipole_z<<setw(12)<<setprecision(8)<<this->m_dmz_neut[i*this->m_n_states_neut+j][k];
         }
      }
      dipole_x<<endl;
      dipole_y<<endl;
      dipole_z<<endl;
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::sphere_dist_save(std::string dist_file_s)
{
   using namespace std;
   ofstream dist_file;
   int n_points_sphere((this->m_n_states_cont)/this->m_n_k);

   dist_file.open(dist_file_s.c_str(),ios_base::trunc);

   if(!dist_file.is_open())
   {
      std::cout<<"ERROR WHILE OPENING ANGULAR DISTRIBUTION FILE"<<std::endl<<dist_file_s.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int i=0;i!=n_points_sphere;i++)
      {
         dist_file<<setprecision(15)<<this->k_orientation[0][i]<<"  "<<this->k_orientation[1][i]<<std::endl;
      }
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::sphere_dist_read(std::string dist_file_s)
{
   using namespace std;
   ifstream dist_file;
   int n_points_sphere((this->m_n_states_cont)/this->m_n_k);

   dist_file.open(dist_file_s.c_str());

   if(!dist_file.is_open())
   {
      std::cout<<"ERROR WHILE OPENING ANGULAR DISTRIBUTION FILE"<<std::endl<<dist_file_s.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int i=0;i!=n_points_sphere;i++)
      {
         dist_file>>this->k_orientation[0][i];
         dist_file>>this->k_orientation[1][i];
      }
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::sphere_dist_gen(bool randiso,int n_phi)
{
   if(this->m_n_states_cont == 0 || this->m_n_k ==0)
      return ;
   double random;
   double random2;
   double random3;
   double temp;
   double n_theta(0);
   double theta(0);
   double phi(0);
   double Pi=acos(-1);
   int n_points_sphere((this->m_n_states_cont)/this->m_n_k);
   srand(time(0));

   if(randiso)
   {
      for(int i=0;i!=n_points_sphere;i++)
      {
         random=double(rand()%2000)-1000;
         random2=double(rand()%2000)-1000;
         random3=double(rand()%2000)-1000;
         temp=sqrt(pow(random,2)+pow(random2,2)+pow(random3,2));
         random/=temp;
         random2/=temp;
         random3/=temp;
         this->k_orientation[0][i]=acos(random3);
         if(random>=0 && random2>=0)
            this->k_orientation[1][i]=atan(random2/random);

         else if(random<0 && random2>=0)
            this->k_orientation[1][i]=Pi+atan(random2/random);

         else if(random<0 && random2<0)
            this->k_orientation[1][i]=Pi+atan(random2/random);

         else if(random>0 && random2 <0)
         this->k_orientation[1][i]=2*Pi+atan(random2/random);

         //DEBOGAGE
//         std::cout<<this->k_orientation[0][i]<<"   "<<this->k_orientation[1][i]<<std::endl;
      }
//      exit(EXIT_SUCCESS);
   }
   else
   {
      if(n_phi != 0)
         n_theta=n_points_sphere/n_phi;
      else
      {
         std::cout<<"CANNOT GENERATE A REGULAR SPHERICAL DISTRIBUTION WITH ZERO AZIMUTHAL ANGLE"<<std::endl;
         exit(EXIT_FAILURE);
      }
      for(int i=0;i!=n_theta;i++)
      {
         theta=i*acos(-1)/n_theta;
         for(int j=0;j!=n_phi;j++)
         {
            phi=j*2*acos(-1)/n_phi;
            this->k_orientation[0][i*n_phi+j]=theta;
            this->k_orientation[1][i*n_phi+j]=phi;
         }
      }
   }
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::k_mod_val(int k_index) const
{
   return this->k_modulus[k_index];
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::show_dm_neut(int state_index_1,int state_index_2,int grid_index,int component)
{
   switch(component)
   {
      case 0:
         return this->m_dmx_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index];
      case 1:
         return this->m_dmy_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index];
      case 2:
         return this->m_dmz_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index];
      default:
         std::cout<<"ERROR COMPONENT OF DIPOLE MOMENT NOT RECOGNIZED IN HAMILTON_MATRIX::SHOW_DM_NEUT === > EXIT"<<std::endl;
         exit(EXIT_FAILURE);
   }
   return 0;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::show_dm_cat(int state_index_1,int state_index_2,int grid_index,int component)
{
   switch(component)
   {
      case 0:
         return this->m_dmx_cat[state_index_1*this->m_n_states_neut+state_index_2][grid_index];
      case 1:
         return this->m_dmy_cat[state_index_1*this->m_n_states_neut+state_index_2][grid_index];
      case 2:
         return this->m_dmz_cat[state_index_1*this->m_n_states_neut+state_index_2][grid_index];
      default:
         std::cout<<"ERROR COMPONENT OF DIPOLE MOMENT NOT RECOGNIZED IN HAMILTON_MATRIX::SHOW_DM_CAT === > EXIT"<<std::endl;
         exit(EXIT_FAILURE);
   }
   return 0;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::k_spher_orient(bool component, int index) const
{
   return this->k_orientation[component][index];
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_pice_mapping()
{
//!!!!!!!!!!!!!!WORKS ONLY FOR ELECTRIC FIELD POLARIZED ALONG Z !!!!!!!!!!!!!!!
   double vector[3];
   double mod(0);
   int ratio(0);
   bool test(1);
   int *ratio_list=new int[this->m_n_times];
   std::vector<int> reduced_ratio_list;
   std::vector<double> reduced_potvec_list_x;
   std::vector<double> reduced_potvec_list_y;
   std::vector<double> reduced_potvec_list_z;

   reduced_ratio_list.pushback(0);
   reduced_potvec_list.pushback(0.0);

   for(int t=0;t!=this->m_n_times;t++)
   {
      test=1;
      potential_vector(t,vector);
      mod=(vector[2]/fabs(vector[2]))*sqrt(pow(vector[0],2)+pow(vector[1],2)+pow(vector[2],2));
      ratio=(mod/fabs(mod))*floor(fabs(mod)/H->pot_vec_thresh());
      ratio_list[t]=ratio;
      for(int i=0;i!=reduced_ratio_list.size();i++)
      {
         if(ratio_list[t]==reduced_ratio_list[i])
            test*=0;
      }
      if(test)
      {
         reduced_ratio_list.pushback(ratio_list[t]);
         reduced_potvec_list_x.pushback(vector[0]);
         reduced_potvec_list_y.pushback(vector[1]);
         reduced_potvec_list_z.pushback(vector[2]);
         this->pice_time_mapping[t]=reduced_ratio_list.size()-1;
      }
      else
      {
         for(int i=0;i!=reduced_ratio_list.size();i++)
         {
            if(ratio_list[t]==reduced_ratio_list[i])
            {
               this->pice_time_mapping[t]=i;
            }
         }
      }
   }

   this->mapping_size=reduced_ratio_list.size();
   this->pot_vec_reduced_mapping=new double*[this->mapping_size];

   for(int i=0;i!=this->mapping_size;i++)
   {
      this->pot_vec_reduced_mapping[i]=new double[3];
      this->pot_vec_reduced_mapping[i][0]=reduced_potvec_list_x[i];
      this->pot_vec_reduced_mapping[i][1]=reduced_potvec_list_y[i];
      this->pot_vec_reduced_mapping[i][2]=reduced_potvec_list_z[i];
   }

   this->m_PICE_sto_x=new std::complex<double> **[reduced_ratio_list.size()];
   this->m_PICE_sto_y=new std::complex<double> **[reduced_ratio_list.size()];
   this->m_PICE_sto_z=new std::complex<double> **[reduced_ratio_list.size()];

   for(int i=0;i!=reduced_ratio_list.size();i++)
   {
      this->m_PICE_sto_x[i]=new *std::complex<double>[n_states_neut*n_states_cat*this->m_n_states_cont];
      this->m_PICE_sto_y[i]=new *std::complex<double>[n_states_neut*n_states_cat*this->m_n_states_cont];
      this->m_PICE_sto_z[i]=new *std::complex<double>[n_states_neut*n_states_cat*this->m_n_states_cont];

      for(int t=0;t!=n_states_neut*n_states_cat*this->m_n_states_cont;t++)
      {
         this->m_PICE_sto_x[i][t]=new std::complex<double>[tgsize_x];
         this->m_PICE_sto_y[i][t]=new std::complex<double>[tgsize_x];
         this->m_PICE_sto_z[i][t]=new std::complex<double>[tgsize_x];
      }
   }
   delete [] ratio_list;
}
//##########################################################################
//
//END OF HAMILTON MATRIX OBJECT INITIALIZATION AND SETTINGS
