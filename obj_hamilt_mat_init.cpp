

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


//HAMILTON MATRIX OBJECT INITIALIZATION AND SETTINGS
//##########################################################################
//
//##########################################################################
hamilton_matrix::hamilton_matrix(int gsize_x,int tgsize_x,int small_gsize_x,int n_states_neut,int n_states_cat,int n_k,int n_angles,double kmin,double kmax,double xmin,double xmax,double mass,int n_times,double h,double efield_thresh,double pot_vec_thresh,std::string pice_data_loc)
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
   double delta_x((xmax-xmin)/gsize_x);
   std::cout<<"Size of the pixel is "<<delta_x<<std::endl;

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
/*   //MANUAL PICE SIGMN CORRECTION PARAMETRING
   ifstream sign_corr_str;
   stringstream sign_corr_ss;
   string sign_corr_s;
   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      sign_corr_ss.str("");
      sign_corr_ss<<"/data1/home/stephan/LiH_512_points_pice_16_05_18/sign_pice_corr_"<<i<<"_0.txt";
      sign_corr_s=sign_corr_ss.str();
      sign_corr_str.open(sign_corr_s.c_str());
      if(!sign_corr_str.is_open())
      {
         std::cout<<"SIGN CORRECTION FILE NOT FOUND. IGNORING"<<std::endl;
      }
      else
      {
         for(int j=0;j!=this->m_small_gsize_x;j++)
         {
            sign_corr_str>>this->sign_corr[i][j];
         }
      }
   }
   //END OF MANUAL PICE SIGMN CORRECTION PARAMETRING
*/
   //initialize potential energy surfaces arrays
   std::cout<<"initializing PES arrays...";
   this->m_pot_neut=new double [tgsize_x*n_states_neut];
   std::cout<<" ...";
   this->m_pot_cat=new double [tgsize_x*n_states_cat*this->m_n_states_cont];
   std::cout<<"PES arrays initialized!"<<std::endl;
   //initialize dipole moment surfaces arrays
   //of the neutral
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
   //of the cation
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
   //initialize photoionization coupling elements surfaces arrays
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
   for(int i=0;i!=this->m_n_k;i++)
   {
      this->k_modulus[i]=kmin+i*(kmax-kmin)/this->m_n_k;
      for(int j=0;j!=this->m_n_angles;j++)
      {
          m_dk_vec[i*this->m_n_angles+j]=2*4*acos(-1)*this->k_modulus[i]*(kmax-kmin)/(pow(2*acos(-1),3)*this->m_n_k*this->m_n_angles);
      }
   }
   std::cout<<"momentum vectors arrays initialized!"<<std::endl;
   //initialize Non-adiabatic coupling surfaces arrays
   std::cout<<"initializing NAC arrays...";
   this->m_NAC=new double*[n_states_neut*n_states_neut];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_neut;i++)
   {
      this->m_NAC[i]=new double [tgsize_x];
   }
   std::cout<<" ...";
   std::cout<<" NAC initialized!"<<std::endl;
   //Initialize and set up kinetic energy matrix
   std::cout<<"initializing kinetic energy array...";
   this->kinetic_energy=new double[tgsize_x*tgsize_x];
   this->derivative_matrix=new double[tgsize_x*tgsize_x];
   //FINITE DIFFERENCE METHOD IMPLEMENTATION OF KINETIC ENERGY (ABRAMOWITZ EQ. 25.3.24)
   for(int i=0;i!=tgsize_x;i++)
   {
      for(int j=0;j!=tgsize_x;j++)
      {
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
//##########################################################################
void hamilton_matrix::set_pot_neut(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   double var(0);
   int dgsize (this->m_tgsize_x-this->m_gsize_x);

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
         for(int j=dgsize;j!=this->m_tgsize_x;j++)
         {
            input_file>>this->m_pot_neut[this->m_tgsize_x*i+j];
           // cout<<this->m_pot_neut[this->m_gsize_x*i+j]<<endl;
         }
         input_file.close();
         for(int j=1;j<=dgsize;j++)
         {
            this->m_pot_neut[this->m_tgsize_x*i+(dgsize-j)]=this->m_pot_neut[this->m_tgsize_x*i+dgsize]+j*(this->m_pot_neut[this->m_tgsize_x*i+dgsize]-this->m_pot_neut[this->m_tgsize_x*i+dgsize+1])+5e-3*j*j;
         }
      }
      else
      {
         cout<<"ERROR POTENTIAL FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
         exit;
      }
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_pot_cat(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   int dgsize (this->m_tgsize_x-this->m_gsize_x);

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
         for(int j=dgsize;j!=this->m_tgsize_x;j++)
         {
            input_file>>this->m_pot_cat[this->m_tgsize_x*i+j];
            //cout<<this->m_pot_cat[this->m_gsize_x*i+j]<<std::endl;
         }
         input_file.close();
         for(int j=1;j<=dgsize;j++)
         {
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
   int dgsize(this->m_tgsize_x-this->m_gsize_x);

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
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    input_file>>temp;
                    this->m_dmx_neut[i*this->m_n_states_neut+j][k]=temp;
                    this->m_dmx_neut[j*this->m_n_states_neut+i][k] = this->m_dmx_neut[i*this->m_n_states_neut+j][k];
                }
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmx_neut[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmx_neut[i*this->m_n_states_neut+j][dgsize];
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
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    input_file>>temp;
                    this->m_dmy_neut[i*this->m_n_states_neut+j][k]=temp;
                    this->m_dmy_neut[j*this->m_n_states_neut+i][k]=this->m_dmy_neut[i*this->m_n_states_neut+j][k];
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmy_neut[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmy_neut[i*this->m_n_states_neut+j][dgsize];
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
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    input_file>>temp;
                    this->m_dmz_neut[i*this->m_n_states_neut+j][k]=temp;
                    this->m_dmz_neut[j*this->m_n_states_neut+i][k]=this->m_dmz_neut[i*this->m_n_states_neut+j][k];
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmz_neut[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmz_neut[i*this->m_n_states_neut+j][dgsize];
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
   int dgsize(this->m_tgsize_x-this->m_gsize_x);

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
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    input_file>>this->m_dmx_cat[i*this->m_n_states_cat+j][k];
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmx_cat[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmx_cat[i*this->m_n_states_neut+j][dgsize];
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
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    input_file>>this->m_dmy_cat[i*this->m_n_states_cat+j][k];
                }
                input_file.close();
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmy_cat[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmy_cat[i*this->m_n_states_neut+j][dgsize];
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
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                    input_file>>this->m_dmz_cat[i*this->m_n_states_cat+j][k];
                }
                for(int k=1;k<=dgsize;k++)
                {
                    this->m_dmz_cat[i*this->m_n_states_neut+j][dgsize-k]=this->m_dmz_cat[i*this->m_n_states_neut+j][dgsize];
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

void hamilton_matrix::set_PICE(double* pot_vec)
{
   const int nstneut=this->m_n_states_neut;
   const int nstcat=this->m_n_states_cat;
   const int nk=this->m_n_k;
   const int ndist=this->m_n_angles;
   //std::cout<<"probe"<<std::endl;

   std::complex<double> pice_x;
   std::complex<double> pice_y;
   std::complex<double> pice_z;

   using namespace std;
   double duration;
   clock_t begin;
   clock_t end;

   int ratio(this->m_gsize_x/this->m_small_gsize_x);

   int k(0);
   int l(0);
   for(int x=0;x!=this->m_gsize_x;x++)
   {
      std::cout<<x<<" position "<<std::endl;
      for( int i=0;i!=nstneut;i++)
      {
         std::cout<<"stneut"<<i<<std::endl;
         for(int j=0;j!=nstcat;j++)
         {
         std::cout<<"stcat"<<j<<std::endl;
         begin = clock();
            #pragma omp parallel for private(k,l,pice_x,pice_y,pice_z) shared(x,i,j,nk,ndist,ratio)
            for(k=0;k<nk;k++)
            {
               for(l=0;l<ndist;l++)
               {
                   if(x%ratio == 0 || x==0)
                   {
                       this->pice_data->fill_pice(&pice_x,&pice_y,&pice_z,x,i,j,this->k_orientation[0][l],this->k_orientation[1][l],this->k_modulus[k],pot_vec);
                       this->m_PICE_x[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]=pice_x;
                       this->m_PICE_y[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]=pice_y;
                       this->m_PICE_z[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]=pice_z;
                   }
                   else
                   {
                       this->m_PICE_x[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]=this->m_PICE_x[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x-1];
                       this->m_PICE_y[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]=this->m_PICE_y[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x-1];
                       this->m_PICE_z[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]=this->m_PICE_z[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x-1];
                   }
//                   std::cout<<this->m_PICE_x[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]<<","<<this->m_PICE_y[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]<<","<<this->m_PICE_z[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]<<std::endl;
//                   exit(EXIT_SUCCESS);
                   //std::cout<<"probe "<<x<<","<<i<<","<<j<<","<<k<<","<<l<<std::endl;
               }
            }
            end = clock();
           duration=double(end-begin)/CLOCKS_PER_SEC;
           std::cout<<duration/16<<"s"<<std::endl;
         }
      }

//      int i=18;
//      int j=0;
//      int k=2;
//      int l=1;
//      std::cout<<this->m_PICE_x[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]<<","<<this->m_PICE_y[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]<<","<<this->m_PICE_z[i*nstcat*nk*ndist+j*nk*ndist+k*ndist+l][x]<<std::endl;
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_NAC(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   int dgsize(this->m_tgsize_x-this->m_gsize_x);

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
                for (int k=dgsize; k!=this->m_tgsize_x; k++)
                {
                   input_file>>this->m_NAC[i*this->m_n_states_neut+j][k];
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
   }
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
   }
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
/*
//##########################################################################
//
//END OF HAMILTON MATRIX OBJECT INITIALIZATION AND SETTINGS
//##########################################################################
std::complex<double> hamilton_matrix::show_pice()
{
   return this->m_PICE_z[i*this->m_n_states_cat*this->m_n_states_cont+j*this->m_n_states_cont+k*ndist+l][x];
}
void hamilton_matrix::spherical_extract_from_cube(int neut_state,int cat_state,int r_index,int component,double *Recube,double *Imcube,double *pot_vec)
{
   double xsphere(0);
   double ysphere(0);
   double zsphere(0);
   int x_index(0);
   int y_index(0);
   int z_index(0);
   double k(0);
   double deltak(0);
   int n_points_sphere(this->m_n_angles);
   int box_length_x(2*acos(-1)/((this->m_kxmax-m_kxmin)/this->m_n_kx));
   int box_length_y(2*acos(-1)/((this->m_kymax-m_kymin)/this->m_n_ky));
   int box_length_z(2*acos(-1)/((this->m_kzmax-m_kzmin)/this->m_n_kz));
   int dgsize=this->m_tgsize_x-this->m_gsize_x;

//   std::cout<<"box lengths parameters: "<<box_length_x<<"  "<<box_length_y<<"   "<<box_length_z<<std::endl;

  // std::cout<<"Extracting spherical distribution PICE from cube file...";

   if(pot_vec==NULL)
   {
      for(int j=0;j!=this->m_n_k;j++)
      {
         k=this->k_modulus[j];
         for(int i=0;i!=n_points_sphere;i++)
         {
            xsphere=k*sin(this->k_orientation[0][i])*cos(this->k_orientation[1][i]);
            ysphere=k*sin(this->k_orientation[0][i])*sin(this->k_orientation[1][i]);
            zsphere=k*cos(this->k_orientation[0][i]);
 //           std::cout<<"probe sphere"<<xsphere<<"   "<<ysphere<<"   "<<zsphere<<"   "<<std::endl;

            x_index=int(round((xsphere-this->m_kxmin)*this->m_n_kx/(this->m_kxmax-this->m_kxmin)));
            y_index=int(round((ysphere-this->m_kymin)*this->m_n_ky/(this->m_kymax-this->m_kymin)));
            z_index=int(round((zsphere-this->m_kzmin)*this->m_n_kz/(this->m_kzmax-this->m_kzmin)));
 //           std::cout<<"probe cube"<<x_index<<"   "<<y_index<<"   "<<z_index<<"   "<<std::endl;
 
            / *if(j==19)
            {
               std::cout<<x_index<<","<<y_index<<","<<z_index<<std::endl;
            } * /

            if(component==0)
                this->m_PICE_x[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+j*n_points_sphere+i][r_index]=std::complex<double>(Recube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index],Imcube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index])*this->sign_corr[neut_state][r_index-dgsize];
            else if(component==1)
                this->m_PICE_y[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+j*n_points_sphere+i][r_index]=std::complex<double>(Recube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index],Imcube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index])*this->sign_corr[neut_state][r_index-dgsize];
            else if(component=2)
                this->m_PICE_z[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+j*n_points_sphere+i][r_index]=std::complex<double>(Recube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index],Imcube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index])*this->sign_corr[neut_state][r_index-dgsize];
            else
            {
               std::cout<<"ERROR ASKED FOR AN OUT OF RANGE COMPONENT OF A 3D VECTOR"<<std::endl; 
               exit(EXIT_FAILURE);
            }
         }
      }
      deltak=(this->k_modulus[this->m_n_k-1]-this->k_modulus[0])/this->m_n_k;
      for(int i=0;i!=this->m_n_states_cont;i++)
      {
         this->m_dk_vec[i]=deltak*(box_length_x*box_length_y*box_length_z)*this->k_modulus[(i-i%this->m_n_angles)/this->m_n_angles]*this->m_n_k/(pow(acos(-1),2)*2*this->m_n_states_cont);
      }
   }
      else
      {
      for(int j=0;j!=this->m_n_k;j++)
      {
         k=this->k_modulus[j];
         for(int i=0;i!=n_points_sphere;i++)
         {
            xsphere=k*sin(this->k_orientation[0][i])*cos(this->k_orientation[1][i])-pot_vec[0];
            ysphere=k*sin(this->k_orientation[0][i])*sin(this->k_orientation[1][i])-pot_vec[1];
            zsphere=k*cos(this->k_orientation[0][i])-pot_vec[2];

            x_index=int(round((xsphere-this->m_kxmin)*this->m_n_kx/(this->m_kxmax-this->m_kxmin)));
            y_index=int(round((ysphere-this->m_kymin)*this->m_n_ky/(this->m_kymax-this->m_kymin)));
            z_index=int(round((zsphere-this->m_kzmin)*this->m_n_kz/(this->m_kzmax-this->m_kzmin)));

            if(component==0)
                this->m_PICE_x[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+j*n_points_sphere+i][r_index]=std::complex<double>(Recube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index],Imcube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index])*this->sign_corr[neut_state][r_index-dgsize];
            else if(component==1)
                this->m_PICE_y[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+j*n_points_sphere+i][r_index]=std::complex<double>(Recube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index],Imcube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index])*this->sign_corr[neut_state][r_index-dgsize];
            else if(component=2)
                this->m_PICE_z[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+j*n_points_sphere+i][r_index]=std::complex<double>(Recube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index],Imcube[x_index*this->m_n_ky*this->m_n_kz+y_index*this->m_n_kz+z_index])*this->sign_corr[neut_state][r_index-dgsize];
            else
            {
               std::cout<<"ERROR ASKED FOR AN OUT OF RANGE COMPONENT OF A 3D VECTOR"<<std::endl; 
               exit(EXIT_FAILURE);
            }
         }
      }
      deltak=(this->k_modulus[this->m_n_k-1]-this->k_modulus[0])/this->m_n_k;
      for(int i=0;i!=this->m_n_states_cont;i++)
      {
         this->m_dk_vec[i]=deltak*(box_length_x*box_length_y*box_length_z)*this->k_modulus[(i-i%this->m_n_angles)/this->m_n_angles]*this->m_n_k/(pow(acos(-1),2)*2*this->m_n_states_cont);
      }

    //  std::cout<<"Done!"<<std::endl;
      }
}
//##########################################################################
//
//##########################################################################
bool hamilton_matrix::cube_reader(std::string MO_cube_loc,double *cube_array,bool extract_dimensions)
{
   using namespace std;
   string temp;
   double dtemp;
   int num_of_nucl(0);
   ifstream MO_cube_out;
   MO_cube_out.open(MO_cube_loc.c_str());
         if (!MO_cube_out.is_open())
         {
             cout<<"ERROR: CANNOT OPEN CUBE FILE "<<MO_cube_loc.c_str()<<endl;
             return 0;
         }

//         std::cout<<"Reading cube file "<<MO_cube_loc.c_str()<<std::endl;;
   if(extract_dimensions)
   {
         getline(MO_cube_out,temp);
         getline(MO_cube_out,temp);
         MO_cube_out>>temp;
   //      std::cout<<temp<<std::endl;
         num_of_nucl=fabs(atoi(temp.c_str()));
         MO_cube_out>>this->m_kxmin;
         MO_cube_out>>this->m_kymin;
         MO_cube_out>>this->m_kzmin;
         MO_cube_out>>this->m_n_kx;
         MO_cube_out>>dtemp;
         std::cout<<dtemp<<std::endl;
         this->m_kxmax=this->m_kxmin+dtemp*this->m_n_kx;
         MO_cube_out>>dtemp;
         MO_cube_out>>dtemp;
         MO_cube_out>>this->m_n_ky;
         MO_cube_out>>dtemp;
         MO_cube_out>>dtemp;
         this->m_kymax=this->m_kymin+dtemp*this->m_n_ky;
         MO_cube_out>>dtemp;
         MO_cube_out>>this->m_n_kz;
         MO_cube_out>>dtemp;
         MO_cube_out>>dtemp;
         MO_cube_out>>dtemp;
         this->m_kzmax=this->m_kzmin+dtemp*this->m_n_kz;
         for(int i=0;i!=this->m_n_k;i++)
         {
             this->k_modulus[i]=(this->m_kzmax-0.1)*(i+1)/this->m_n_k;
         }

         
         std::cout<<"kxmin = "<<this->m_kxmin<<std::endl;
         std::cout<<"kymin = "<<this->m_kymin<<std::endl;
         std::cout<<"kzmin = "<<this->m_kzmin<<std::endl;
         std::cout<<"kxmax = "<<this->m_kxmax<<std::endl;
         std::cout<<"kymax = "<<this->m_kymax<<std::endl;
         std::cout<<"kzmax = "<<this->m_kzmax<<std::endl;
         std::cout<<"nkx = "<<this->m_n_kx<<std::endl;
         std::cout<<"nky = "<<this->m_n_ky<<std::endl;
         std::cout<<"nkz = "<<this->m_n_kz<<std::endl;
         
         //std::cout<<"read number of atoms: "<<num_of_nucl<<std::endl;

     //   std::cout<<"Done!"<<std::endl;
   }
   else
   {
         getline(MO_cube_out,temp);
         getline(MO_cube_out,temp);
         MO_cube_out>>temp;
         num_of_nucl=fabs(atoi(temp.c_str()));
         //std::cout<<"read number of atoms: "<<num_of_nucl<<std::endl;

         for(int i=0;i!=15+fabs(num_of_nucl)*5+2;i++)
         {
            MO_cube_out>>temp;
            //std::cout<<temp<<std::endl;
         }

        for(int i=0;i!=this->m_n_kx*this->m_n_ky*this->m_n_kz;i++)
        {
            MO_cube_out>>cube_array[i];
        }
       // std::cout<<"Done!"<<std::endl;
   }

   MO_cube_out.close();
   return 1;
}
//##########################################################################
//
//##########################################################################

//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_PICE(std::string file_address,double* pot_vec)
{
   using namespace std;
   stringstream name_indenter2;
   stringstream name_indenter;
   double *position=new double[this->m_small_gsize_x];
   bool test(0);
   bool test2(0);
   double temp2(0);
   int index_pos=0;
   string filename;
   int dgsize(this->m_tgsize_x-this->m_gsize_x);
   double temp;
   double pos;
   double r_val(0);
   int x(0);
   
   



      ifstream input_file;
      name_indenter2.str("");
      name_indenter2<<file_address<<"coordinates.input";
      filename=name_indenter2.str();
      input_file.open(filename.c_str());
      if(!input_file.is_open())
      {
         std::cout<<"POSITION INPUT SCRIPT FILE CANNOT BE FOUND "<<std::endl;
         exit(EXIT_FAILURE);
      }
      for(int i=0;i!=this->m_small_gsize_x;i++)
      {
         input_file>>position[i];
      }
      input_file.close();

      pos=position[0];
      name_indenter2.str("");
      name_indenter2<<file_address<<"RePICE_"<<pos<<"_X_"<<0<<"_"<<0<<".txt";
      filename=name_indenter2.str();
      if(!this->cube_reader(filename,NULL,1))
      {
         cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
         exit(EXIT_FAILURE);
      }

      std::cout<<"momentum cube parameters read !"<<std::endl;

   double *Retemp_cube_x=new double [this->m_n_kx*this->m_n_ky*this->m_n_kz];
   double *Imtemp_cube_x=new double [this->m_n_kx*this->m_n_ky*this->m_n_kz];
   double *Retemp_cube_y=new double [this->m_n_kx*this->m_n_ky*this->m_n_kz];
   double *Imtemp_cube_y=new double [this->m_n_kx*this->m_n_ky*this->m_n_kz];
   double *Retemp_cube_z=new double [this->m_n_kx*this->m_n_ky*this->m_n_kz];
   double *Imtemp_cube_z=new double [this->m_n_kx*this->m_n_ky*this->m_n_kz];


   if(file_address != "")//If we give the files address, it means that it is the first time the routine is called => We generate PICE for zero external field and we generate a random angular distribution.
   {
      this->m_PICE_address=file_address;

      for(int i=0;i!=this->m_n_states_neut;i++)
      {
         for(int j=0;j!=this->m_n_states_cat;j++)
         {
            std::cout<<"Getting PICE for states "<<i<<" - "<<j<<std::endl;

            for(int k=dgsize;k!=this->m_tgsize_x;k++)
            {
               test2=0;

               r_val=0.529*(this->m_xmin+(k-dgsize)*(this->m_xmax-this->m_xmin)/this->m_gsize_x);

               std::cout<<"Loop "<<k-dgsize<<", position "<<r_val<<std::endl;
/ *
               for(int l=0;l!=this->m_small_gsize_x;l++)
               {
//                  std::cout<<" now comparing : "<<position[l]<<" <= "<<r_val<<" < "<<position[l+1]<<"..."<<std::endl;
                  if( r_val < position[l+1] && r_val >= position[l] )
                  {
//                     std::cout<<"TRUE !"<<std::endl;
                     x=l;
                     if(l==m_small_gsize_x-1 || l==0)
                     {
//                        std::cout<<" Warning! out of bound for interpolated data! "<<std::endl;
                     }
                     else if(temp2 == position[l] && l!=0)
                     {
//                        std::cout<<"Checking test2 ! "<<temp2<<" ; "<<position[l]<<std::endl;
                        test2=1;
                     }
//                     break;
                  }
               }
               if(test2)
               {
//                  std::cout<<"test 2 checked "<<r_val<<std::endl;
                  this->spherical_extract_from_cube(i,j,k,0,Retemp_cube_x,Imtemp_cube_x,NULL);
                  this->spherical_extract_from_cube(i,j,k,1,Retemp_cube_y,Imtemp_cube_y,NULL);
                  this->spherical_extract_from_cube(i,j,k,2,Retemp_cube_z,Imtemp_cube_z,NULL);
                  continue;
               }
* /
               x=k-dgsize;
               temp2=position[x];
              // std::cout<<"Loading new PICE file for R = "<<position[x]<<". New ref value = "<<temp2<<std::endl;

               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"RePICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Retemp_cube_x))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               std::cout<<x<<"/"<<this->m_tgsize_x<<"..."<<std::endl;
               pos=position[x-dgsize];
               name_indenter.str("");
               name_indenter<<file_address<<"RePICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Retemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"RePICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Retemp_cube,1))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"ImPICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Imtemp_cube_x))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"ImPICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Imtemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"ImPICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Imtemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }* /
               this->spherical_extract_from_cube(i,j,k,0,Retemp_cube_x,Imtemp_cube_x,NULL);

               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"RePICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Retemp_cube_y))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"RePICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Retemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"RePICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Retemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"ImPICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Imtemp_cube_y))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"ImPICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Imtemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"ImPICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Imtemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               this->spherical_extract_from_cube(i,j,k,1,Retemp_cube_y,Imtemp_cube_y,NULL);

               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"RePICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Retemp_cube_z))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"RePICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Retemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"RePICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Retemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"ImPICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Imtemp_cube_z))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"ImPICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Imtemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"ImPICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Imtemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /

               //std::cout<<"pice "<<i<<"   "<<j<<"   "<<k<<"   "<<Retemp_cube_z[49*this->m_n_ky*this->m_n_kz+51*this->m_n_kz+49]<<"   "<<Imtemp_cube_z[49*this->m_n_ky*this->m_n_kz+51*this->m_n_kz+49]<<std::endl;
               this->spherical_extract_from_cube(i,j,k,2,Retemp_cube_z,Imtemp_cube_z,NULL);

            }
         }
      }
   }
   else
   {

      for(int i=0;i!=this->m_n_states_neut;i++)
      {
         for(int j=0;j!=this->m_n_states_cat;j++)
         {
            std::cout<<"Getting PICE for states "<<i<<" - "<<j<<std::endl;
            std::cout<<"point...";
            for(int k=dgsize;k!=this->m_tgsize_x;k++)
            {

/ *               test2=0;

               r_val=0.529*(this->m_xmin+(k-dgsize)*(this->m_xmax-this->m_xmin)/this->m_gsize_x);

               std::cout<<"Loop "<<k<<", position "<<r_val<<std::endl;

               for(int l=0;l!=this->m_small_gsize_x;l++)
               {
//                  std::cout<<" now comparing : "<<position[l]<<" <= "<<r_val<<" < "<<position[l+1]<<"..."<<std::endl;
                  if( r_val < position[l+1] && r_val >= position[l] )
                  {
//                     std::cout<<"TRUE !"<<std::endl;
                     x=l;
                     if(l==m_small_gsize_x-1 || l==0)
                     {
//                        std::cout<<" Warning! out of bound for interpolated data! "<<std::endl;
                     }
                     else if(temp2 == position[l] && l!=0)
                     {
//                        std::cout<<"Checking test2 ! "<<temp2<<" ; "<<position[l]<<std::endl;
                        test2=1;
                     }
//                     break;
                  }
               }
               if(test2)
               {
//                  std::cout<<"test 2 checked "<<r_val<<std::endl;
                  this->spherical_extract_from_cube(i,j,k,0,Retemp_cube_x,Imtemp_cube_x,pot_vec);
                  this->spherical_extract_from_cube(i,j,k,1,Retemp_cube_y,Imtemp_cube_y,pot_vec);
                  this->spherical_extract_from_cube(i,j,k,2,Retemp_cube_z,Imtemp_cube_z,pot_vec);
                  continue;
               }
* /
               x=k-dgsize;
               temp2=position[x];

               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"RePICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Retemp_cube_x))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"RePICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Retemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"RePICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Retemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"ImPICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Imtemp_cube_x))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"ImPICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Imtemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"ImPICE_"<<pos<<"_X_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Imtemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               this->spherical_extract_from_cube(i,j,k,0,Retemp_cube_x,Imtemp_cube_x,pot_vec);

               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"RePICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Retemp_cube_y))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"RePICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Retemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"RePICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Retemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"ImPICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Imtemp_cube_y))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"ImPICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Imtemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"ImPICE_"<<pos<<"_Y_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Imtemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               this->spherical_extract_from_cube(i,j,k,1,Retemp_cube_y,Imtemp_cube_y,pot_vec);

               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"RePICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Retemp_cube_z))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"RePICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Retemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"RePICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Retemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               test=0;
               index_pos=0;
               do
               {
                  pos=position[x-index_pos];
                  name_indenter.str("");
                  name_indenter<<file_address<<"ImPICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();

                  if(!this->cube_reader(filename,Imtemp_cube_z))
                  {
                     if(x==0)
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                     else
                     {
                        index_pos++;
                        cout<<" PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"GO TO R="<<position[x-index_pos]<<endl;
                        continue;
                     }
                  }
                  else
                     test=1;
               }while(test!=1);
               / *
               name_indenter.str("");
               name_indenter<<file_address<<"ImPICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
               filename=name_indenter.str();

               if(!this->cube_reader(filename,Imtemp_cube))
               {
                  if(x==0)
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
                  else
                  {
                     pos=position[x-1-dgsize];
                     name_indenter.str("");
                     name_indenter<<file_address<<"ImPICE_"<<pos<<"_Z_"<<i<<"_"<<j<<".txt";
                     filename=name_indenter.str();
                     if(!this->cube_reader(filename,Imtemp_cube))
                     {
                        cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                        exit(EXIT_FAILURE);
                     }
                  }
               }
               * /
               this->spherical_extract_from_cube(i,j,k,2,Retemp_cube_z,Imtemp_cube_z,pot_vec);
            }
         }
      }
   }
   delete [] Retemp_cube_x;
   delete [] Imtemp_cube_x;
   delete [] Retemp_cube_y;
   delete [] Imtemp_cube_y;
   delete [] Retemp_cube_z;
   delete [] Imtemp_cube_z;
/ *   
   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      for(int j=0;j!=this->m_n_states_cat;j++)
      {
         std::cout<<"Getting PICE for states "<<i<<" - "<<j<<std::endl;
         for(int x=dgsize;x!=this->m_tgsize_x;x++)
         {
            pos=position[x-dgsize];
            name_indenter.str("");
            name_indenter<<file_address<<pos<<"_"<<i<<"_"<<j<<".txt";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(!input_file.is_open())
            {
               if(x==0)
               {
                  cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                  exit(EXIT_FAILURE);
               }
               else
               {
                  pos=position[x-1-dgsize];
                  name_indenter.str("");
                  name_indenter<<file_address<<pos<<"_"<<i<<"_"<<j<<".txt";
                  filename=name_indenter.str();
                  input_file.open(filename.c_str());
                  if(!input_file.is_open())
                  {
                     cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
                     exit(EXIT_FAILURE);
                  }
               }
            }
            {
               for(int k=0;k!=this->m_n_k;k++)
               {
                  for(int l=0;l!=this->m_n_angles;l++)
                  {
                     input_file>>this->k_modulus[k];
                     if(k==50)
                     {
                        input_file>>this->k_orientation[0][l];
                        input_file>>this->k_orientation[1][l];
                     }
                     else
                     {
                        input_file>>temp;
                        input_file>>temp;
                     }
                     input_file>>Re_value;
                     input_file>>Im_value;
                     this->m_PICE_x[i*this->m_n_states_cat*this->m_n_states_cont+j*this->m_n_states_cont+k*this->m_n_angles+l][x]=std::complex<double>(Re_value,Im_value);
                     input_file>>Re_value;
                     input_file>>Im_value;
                     this->m_PICE_y[i*this->m_n_states_cat*this->m_n_states_cont+j*this->m_n_states_cont+k*this->m_n_angles+l][x]=std::complex<double>(Re_value,Im_value);
                     input_file>>Re_value;
                     input_file>>Im_value;
                     this->m_PICE_z[i*this->m_n_states_cat*this->m_n_states_cont+j*this->m_n_states_cont+k*this->m_n_angles+l][x]=std::complex<double>(Re_value,Im_value);

                  }
               }
               input_file.close();
            }
         }
      }
   }
   this->m_dk_vec=new double[this->m_n_states_cont];
   double deltak=(this->k_modulus[this->m_n_k-1]-this->k_modulus[0])/this->m_n_k;
   for(int i=0;i!=this->m_n_states_cont;i++)
   {
      this->m_dk_vec[i]=deltak*(55*55*55)*this->k_modulus[(i-i%this->m_n_angles)/this->m_n_angles]*this->m_n_k/(pow(acos(-1),2)*2*this->m_n_states_cont);
   }

   std::cout<<"Got all PICE ! "<<std::endl<<"NOT Determining closest pair of plane waves in reciprocal space"<<std::endl;
   for(int t=0;t!=this->m_n_times;t++)
   {
      for(int i=0;i!=this->m_n_states_cont;i++)
         this->translation_vector[i][t]=i;
   }
   
      //DETERMINING THE CLOSEST PAIR OF PLANE WAVES IN THE RECIPROCAL SPACE

                               for(int m=0;m!=this->m_n_states_cont;m++)
                               {   
                                  std::cout<<"element "<<m<<"/"<<this->m_n_states_cont<<"checked"<<std::endl;
                                  for(int v=m+1;v!=this->m_n_states_cont;v++)
                                  {   
                                     if(m==0 && v==1)
                                        this->min_distance=sqrt(pow(this->k_modulus[m%this->m_n_angles]*sin(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])*cos(this->k_orientation[1][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*sin(this->k_orientation[0][v-v%this->m_n_k])*cos(this->k_orientation[1][v-v%this->m_n_k]),2)+pow(this->k_modulus[m%this->m_n_angles]*sin(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])*sin(this->k_orientation[1][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*sin(this->k_orientation[0][v-v%this->m_n_k])*sin(this->k_orientation[1][v-v%this->m_n_k]),2)+pow(this->k_modulus[m%this->m_n_angles]*cos(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*cos(this->k_orientation[0][v-v%this->m_n_k]),2));
                                     if(sqrt(pow(this->k_modulus[m%this->m_n_angles]*sin(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])*cos(this->k_orientation[1][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*sin(this->k_orientation[0][v-v%this->m_n_k])*cos(this->k_orientation[1][v-v%this->m_n_k]),2)+pow(this->k_modulus[m%this->m_n_angles]*sin(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])*sin(this->k_orientation[1][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*sin(this->k_orientation[0][v-v%this->m_n_k])*sin(this->k_orientation[1][v-v%this->m_n_k]),2)+pow(this->k_modulus[m%this->m_n_angles]*cos(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*cos(this->k_orientation[0][v-v%this->m_n_k]),2))<min_distance)
                                     {   
                                        this->min_distance=sqrt(pow(this->k_modulus[m%this->m_n_angles]*sin(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])*cos(this->k_orientation[1][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*sin(this->k_orientation[0][v-v%this->m_n_k])*cos(this->k_orientation[1][v-v%this->m_n_k]),2)+pow(k_modulus[m%this->m_n_angles]*sin(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])*sin(this->k_orientation[1][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*sin(this->k_orientation[0][v-v%this->m_n_k])*sin(this->k_orientation[1][v-v%this->m_n_k]),2)+pow(this->k_modulus[m%this->m_n_angles]*cos(this->k_orientation[0][m-m%this->m_n_angles*this->m_n_angles])-this->k_modulus[v%this->m_n_k]*cos(this->k_orientation[0][v-v%this->m_n_k]),2));
                                     }   
                                  }   
                               }   
                               std::cout<<"Got closest pair !"<<std::endl;
   //END OF CLOSEST PAIR DETERMINATION
         for(int i=0;i!=this->m_n_states_cont;i++)
         {
            ref_px[i]=this->k_modulus[i%this->m_n_angles]*sin(this->k_orientation[0][i-i%this->m_n_angles*this->m_n_angles])*cos(this->k_orientation[1][i-i%this->m_n_angles*this->m_n_angles]);
            ref_py[i]=this->k_modulus[i%this->m_n_angles]*sin(this->k_orientation[0][i-i%this->m_n_angles*this->m_n_angles])*sin(this->k_orientation[1][i-i%this->m_n_angles*this->m_n_angles]);
            ref_pz[i]=this->k_modulus[i%this->m_n_angles]*cos(this->k_orientation[0][i-i%this->m_n_angles*this->m_n_angles]);
         }
//#pragma omp parallel for
   for(int t=0;t!=this->m_n_times;t++)
   {
         potential_vector(t,pot_vec);
         std::cout<<"Translating plane waves for time "<<t<<" / "<<this->m_n_times<<std::endl;
//#pragma omp parallel for
      if(sqrt(pow(pot_vec[0],2)+pow(pot_vec[1],2)+pow(pot_vec[2],2))>=min_distance && if (sqrt(pow(pot_vec[0],2)+pow(pot_vec[1],2)+pow(pot_vec[2],2)) >= 1e-5))
      {
//#pragma omp parallel for
         for(int i=0;i!=this->m_n_states_cont;i++)
         {
            new_px[i]=ref_px[i]-pot_vec[0];
            new_py[i]=ref_py[i]-pot_vec[1];
            new_pz[i]=ref_pz[i]-pot_vec[2];
            trial_px[i]=ref_px[i];
            trial_py[i]=ref_py[i];
            trial_pz[i]=ref_pz[i];

//#pragma omp parallel for
            for(int j=0;j!=this->m_n_states_cont;j++)
            {
               if(pow(new_px[i]-ref_px[i],2)+pow(new_py[i]-ref_py[i],2)+pow(new_pz[i]-ref_pz[i],2) < pow(new_px[i]-trial_px[i],2)+pow(new_py[i]-trial_py[i],2)+pow(new_pz[i]-trial_pz[i],2) && i!=j)
               {
                  this->translation_vector[i][t]=j;
                  std::cout<<t<<"   "<<i<<"   "<<this->translation_vector[i][t]<<std::endl;
                  trial_px[i]=ref_px[j];
                  trial_py[i]=ref_py[j];
                  trial_pz[i]=ref_pz[j];
               }
            }
         }
      }
   }
   * /
   

}
//##########################################################################
//
//##########################################################################
 * */
