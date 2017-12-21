
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
//HAMILTON MATRIX OBJECT COMPUTATION


//HAMILTON MATRIX OBJECT INITIALIZATION AND SETTINGS
//##########################################################################
//
//##########################################################################
hamilton_matrix::hamilton_matrix(int gsize_x,int n_states_neut,int n_states_cat,int n_k,int n_angles,double xmin,double xmax,double mass,int n_times,double h,double efield_thresh)
{
   //initialize grid parameters and time settings
   this->m_gsize_x=gsize_x;
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
   double delta_x((xmax-xmin)/gsize_x);
   std::cout<<"Size of the pixel is "<<delta_x<<std::endl;

   //initialize potenital energy surfaces arrays
   std::cout<<"initializing PES arrays...";
   this->m_pot_neut=new double [gsize_x*n_states_neut];
   std::cout<<" ...";
   this->m_pot_cat=new double [gsize_x*n_states_cat*this->m_n_states_cont];
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
      this->m_dmx_neut[i]=new double[gsize_x];
      this->m_dmy_neut[i]=new double[gsize_x];
      this->m_dmz_neut[i]=new double[gsize_x];
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
      this->m_dmx_cat[i]=new double[gsize_x];
      this->m_dmy_cat[i]=new double[gsize_x];
      this->m_dmz_cat[i]=new double[gsize_x];
   }
   std::cout<<"dipole arrays of the cation initialized!"<<std::endl;
   //initialize photoionization coupling elements surfaces arrays
   std::cout<<"initializing PICE arrays...";
   this->m_PICE_x=new std::complex<double> *[n_states_neut*n_states_cat*this->m_n_states_cont];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_cat*this->m_n_states_cont;i++)
   {
      this->m_PICE_x[i]=new std::complex<double>[gsize_x];
   }
   std::cout<<" ...";
   this->m_PICE_y=new std::complex<double> *[n_states_neut*n_states_cat*this->m_n_states_cont];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_cat*this->m_n_states_cont;i++)
   {
      this->m_PICE_y[i]=new std::complex<double>[gsize_x];
   }
   std::cout<<" ...";
   this->m_PICE_z=new std::complex<double> *[n_states_neut*n_states_cat*this->m_n_states_cont];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_cat*this->m_n_states_cont;i++)
   {
      this->m_PICE_z[i]=new std::complex<double>[gsize_x];
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
   std::cout<<" ...";
   this->translation_vector=new int*[this->m_n_states_cont];
   std::cout<<" ...";
   for(int i=0;i!=m_n_states_cont;i++)
   {
      this->translation_vector[i]=new int[m_n_times];
   }
   std::cout<<" ...";
   std::cout<<"momentum vectors arrays initialized!"<<std::endl;
   //initialize Non-adiabatic coupling surfaces arrays
   std::cout<<"initializing NAC arrays...";
   this->m_NAC=new double*[n_states_neut*n_states_neut];
   std::cout<<" ...";
   for(int i=0;i!=n_states_neut*n_states_neut;i++)
   {
      this->m_NAC[i]=new double [gsize_x];
   }
   std::cout<<" ...";
   std::cout<<" NAC initialized!"<<std::endl;
   //Initialize and set up kinetic energy matrix
   std::cout<<"initializing kinetic energy array...";
   this->kinetic_energy=new double[gsize_x*gsize_x];
   this->derivative_matrix=new double[gsize_x*gsize_x];
   //FINITE DIFFERENCE METHOD IMPLEMENTATION OF KINETIC ENERGY (ABRAMOWITZ EQ. 25.3.24)
   for(int i=0;i!=gsize_x;i++)
   {
      for(int j=0;j!=gsize_x;j++)
      {
         if(i==j+2)//(i,i-2)
         {
            this->kinetic_energy[i*gsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(-1);
         }
         else if(i==j+1)//(i,i-1)
         {
            this->kinetic_energy[i*gsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(16);
         }
         else if(i==j)//(i,i)
         {
            this->kinetic_energy[i*gsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(-30);
         }
         else if(i==j-1)//(i,i+1)
         {
            this->kinetic_energy[i*gsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(16);
         }
         else if(i==j-2)//(i,i+2)
         {
            this->kinetic_energy[i*gsize_x+j]=(-1/(2*mass))*(1/(12*delta_x*delta_x))*(-1);
         }
         else
            this->kinetic_energy[i*gsize_x+j]=0;
      }
   }
   std::cout<<"kinetic energy array initialized!"<<std::endl;
   //FINITE DIFFERENCE METHOD IMPLEMENTATION OF THE DERIVATIVE MATRIX
   for(int i=0;i!=gsize_x;i++)
   {
      for(int j=0;j!=gsize_x;j++)
      {
         if(i==j+2)//(i,i-2)
         {
            this->derivative_matrix[i*gsize_x+j]=(1/(12*delta_x))*(1);
         }
         else if(i==j+1)//(i,i-1)
         {
            this->derivative_matrix[i*gsize_x+j]=(1/(12*delta_x))*(-8);
         }
         else if(i==j-1)//(i,i+1)
         {
            this->derivative_matrix[i*gsize_x+j]=(1/(12*delta_x))*(8);
         }
         else if(i==j-2)//(i,i+2)
         {
            this->derivative_matrix[i*gsize_x+j]=(1/(12*delta_x))*(-1);
         }
         else
            this->derivative_matrix[i*gsize_x+j]=0;
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

   ifstream input_file;

   cout<<"gathering PES of the neutral"<<endl;

   for(int i=0;i!=m_n_states_neut;i++)
   {
      name_indenter.str("");
      name_indenter<<file_address<<i+1<<".input";
      filename=name_indenter.str();
      input_file.open(filename.c_str());
      
      cout<<"state "<<i+1<<"  ########################## in "<<filename.c_str()<<endl;
      if(input_file.is_open())
      {
         input_file.seekg(0);
         for(int j=0;j!=this->m_gsize_x;j++)
         {
            input_file>>this->m_pot_neut[this->m_gsize_x*i+j];
            //cout<<this->m_pot_neut[this->m_gsize_x*i+j]<<endl;
         }
         input_file.close();
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
         for(int j=0;j!=this->m_gsize_x;j++)
         {
            input_file>>this->m_pot_cat[this->m_gsize_x*i+j];
            //cout<<this->m_pot_cat[this->m_gsize_x*i+j]<<std::endl;
         }
         input_file.close();
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
void hamilton_matrix::set_dm_neut(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;

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
                for (int k=0; k!=this->m_gsize_x; k++)
                {
                    input_file>>this->m_dmx_neut[i*this->m_n_states_neut+j][k];
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
                for (int k=0; k!=this->m_gsize_x; k++)
                {
                    input_file>>this->m_dmy_neut[i*this->m_n_states_neut+j][k];
                }
                input_file.close();
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
                for (int k=0; k!=this->m_gsize_x; k++)
                {
                    input_file>>this->m_dmz_neut[i*this->m_n_states_neut+j][k];
   //                 cout<<i+1<<"_"<<j+1<<" -- "<<k<<"="<<this->m_dmz_neut[i*this->m_n_states_neut+j][k]<<endl;
                }
                input_file.close();
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
                for (int k=0; k!=this->m_gsize_x; k++)
                {
                    input_file>>this->m_dmx_cat[i*this->m_n_states_cat+j][k];
                }
                input_file.close();
            }
            else
            {
               cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit;
            }
            
            name_indenter.str("");
            name_indenter<<file_address<<"DMY"<<"_"<<i+1<<"_"<<j+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
                for (int k=0; k!=this->m_gsize_x; k++)
                {
                    input_file>>this->m_dmy_cat[i*this->m_n_states_cat+j][k];
                }
                input_file.close();
            }
            else
            {
               cout<<"ERROR DIPOLE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit;
            }
            
            name_indenter.str("");
            name_indenter<<file_address<<"DMZ"<<"_"<<i+1<<"_"<<j+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               input_file.seekg(0);
                for (int k=0; k!=this->m_gsize_x; k++)
                {
                    input_file>>this->m_dmz_cat[i*this->m_n_states_cat+j][k];
                }
                input_file.close();
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
void hamilton_matrix::set_PICE(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;
   double pot_vec[3];
   double elec_field[3];

   double temp;
   double Re_value;
   double Im_value;
   double pos;

   double *new_px=new double[this->m_n_states_cont];
   double *new_py=new double[this->m_n_states_cont];
   double *new_pz=new double[this->m_n_states_cont];
   double *trial_px=new double[this->m_n_states_cont];
   double *trial_py=new double[this->m_n_states_cont];
   double *trial_pz=new double[this->m_n_states_cont];

   ifstream input_file;
   
   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      for(int j=0;j!=this->m_n_states_cat;j++)
      {
         for(int x=0;x!=this->m_gsize_x;x++)
         {
            pos=this->m_xmin+x*(this->m_xmax-this->m_xmin)/this->m_gsize_x;
            name_indenter.str("");
            name_indenter<<file_address<<"pice"<<"_"<<i+1<<"_"<<j+1<<"_"<<x<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
               for(int k=0;k!=this->m_n_k;k++)
               {
                  for(int l=0;l!=this->m_n_angles;l++)
                  {
                     if(k==0)
                     {
                        input_file>>this->k_modulus[k];
                        input_file>>this->k_orientation[0][l];
                        input_file>>this->k_orientation[1][l];
                     }
                     else
                     {
                        input_file>>temp;
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
            else
            {
               cout<<"ERROR PICE FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit;
            }
         }
      }
   }

      //DETERMINING THE CLOSEST PAIR OF PLANE WAVES IN THE RECIPROCAL SPACE

                               for(int m=0;m!=this->m_n_states_cont;m++)
                               {   
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
   //END OF CLOSEST PAIR DETERMINATION
#pragma omp parallel for
   for(int t=0;t!=this->m_n_times;t++)
   {
#pragma omp parallel for
         for(int i=0;i!=this->m_n_states_cont;i++)
         {
            this->translation_vector[i][t]=i;
         }
      if(sqrt(pow(pot_vec[0],2)+pow(pot_vec[1],2)+pow(pot_vec[2],2))>min_distance/2)
      {
#pragma omp parallel for
         for(int i=0;i!=this->m_n_states_cont;i++)
         {
            potential_vector(t,pot_vec);
            ref_px[i]=this->k_modulus[i%this->m_n_angles]*sin(this->k_orientation[0][i-i%this->m_n_angles*this->m_n_angles])*cos(this->k_orientation[1][i-i%this->m_n_angles*this->m_n_angles]);
            ref_py[i]=this->k_modulus[i%this->m_n_angles]*sin(this->k_orientation[0][i-i%this->m_n_angles*this->m_n_angles])*sin(this->k_orientation[1][i-i%this->m_n_angles*this->m_n_angles]);
            ref_pz[i]=this->k_modulus[i%this->m_n_angles]*cos(this->k_orientation[0][i-i%this->m_n_angles*this->m_n_angles]);
            new_px[i]=ref_px[i]-pot_vec[0];
            new_py[i]=ref_py[i]-pot_vec[1];
            new_pz[i]=ref_pz[i]-pot_vec[2];
            trial_px[i]=ref_px[i];
            trial_py[i]=ref_py[i];
            trial_pz[i]=ref_pz[i];

#pragma omp parallel for
         for(int i=0;this->m_n_states_cont;i++)
         {
#pragma omp parallel for
            for(int j=0;j!=this->m_n_states_cont;j++)
            {
               if(pow(new_px[i]-ref_px[i],2)+pow(new_py[i]-ref_py[i],2)+pow(new_pz[i]-ref_pz[i],2) < pow(new_px[i]-trial_px[i],2)+pow(new_py[i]-trial_py[i],2)+pow(new_pz[i]-trial_pz[i],2) && i!=j)
               {
                  this->translation_vector[i][t]=j;
                  trial_px[i]=ref_px[j];
                  trial_py[i]=ref_py[j];
                  trial_pz[i]=ref_pz[j];
               }
            }
         }
      }
   }
   delete [] trial_px;
   delete [] trial_py;
   delete [] trial_pz;
   delete [] new_px;
   delete [] new_py;
   delete [] new_pz;
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_NAC(std::string file_address)
{
   using namespace std;
   stringstream name_indenter;
   string filename;

   ifstream input_file;

   for(int i=0;i!=m_n_states_neut;i++)
   {
      for(int j=i;j!=m_n_states_neut;j++)
      {   
            name_indenter.str("");
            name_indenter<<file_address<<"nac"<<"_"<<i+1<<"_"<<j+1<<".input";
            filename=name_indenter.str();
            input_file.open(filename.c_str());
            if(input_file.is_open())
            {
                for (int k=0; k!=this->m_gsize_x; k++)
                {
                   input_file>>m_NAC[i*this->m_n_states_neut+j][k];
                }
            }
            else
            {
               cout<<"ERROR NAC FILE NOT FOUND:"<<filename.c_str()<<endl<<"EXIT"<<endl;
               exit;
            }
      }
   }

}
double hamilton_matrix::kinetic_energy_matrix(int i,int j)
{
   return kinetic_energy[this->m_gsize_x*i+j];
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::pot_neut(int state_index,int grid_index)
{
   return this->m_pot_neut[state_index*(this->m_gsize_x)+grid_index];
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::rescale_pot(double min_pot)
{
   for(int m=0;m!=this->m_n_states_neut;m++)
   {
      for(int g=0;g!=this->m_gsize_x;g++)
      {
         this->m_pot_neut[m*(this->m_gsize_x)+g]-=min_pot;
      }   
   }
   for(int m=0;m!=this->m_n_states_cat;m++)
   {
      for(int g=0;g!=this->m_gsize_x;g++)
      {
         this->m_pot_cat[m*(this->m_gsize_x)+g]-=min_pot;
      }   
   }
}
//##########################################################################
//
//##########################################################################
//END OF HAMILTON MATRIX OBJECT INITIALIZATION AND SETTINGS

