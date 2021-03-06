
//WAVE FUNCTION OBJECT
//constructor of the wavefunction
//##########################################################################
//
//##########################################################################
wavefunction::wavefunction(int gsize_x,int tgsize_x,int n_states_neut,int n_states_cat,int n_states_cont)
{

   try
   {
      this->m_gsize_x=gsize_x;
      this->m_tgsize_x=tgsize_x;
      this->m_n_states_neut=n_states_neut;
      this->m_n_states_cat=n_states_cat;
      this->m_n_states_cont=n_states_cont;
      this->m_neut_part=new std::complex<double>[tgsize_x*n_states_neut];
      if(n_states_cat!=0 && n_states_cont != 0)
         this->m_cat_part=new std::complex<double>[tgsize_x*n_states_cat*n_states_cont];
      this->m_dipole_neut=new double[3];
      this->m_dipole_cat=new double[3];
   }
   catch(std::exception const& e)
   {
      std::cout<<" ERROR IN WAVE FUNCTION INITIATOR! : "<< e.what() <<std::endl;
      exit(EXIT_FAILURE);
   }
}
//##########################################################################
//
//##########################################################################
//Destructor of the wavefunction
wavefunction::~wavefunction()
{
   delete [] this->m_neut_part;
   if(this->m_n_states_cat!=0 && this->m_n_states_cont != 0)
      delete [] this->m_cat_part;
   delete [] this->m_dipole_neut;
   delete [] this->m_dipole_cat;

}
//##########################################################################
//
//##########################################################################
//set the neutral part of the wavefunction
void wavefunction::set_neut_psi(int state_index,int grid_index,std::complex<double> value)
{
   if(this->m_neut_part==NULL)
   {
      std::cout<<"FATAL ERROR: TRYING TO ADDRESS UNINITIALIZED WAVE FUNCTION VECTOR. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(grid_index < 0 || grid_index >= this->m_tgsize_x)
   {
      std::cout<<"FATAL ERROR: TRYING TO ADDRESS A GRID INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(state_index < 0 || state_index >= this->m_n_states_neut)
   {
      std::cout<<"FATAL ERROR: TRYING TO ADDRESS AN ELECTRONIC STATE INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else
   this->m_neut_part[state_index*this->m_tgsize_x+grid_index]=value;
   //The direct dimension of the neutral part is the grid size. The second dimension is the state index.
}
//##########################################################################
//
//##########################################################################
//set the cation part of the wavefunction
void wavefunction::set_cat_psi(int state_cat_index,int state_cont_index,int grid_index,std::complex<double> value)
{
   if(this->m_cat_part==NULL)
   {
      std::cout<<"FATAL ERROR: TRYING TO ADDRESS UNINITIALIZED WAVE FUNCTION VECTOR. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(grid_index < 0 || grid_index >= this->m_tgsize_x)
   {
      std::cout<<"FATAL ERROR: TRYING TO ADDRESS A GRID INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(state_cat_index < 0 || state_cat_index >= this->m_n_states_cat)
   {
      std::cout<<"FATAL ERROR: TRYING TO ADDRESS AN ELECTRONIC STATE INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(state_cont_index < 0 || state_cont_index >= this->m_n_states_cont)
   {
      std::cout<<"FATAL ERROR: TRYING TO ADDRESS AN CONTINUUM STATE INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else
   this->m_cat_part[state_cat_index*this->m_n_states_cont*this->m_tgsize_x+state_cont_index*this->m_tgsize_x+grid_index]=value;
   //The direct dimension of the cation part is the grid size.
   //The second dimension is the continuum state index
   //The third dimension is the electronic state index.
}
//##########################################################################
//
//##########################################################################
std::complex<double> wavefunction::show_neut_psi(int grid_index,int state_index)
{
/*   if(this->m_neut_part==NULL)
   {
      std::cout<<"FATAL ERROR: TRYING TO SHOW UNINITIALIZED WAVE FUNCTION VECTOR. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(grid_index < 0 || grid_index >= this->m_tgsize_x)
   {
      std::cout<<"FATAL ERROR: TRYING TO SHOW A GRID INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(state_index < 0 || state_index >= this->m_n_states_neut)
   {
      std::cout<<"FATAL ERROR: TRYING TO SHOW AN ELECTRONIC STATE INDEX OUT OF RANGE (NEUTRAL)."<<std::endl<<"STATE "<<state_index<<std::endl<<"EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);

   }*/ //DEBUG
   //else
   {
      //std::cout<<"probe show neut psi index"<<state_index*this->m_tgsize_x+grid_index<<std::endl;//DEBOGAGE
      return this->m_neut_part[state_index*this->m_tgsize_x+grid_index];
   }
   //The direct dimension of the neutral part is the grid size. The second dimension is the state index.
}
//##########################################################################
//
//##########################################################################
//set the cation part of the wavefunction
std::complex<double> wavefunction::show_cat_psi(int grid_index,int state_cat_index,int state_cont_index)
{
/*   if(this->m_cat_part==NULL)
   {
      std::cout<<"FATAL ERROR: TRYING TO SHOW UNINITIALIZED WAVE FUNCTION VECTOR. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(grid_index < 0 || grid_index >= this->m_tgsize_x)
   {
      std::cout<<"FATAL ERROR: TRYING TO SHOW A GRID INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(state_cat_index < 0 || state_cat_index >= this->m_n_states_cat)
   {
      std::cout<<"FATAL ERROR: TRYING TO SHOW AN ELECTRONIC STATE INDEX OUT OF RANGE (CATION). EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else if(state_cont_index < 0 || state_cont_index >= this->m_n_states_cont)
   {
      std::cout<<"FATAL ERROR: TRYING TO SHOW AN CONTINUUM STATE INDEX OUT OF RANGE. EXIT"<<std::endl;
      std::exit(EXIT_FAILURE);
   }
   else*/ // DEBUG
   return this->m_cat_part[state_cat_index*this->m_n_states_cont*this->m_tgsize_x+state_cont_index*this->m_tgsize_x+grid_index];
   //The direct dimension of the cation part is the grid size.
   //The second dimension is the electronic state index.
   //The third dimension is the continuum state index
}
//##########################################################################
//
//##########################################################################
int wavefunction::gsize_x()
{
   return this->m_gsize_x;
}
//##########################################################################
//
//##########################################################################
int wavefunction::tgsize_x()
{
   return this->m_tgsize_x;
}
//##########################################################################
//
//##########################################################################
int wavefunction::n_states_neut()
{
   return this->m_n_states_neut;
}
//##########################################################################
//
//##########################################################################
int wavefunction::n_states_cat()
{
   return this->m_n_states_cat;
}
//##########################################################################
//
//##########################################################################
int wavefunction::n_states_cont()
{
   return this->m_n_states_cont;
}
//##########################################################################
//
//##########################################################################
void wavefunction::initialize(hamilton_matrix* H)
{
   double *H_mat_gs=new double[(this->m_tgsize_x)*(this->m_tgsize_x)];
   double *d=new double[(this->m_tgsize_x)];
   double *e=new double [(this->m_tgsize_x)-1];
   std::complex<double> *tau=new std::complex<double> [(this->m_tgsize_x)-1];
   std::complex<double> *cmatrix=new std::complex<double> [(this->m_tgsize_x)*(this->m_tgsize_x)];
   double min_pot(10000);
   for(int i=0;i!=this->m_tgsize_x;i++)
   {
      if(H->pot_neut(0,i) < min_pot)
         min_pot=H->pot_neut(0,i);
   }
   std::cout<<"GS PES minimum is at "<<min_pot<<std::endl;
   H->rescale_pot(min_pot);

   for(int i=0;i!=this->m_tgsize_x;i++)
   {
      for(int j=0;j!=this->m_tgsize_x;j++)
      {
         H_mat_gs[i*(this->m_tgsize_x)+j]=0;
         if(i==j)
         H_mat_gs[i*(this->m_tgsize_x)+i]=H->pot_neut(0,i);
         H_mat_gs[i*(this->m_tgsize_x)+j]+=H->kinetic_energy_matrix(i,j);
         cmatrix[i*(this->m_tgsize_x)+j]=std::complex<double>(H_mat_gs[i*(this->m_tgsize_x)+j],0);
//         std::cout<<H_mat_gs[i*(this->m_tgsize_x)+j]<<"   ";
      }
//      std::cout<<std::endl;
   }

/*   std::cout<<"PES OF THE GS "<<std::endl;
   for(int i=0;i!=this->m_tgsize_x;i++)
   {
      std::cout<<H->pot_neut(0,i)<<std::endl;
   }*/

   
   std::cout<<"#############################"<<std::endl;

   std::cout<<"LAPACKE zhetrd returns "<<LAPACKE_zhetrd(LAPACK_ROW_MAJOR,'U',(this->m_tgsize_x),cmatrix,(this->m_tgsize_x),d,e,tau)<<std::endl;
   std::cout<<"LAPACKE zungtr returns "<<LAPACKE_zungtr(LAPACK_ROW_MAJOR,'U',(this->m_tgsize_x),cmatrix,(this->m_tgsize_x),tau)<<std::endl;
   std::cout<<"LAPACKE zstedc returns "<<LAPACKE_zstedc(LAPACK_ROW_MAJOR,'V',(this->m_tgsize_x),d,e,cmatrix,(this->m_tgsize_x))<<std::endl<<std::endl;

   std::cout<<"EIGENVALUES OF THE GS PES"<<std::endl;
   for(int i=0;i!=this->m_tgsize_x;i++)
   {
      std::cout<<d[i]<<std::endl;
   }
   std::cout<<"#############################"<<std::endl;

   for(int i=0;i!=this->m_tgsize_x;i++)
   {
      this->set_neut_psi(0,i,cmatrix[i*(this->m_tgsize_x)]);
   }
   /*for(int m=1;m!=this->m_n_states_neut;m++)
   {
      for(int i=0;i!=this->m_gsize_x;i++)
      {
         this->set_neut_psi(m,i,0);
      }
   }*/


   delete [] H_mat_gs;
   std::cout<<"Wavefunction Initialized ! "<<std::endl;
}
//##########################################################################
//
//##########################################################################
void wavefunction::projection_eigenstates(hamilton_matrix *H,int direction)
{
   std::complex<double> ctemp;
   double *reprojected_neut=new double[this->m_n_states_neut*this->m_tgsize_x];
   double *improjected_neut=new double[this->m_n_states_neut*this->m_tgsize_x];

   double *repartial_neut=new double[this->m_n_states_neut*this->m_tgsize_x];
   double *impartial_neut=new double[this->m_n_states_neut*this->m_tgsize_x];

   double *repartial_cat=new double[this->m_n_states_cat*this->m_tgsize_x];
   double *impartial_cat=new double[this->m_n_states_cat*this->m_tgsize_x];

   double *reprojected_cat=new double[this->m_n_states_cat*this->m_tgsize_x];
   double *improjected_cat=new double[this->m_n_states_cat*this->m_tgsize_x];

   double *eigenmat_neut=new double[this->m_n_states_neut*this->m_n_states_neut*this->m_tgsize_x*this->m_tgsize_x];
   double *eigenmat_cat=new double[this->m_n_states_cat*this->m_n_states_cat*this->m_tgsize_x*this->m_tgsize_x];


   if(direction==1) // THIS IS FOR PROJECTION OF WAVEFUNCTION ON EIGENSTATE BASIS SET
   {
//      std::cout<<"probe_neut"<<std::endl;
      wavefunction *proj_eigenstate=new wavefunction(this->m_gsize_x,this->m_tgsize_x,this->m_n_states_neut,this->m_n_states_cat,this->m_n_states_cont);

      H->eigenstates_matrix(0,eigenmat_neut);
      for(int n=0;n!=this->m_n_states_neut;n++)
      {
         for(int g=0;g!=this->m_tgsize_x;g++)
         {
            repartial_neut[n*this->m_tgsize_x+g]=real(this->m_neut_part[n*this->m_tgsize_x+g]);
            impartial_neut[n*this->m_tgsize_x+g]=imag(this->m_neut_part[n*this->m_tgsize_x+g]);
         }
      }

   
      cblas_dgemv(CblasRowMajor,CblasNoTrans,this->m_n_states_neut*this->m_tgsize_x,this->m_n_states_neut*this->m_tgsize_x,1,eigenmat_neut,this->m_n_states_neut*this->m_tgsize_x,repartial_neut,1,0,reprojected_neut,1);
      cblas_dgemv(CblasRowMajor,CblasNoTrans,this->m_n_states_neut*this->m_tgsize_x,this->m_n_states_neut*this->m_tgsize_x,1,eigenmat_neut,this->m_n_states_neut*this->m_tgsize_x,impartial_neut,1,0,improjected_neut,1);

      for(int n=0;n!=this->m_n_states_neut;n++)
      {
         for(int g=0;g!=this->m_tgsize_x;g++)
         {
            proj_eigenstate->set_neut_psi(n,g,std::complex<double>(reprojected_neut[n*this->m_tgsize_x+g],improjected_neut[n*this->m_tgsize_x+g]));
         }
      }
//      std::cout<<"probe_cat"<<std::endl;
//#pragma omp parallel for
      H->eigenstates_matrix(1,eigenmat_cat);
      //The eigenmat_cat only represents the eigenvectors of the cation PEC => m_tgsize_x*m_n_states_cat
      for(int k=0;k<this->m_n_states_cont;k++)
      {
         for(int n=0;n<this->m_n_states_cat;n++)
         {
            for(int g=0;g<this->m_tgsize_x;g++)
            {
               repartial_cat[n*this->m_tgsize_x+g]=real(this->m_cat_part[n*this->m_tgsize_x*this->m_n_states_cont+k*this->m_tgsize_x+g]);
               impartial_cat[n*this->m_tgsize_x+g]=imag(this->m_cat_part[n*this->m_tgsize_x*this->m_n_states_cont+k*this->m_tgsize_x+g]);
            }
         }
         cblas_dgemv(CblasRowMajor,CblasNoTrans,this->m_n_states_cat*this->m_tgsize_x,this->m_n_states_cat*this->m_tgsize_x,1,eigenmat_cat,this->m_n_states_cat*this->m_tgsize_x,repartial_cat,1,0,reprojected_cat,1);
         cblas_dgemv(CblasRowMajor,CblasNoTrans,this->m_n_states_cat*this->m_tgsize_x,this->m_n_states_cat*this->m_tgsize_x,1,eigenmat_cat,this->m_n_states_cat*this->m_tgsize_x,impartial_cat,1,0,improjected_cat,1);
         for(int n=0;n<this->m_n_states_cat;n++)
         {
            for(int g=0;g<this->m_tgsize_x;g++)
            {
               //std::cout<<n<<","<<k<<","<<g<<std::endl;
               //H->eigenstate(g,n,k,temp);
               proj_eigenstate->set_cat_psi(n,k,g,std::complex<double>(reprojected_cat[n*this->m_tgsize_x+g],improjected_cat[n*this->m_tgsize_x+g]));
            }
         }
      }
//      std::cout<<"probe_end"<<std::endl;
      this->set_wf(proj_eigenstate,1);
   delete proj_eigenstate;

   }
   
   else if(direction==-1)
   {

      wavefunction *proj_position=new wavefunction(this->m_gsize_x,this->m_tgsize_x,this->m_n_states_neut,this->m_n_states_cat,this->m_n_states_cont);

      H->eigenstates_matrix(0,eigenmat_neut);

      for(int n=0;n!=this->m_n_states_neut;n++)
      {
         for(int g=0;g!=this->m_tgsize_x;g++)
         {
            repartial_neut[n*this->m_tgsize_x+g]=real(this->m_neut_part[n*this->m_tgsize_x+g]);
            impartial_neut[n*this->m_tgsize_x+g]=imag(this->m_neut_part[n*this->m_tgsize_x+g]);
         }
      }

      cblas_dgemv(CblasRowMajor,CblasTrans,this->m_n_states_neut*this->m_tgsize_x,this->m_n_states_neut*this->m_tgsize_x,1,eigenmat_neut,this->m_n_states_neut*this->m_tgsize_x,repartial_neut,1,0,reprojected_neut,1);
      cblas_dgemv(CblasRowMajor,CblasTrans,this->m_n_states_neut*this->m_tgsize_x,this->m_n_states_neut*this->m_tgsize_x,1,eigenmat_neut,this->m_n_states_neut*this->m_tgsize_x,impartial_neut,1,0,improjected_neut,1);
      for(int n=0;n!=this->m_n_states_neut;n++)
      {
         for(int g=0;g!=this->m_tgsize_x;g++)
         {
            proj_position->set_neut_psi(n,g,std::complex<double>(reprojected_neut[n*this->m_tgsize_x+g],improjected_neut[n*this->m_tgsize_x+g]));
         }
      }

      H->eigenstates_matrix(1,eigenmat_cat);
      for(int k=0;k!=this->m_n_states_cont;k++)
      {
         for(int n=0;n<this->m_n_states_cat;n++)
         {
            for(int g=0;g<this->m_tgsize_x;g++)
            {
               repartial_cat[n*this->m_tgsize_x+g]=real(this->m_cat_part[n*this->m_tgsize_x*this->m_n_states_cont+k*this->m_tgsize_x+g]);
               impartial_cat[n*this->m_tgsize_x+g]=imag(this->m_cat_part[n*this->m_tgsize_x*this->m_n_states_cont+k*this->m_tgsize_x+g]);
            }
         }
         cblas_dgemv(CblasRowMajor,CblasTrans,this->m_n_states_cat*this->m_tgsize_x,this->m_n_states_cat*this->m_tgsize_x,1,eigenmat_cat,this->m_n_states_cat*this->m_tgsize_x,repartial_cat,1,0,reprojected_cat,1);
         cblas_dgemv(CblasRowMajor,CblasTrans,this->m_n_states_cat*this->m_tgsize_x,this->m_n_states_cat*this->m_tgsize_x,1,eigenmat_cat,this->m_n_states_cat*this->m_tgsize_x,impartial_cat,1,0,improjected_cat,1);
         for(int n=0;n<this->m_n_states_cat;n++)
         {
            for(int g=0;g<this->m_tgsize_x;g++)
            {
               //std::cout<<n<<","<<k<<","<<g<<std::endl;
               //H->eigenstate(g,n,k,temp);
               proj_position->set_cat_psi(n,k,g,std::complex<double>(reprojected_cat[n*this->m_tgsize_x+g],improjected_cat[n*this->m_tgsize_x+g]));
            }
         }
      }

      this->set_wf(proj_position,1);
      delete proj_position;
   }
   else
   {
      std::cout<<"ERROR. PROJECTION DIRECTION NOT SPECIFIED CORRECTLY. EXIT"<<std::endl;
      exit(EXIT_SUCCESS);
   }

   delete [] reprojected_neut;
   delete [] improjected_neut;
   delete [] reprojected_cat;
   delete [] improjected_cat;
   delete [] repartial_cat;
   delete [] impartial_cat;
   delete [] repartial_neut;
   delete [] impartial_neut;
   delete [] eigenmat_neut;
   delete [] eigenmat_cat;

}
//##########################################################################
//
//##########################################################################
/*double wavefunction::norm(hamilton_matrix *H)
{
   this->m_norm=real(this->dot_prod(this,H));
   return this->m_norm;
}*/
//##########################################################################
//
//##########################################################################
void wavefunction::set_norm(double value)
{
   this->m_norm=value;
}
//##########################################################################
//
//##########################################################################
void wavefunction::set_dipole(hamilton_matrix *H)
{
   std::complex<double> *vector=new std::complex<double> [3];
   int state_index;
   int grid_index;
   int state_index_cont;

   for(int i=0;i!=this->m_tgsize_x*this->m_n_states_neut;i++)
   {
      H->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);
      for(int s=0;s!=this->m_n_states_neut;s++)
      {
         vector[0]+=H->show_dm_neut(state_index,s,grid_index,0)*std::conj(this->m_neut_part[state_index*this->m_tgsize_x+grid_index])*this->m_neut_part[s*this->m_tgsize_x+grid_index];
         vector[1]+=H->show_dm_neut(state_index,s,grid_index,1)*std::conj(this->m_neut_part[state_index*this->m_tgsize_x+grid_index])*this->m_neut_part[s*this->m_tgsize_x+grid_index];
         vector[2]+=H->show_dm_neut(state_index,s,grid_index,2)*std::conj(this->m_neut_part[state_index*this->m_tgsize_x+grid_index])*this->m_neut_part[s*this->m_tgsize_x+grid_index];
      }
   }
   this->m_dipole_neut[0]=real(vector[0]);
   this->m_dipole_neut[1]=real(vector[1]);
   this->m_dipole_neut[2]=real(vector[2]);
   delete [] vector;
}
//##########################################################################
//
//##########################################################################
void wavefunction::show_dipole(double* vector,bool species)
{
   if(!species)
   {
      vector[0]=this->m_dipole_neut[0];
      vector[1]=this->m_dipole_neut[1];
      vector[2]=this->m_dipole_neut[2];
   }
   else
   {
      vector[0]=this->m_dipole_cat[0];
      vector[1]=this->m_dipole_cat[1];
      vector[2]=this->m_dipole_cat[2];
   }
}
//##########################################################################
//
//##########################################################################
void wavefunction::set_wf(wavefunction* x,bool cat)
{
   std::complex<double>* x_neut = new std::complex<double>[this->m_n_states_neut*this->m_tgsize_x];
   int i(0);
   if(cat)
   {
      std::complex<double>* x_cat = new std::complex<double>[this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x];
      x->wf_vec(x_neut,x_cat);
#pragma omp parallel for private(i)
      for(i=0;i<this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x;i++)
      {
         this->m_cat_part[i]=x_cat[i];
      }
//      cblas_zcopy(this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x,x_cat,1,this->m_cat_part,1);
      delete [] x_cat;
   }
   else
      x->wf_vec(x_neut,NULL);
#pragma omp parallel for private(i)
   for(i=0;i<this->m_n_states_neut*this->m_tgsize_x;i++)
   {
      this->m_neut_part[i]=x_neut[i];
   }
//   cblas_zcopy(this->m_n_states_neut*this->m_tgsize_x,x_neut,1,this->m_neut_part,1);

   delete [] x_neut;
}
//##########################################################################
//
//##########################################################################
void wavefunction::add_wf(std::complex<double>* a,wavefunction* y,bool cat)
{
   int i(0);
   if( y != NULL)
   {
      std::complex<double>* x_neut = new std::complex<double>[this->m_n_states_neut*this->m_tgsize_x];
      if(cat)
      {
         std::complex<double>* x_cat = new std::complex<double>[this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x];
         y->wf_vec(x_neut,x_cat);
#pragma omp parallel for private(i)
      for(i=0;i<this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x;i++)
      {
         this->m_cat_part[i]+=*a*x_cat[i];
      }
//         cblas_zaxpy(this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x,a,x_cat,1,this->m_cat_part,1);
         delete [] x_cat;
      }
      else
         y->wf_vec(x_neut,NULL);
#pragma omp parallel for private(i)
      for(i=0;i<this->m_n_states_neut*this->m_tgsize_x;i++)
      {
         this->m_neut_part[i]+=*a*x_neut[i];
      }
//      cblas_zaxpy(this->m_n_states_neut*this->m_tgsize_x,a,x_neut,1,this->m_neut_part,1);
      delete [] x_neut;
   }
   else
   {
      cblas_zscal(this->m_n_states_neut*this->m_tgsize_x,a,this->m_neut_part,1);
      if(cat)
         cblas_zscal(this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x,a,this->m_cat_part,1);
   }
}
//##########################################################################
//
//##########################################################################
void wavefunction::set_psi_elwise(int i,std::complex<double> val)
{
   int state_index_1(0);
   int grid_index_1(0);
   int state_index_cont_1(-1);
   if(i==0)
   {
      state_index_1 = 0;
      grid_index_1 = 0;
      state_index_cont_1 = -1;
   }
   else if(i < (this->n_states_neut())*(this->tgsize_x()))
   {
      grid_index_1 = int(i%this->tgsize_x());
      state_index_1 = int((i-grid_index_1)/this->tgsize_x());
      state_index_cont_1 = -1;
   }
   else
   {
      grid_index_1=int((i-this->n_states_neut()*this->tgsize_x())%this->tgsize_x());
      state_index_cont_1=int(((i-this->n_states_neut()*this->tgsize_x()-grid_index_1))/(this->tgsize_x()))%this->n_states_cont();
      state_index_1=int(((i-this->n_states_neut()*this->tgsize_x()-grid_index_1-(this->tgsize_x())*state_index_cont_1))/(this->tgsize_x()*this->n_states_cont()));
   }
//   std::cout<<grid_index_1<<","<<state_index_cont_1<<","<<state_index_1<<std::endl;

   if(state_index_cont_1 == -1)
      this->set_neut_psi(state_index_1,grid_index_1,val);
   else
      this->set_cat_psi(state_index_1,state_index_cont_1,grid_index_1,val);
   
}
//##########################################################################
//
//##########################################################################
std::complex<double> wavefunction::dot_prod(wavefunction* Bra, hamilton_matrix *H)
{
   std::complex<double>* Bra_neut=new std::complex<double> [this->m_tgsize_x*this->m_n_states_neut];
   std::complex<double> result_neut(0);
 std::complex<double> result_cat(0);
   if(this->m_n_states_cat != 0)
   {
      std::complex<double>* Bra_cat=new std::complex<double> [this->m_tgsize_x*this->m_n_states_cat*this->m_n_states_cont];
      Bra->wf_vec(Bra_neut,Bra_cat);
      /*
      for(int i=0;i!=this->m_n_states_cat;i++)
      {
         for(int j=0;j!=this->m_n_states_cont;j++)
         {
            for(int k=0;k!=this->m_tgsize_x;k++)
            {
               Bra_cat[i*this->m_n_states_cont*this->m_tgsize_x+j*this->m_tgsize_x+k] *= H->dk(j);
            }
         }
      }*/
      cblas_zdotc_sub(this->m_tgsize_x*this->m_n_states_cat*this->m_n_states_cont,Bra_cat,1,this->m_cat_part,1,&result_cat);
      delete [] Bra_cat;
   }
   else
      Bra->wf_vec(Bra_neut,NULL);

   cblas_zdotc_sub(this->m_tgsize_x*this->m_n_states_neut,Bra_neut,1,this->m_neut_part,1,&result_neut);

   delete [] Bra_neut;
   return result_neut+result_cat;
}
//##########################################################################
//
//##########################################################################
void wavefunction::wf_vec(std::complex<double>* neut_vec,std::complex<double>* cat_vec)
{
   cblas_zcopy(this->m_n_states_neut*this->m_tgsize_x,this->m_neut_part,1,neut_vec,1);
   if(cat_vec != NULL)
      cblas_zcopy(this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x,this->m_cat_part,1,cat_vec,1);
}
//##########################################################################
//
//##########################################################################
void wavefunction::matrix_prod(wavefunction** mat,wavefunction* ket,hamilton_matrix *H)
{
   for(int i=0;i!=this->m_n_states_neut*this->m_tgsize_x;i++)
   {
      this->set_psi_elwise(i,ket->dot_prod(mat[i],H));
   }
}
//##########################################################################
//
//##########################################################################
double wavefunction::state_pop(bool species,int state_index,hamilton_matrix* H)
{
   std::complex<double> value;

   if(!species)
   {
      cblas_zdotc_sub(this->m_tgsize_x,&this->m_neut_part[this->m_tgsize_x*state_index],1,&this->m_neut_part[this->m_tgsize_x*state_index],1,&value);
      return value.real();
   }
   else
   {
      if(H==NULL)
      {
         std::cout<<"ERROR ! Trying to compute cation state pop without giving the dk differential"<<std::endl;
         exit(EXIT_FAILURE);
      }
/*      std::complex<double> *vector=new std::complex<double>[this->m_n_states_cont*this->m_tgsize_x];

//#pragma omp parallel for private(state_index)
      for(int i=0;i<this->m_n_states_cont;i++)
      {
         for(int k=0;k<this->m_tgsize_x;k++)
         {
            vector[i*this->m_tgsize_x+k]=H->dk(i)*this->m_cat_part[this->m_tgsize_x*this->m_n_states_cont*state_index+this->m_tgsize_x*i+k];
         }
      }
      cblas_zdotc_sub(this->m_tgsize_x*this->m_n_states_cont,&this->m_cat_part[this->m_tgsize_x*this->m_n_states_cont*state_index],1,vector,1,&value);
  */    
      cblas_zdotc_sub(this->m_tgsize_x*this->m_n_states_cont,&this->m_cat_part[this->m_tgsize_x*this->m_n_states_cont*state_index],1,&this->m_cat_part[this->m_tgsize_x*this->m_n_states_cont*state_index],1,&value);
//      delete [] vector;

      return value.real();
   }
}
//##########################################################################
//
//##########################################################################
void wavefunction::save_wf(std::string file_loc)
{
    using namespace std;
    ofstream savestream;
    savestream.open(file_loc.c_str(),ios_base::trunc);
    savestream<<this->m_tgsize_x<<"\n"
    <<this->m_gsize_x<<"\n"
    <<this->m_n_states_neut<<"\n"
    <<this->m_n_states_cat<<"\n"
    <<this->m_n_states_cont<<"\n";

    for(int i=0;i!=this->m_tgsize_x*this->m_n_states_neut;i++)
    {
        savestream<<setprecision(15)<<this->m_neut_part[i]<<"\n";
    }
    if(this->m_n_states_cat!=0 && this->m_n_states_cont != 0)
    {
        for(int i=0;i!=this->m_tgsize_x*this->m_n_states_cat*this->m_n_states_cont;i++)
        {
            savestream<<setprecision(15)<<this->m_cat_part[i]<<"\n";
        }
    }
    savestream.close();
}
//##########################################################################
//
//##########################################################################
bool wavefunction::load_wf(std::string file_loc)
{
    using namespace std;
    ifstream loadstream;
    loadstream.open(file_loc.c_str());
    if(!loadstream.is_open())
    {
       std::cout<<"ERROR. RETART FILE CANNOT BE OPENED: "<<file_loc.c_str()<<std::endl<<"!!! EXIT"<<std::endl;
       exit(EXIT_FAILURE);
    }

    int tgsize_x(0);
    int gsize_x(0);
    int n_states_neut(0);
    int n_states_cat(0);
    int n_states_cont(0);
    loadstream>>tgsize_x;
    loadstream>>gsize_x;
    loadstream>>n_states_neut;
    loadstream>>n_states_cat;
    loadstream>>n_states_cont;

    if(tgsize_x != this->m_tgsize_x || gsize_x != this->m_gsize_x || n_states_neut != this->m_n_states_neut || n_states_cat != this->m_n_states_cat || n_states_cont != this->m_n_states_cont)
    {
       std::cout<<"ERROR WHEN LOADING WAVEFUNCTION UPON RESTART. RESTARTED WF PARAMETERS DO NOT CORRESPOND TO SIMULATION PARAMETERS. EXIT"<<std::endl;
//       exit(EXIT_FAILURE);
    }
//    else
    {
       for(int i=0;i!=this->m_tgsize_x*this->m_n_states_neut;i++)
       {
           loadstream>>this->m_neut_part[i];
       }
       if(this->m_n_states_cat!=0 && this->m_n_states_cont != 0)
       {
           for(int i=0;i!=this->m_tgsize_x*this->m_n_states_cat*this->m_n_states_cont;i++)
           {
               loadstream>>this->m_cat_part[i];
           }
       }
       loadstream.close();
    }
    return 1;
}
//##########################################################################
//
//##########################################################################
void wavefunction::analytic_propagation(hamilton_matrix *H,int timestep_number)
{
   wavefunction *new_state=new wavefunction(this->m_gsize_x,this->m_tgsize_x,this->m_n_states_neut,this->m_n_states_cat,this->m_n_states_cont);

   //NEUTRAL PART
   for( int n = 0;n != this->m_n_states_neut;n++ )
   {
      for(int g=0;g!=this->m_tgsize_x;g++)
      {
         new_state->set_neut_psi(n,g,this->m_neut_part[n*this->m_tgsize_x+g]*exp(std::complex<double>(0,-H->eigenvalue_neut(n,g)*timestep_number*H->h())));
      }
   }
   //CATION PART
   for(int n = 0;n != this->m_n_states_cat;n++ )
   {
      for(int k=0;k!=this->m_n_states_cont;k++)
      {
         for(int g=0;g!=this->m_tgsize_x;g++)
         {
            new_state->set_cat_psi(n,k,g,this->m_cat_part[n*this->m_n_states_cont*this->m_tgsize_x+k*this->m_tgsize_x+g ]*exp(std::complex<double>(0,-H->eigenvalue_cat(n,k,g)*timestep_number*H->h())));
         }
      }
   }
   this->set_wf(new_state);

   delete new_state;

}
//##########################################################################
//
//##########################################################################
void wavefunction::photoelectron_density(hamilton_matrix *H,double *cube_density,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double time_index)
{

   std::cout<<"Entering photoelectron density routine"<<std::endl;
   std::complex<double> *photoelectron_wavefunction=new std::complex<double>[nx*ny*nz]; 
   std::complex<double> eikr;
   double x(0);
   double y(0); 
   double z(0);

   int k_mod_index(0);
   int k_angle_index(0);

   double *pot_vec=new double [3];
   std::complex<double> basis_func_val(std::complex<double>(0,0));
//   std::complex<double> orthog(std::complex<double>(0,0));
   double *mo_val=new double [this->m_tgsize_x*H->pice_data_n_occ()];


   int g(0);
   int mm(0);
   int m(0);
   int n(0);
   int o(0);
   int k(0);
   int i(0);

   double *kx =new double[this->m_n_states_cont];
   double *ky =new double[this->m_n_states_cont];
   double *kz =new double[this->m_n_states_cont];
   double *kp =new double[this->m_n_states_cont];
   double *thet =new double[this->m_n_states_cont];
   double *phi =new double[this->m_n_states_cont];
   std::complex<double> *mo_ft=new std::complex<double> [this->m_n_states_cont*H->pice_data_n_occ()*this->m_tgsize_x];

   H->potential_vector(time_index,pot_vec);

   //FIRST, EVALUATE THE VARIABLES IN RECIPROCAL SPACE THAT DO NOT DEPEND ON THE POSITION OF THE PHOTOELECTRON IN DIRECT SPACE.

   try
   {
      for(k=0;k<this->m_n_states_cont;k++)
      {
         k_angle_index=int(k%H->n_angles());
         k_mod_index=int((k-k_angle_index)/H->n_angles());

//      std::cout<<" k = "<<k<<"==>"<<k_mod_index<<","<<k_angle_index<<std::endl; 
//      std::cout<<" pot vec =  "<<pot_vec[0]<<","<<pot_vec[1]<<","<<pot_vec[2]<<std::endl; 

         kx[k]=(H->k_mod_val(k_mod_index))*sin(H->k_spher_orient(0,k_angle_index))*cos(H->k_spher_orient(1,k_angle_index))+pot_vec[0];
         ky[k]=(H->k_mod_val(k_mod_index))*sin(H->k_spher_orient(0,k_angle_index))*sin(H->k_spher_orient(1,k_angle_index))+pot_vec[1];
         kz[k]=(H->k_mod_val(k_mod_index))*cos(H->k_spher_orient(0,k_angle_index))+pot_vec[2];
//      std::cout<<" kx = "<<kx[k]<<" ky = "<<ky[k]<<" kz = "<<kz[k]<<std::endl; 
         kp[k]=sqrt(pow(kx[k],2)+pow(ky[k],2)+pow(kz[k],2));
         if(kp[k]==0)
         {
            thet[k]=0;
            phi[k]=0;
         }
         else
         {
            thet[k]=acos(kz[k]/kp[k]);
            if( kx[k] == 0 && kx[k] >= 0 )
            {
               phi[k]=acos(-1)/2;
            }
            else if(kx[k] == 0 && ky[k] < 0)
            {
               phi[k]=3.*acos(-1)/2; 
            }
            else if(kx[k] < 0 && ky[k] == 0)
            {
               phi[k]=acos(-1); 
            }
            else
            {
               phi[k]=atan2(ky[k],kx[k]);
               if(phi[k]<0)
                  phi[k]+=2*acos(-1);
            }
         }
//      std::cout<<"SPHERICAL COORD: "<<kp[k]<<"::"<<thet[k]<<","<<phi[k]<<std::endl;
         #pragma omp parallel for private(g,mm) shared(k)
         for(g=0;g < this->m_tgsize_x;g++)
         {
            for(mm=0;mm < H->pice_data_n_occ();mm++)
            {
               mo_ft[k*H->pice_data_n_occ()*this->m_tgsize_x+g*H->pice_data_n_occ()+mm]=H->pice_data_mo_ft_value(kp[k],thet[k],phi[k],mm,g);
//            std::cout<<"MO FT VALUES "<<k<<"::"<<g<<","<<mm<<" == "<<mo_ft[k*H->pice_data_n_occ()*this->m_tgsize_x+g*H->pice_data_n_occ()+mm]<<std::endl;
               if(isnan(real(mo_ft[k*H->pice_data_n_occ()*this->m_tgsize_x+g*H->pice_data_n_occ()+mm]))||isnan(imag(mo_ft[k*H->pice_data_n_occ()*this->m_tgsize_x+g*H->pice_data_n_occ()+mm])))
               {
                  std::cout<<"error in the reading of MO FT values !! EXIT"<<std::endl;
                  exit(EXIT_FAILURE);
               }
            }
         }
      }
   }
   catch(std::exception const& e)
   {
      std::cout<<" ERROR IN COMPUTATION OF MO FT : "<< e.what() <<std::endl;
      exit(EXIT_FAILURE);
   }
   std::cout<<"We made it... Now evaluating mo values in direct space"<<std::endl;
//      exit(EXIT_SUCCESS);

   //THEN, USE THEM IN THE COMPUTATION OF THE DIRECT SPACE WAVEFUNCTION. 
   try
   {
      for(m=0;m<nx;m++)
      {
         x=xmin+m*(xmax-xmin)/nx;
         for(n=0;n<ny;n++)
         {
            y=ymin+n*(ymax-ymin)/ny;
            for(o=0;o<nz;o++)
            {
               z=zmin+o*(zmax-zmin)/nz;
 //              std::cout<<" m = "<<m<<"/"<<nx<<" ; n = "<<n<<"/"<<ny<<" ; o = "<<o<<"/"<<nz<<std::endl; 

               photoelectron_wavefunction[m*ny*nz+n*nz+o]=0;
               //FOR EVERY POINT IN SPACE, BUILD AN ARRAY THAT CONTAINS THE VALUE OF EVERY MO
//#pragma omp parallel for private(g,mm) shared(m,n,o)
               for(g=0;g < this->m_tgsize_x;g++)
               {
                  for(mm=0;mm < H->pice_data_n_occ();mm++)
                  {
                      mo_val[g*H->pice_data_n_occ()+mm]=H->pice_data_mo_value(x,y,z,mm,g);

                      if(isnan(mo_val[g*H->pice_data_n_occ()+mm]))
                      {
                         std::cout<<"error in the reading of MO values !! EXIT"<<std::endl;
                         exit(EXIT_FAILURE);
                      }
                  }
               }
               for(i = 0;i < this->m_n_states_cat;i++ )
               {
                  for(k = 0 ; k < this->m_n_states_cont;k++)
                  {
//                     std::cout<<" k = "<<k<<std::endl; 
                     eikr=exp(std::complex<double>(0,-(x*kx[k]+y*ky[k]+z*kz[k])));
#pragma omp parallel for private(g,mm,basis_func_val) shared(eikr,m,n,o,i,k)
                     for(g=0;g < this->m_tgsize_x;g++)
                     {
//                        orthog=std::complex<double>(0,0);
                        basis_func_val=eikr;
                        
                        for(mm=0;mm < H->pice_data_n_occ();mm++)
                        {
                           //orthog+=H->pice_data_mo_value(x,y,z,mm,g)*mo_ft[k*H->pice_data_n_occ()*this->m_tgsize_x+g*H->pice_data_n_occ()+mm];
                           basis_func_val-=mo_val[g*H->pice_data_n_occ()+mm]*mo_ft[k*H->pice_data_n_occ()*this->m_tgsize_x+g*H->pice_data_n_occ()+mm];
                        }
//                        basis_func_val=eikr-orthog;/*exp(-ikr)-SUM_i^MO[FT_MO_i*MO_i]*/
                        photoelectron_wavefunction[m*ny*nz+n*nz+o]+=this->m_cat_part[i*this->m_n_states_cont*this->m_tgsize_x+k*this->m_tgsize_x+g]*basis_func_val*H->dk(k);
//                        std::cout<<k<<","<<g<<"-"<<". vals :  "<<this->m_cat_part[i*this->m_n_states_cont*this->m_tgsize_x+k*this->m_tgsize_x+g]<<" * "<<basis_func_val<<std::endl;

                     }
                  }
               }
               cube_density[m*ny*nz+n*nz+o]=pow(abs(photoelectron_wavefunction[m*ny*nz+n*nz+o]),2);
               if(cube_density[m*ny*nz+n*nz+o] < 0)
               {
                  std::cout<<"FATAL ERROR !!! cube density > 0 :"<<cube_density[m*ny*nz+n*nz+o]<<std::endl;
                  exit(EXIT_FAILURE);
               }
               else if(isnan(cube_density[m*ny*nz+n*nz+o]))
               {
                  std::cout<<"FATAL ERROR !!! Cube density is Not A Number! "<<std::endl;
                  exit(EXIT_FAILURE);
               }
//               std::cout<<"val : "<<cube_density[m*ny*nz+n*nz+o]<<std::endl;
            }
         }
      }
   }
   catch(std::exception const& e)
   {
      std::cout<<" ERROR IN COMPUTATION OF MO VALUE : "<< e.what() <<std::endl;
      std::cout<<" INDEXES ARE m = "<<m<<"; n = "<<n<<"; o = "<<o<<" ; mm = "<<mm<<"; i = "<<i<<" ; g = "<<g<<"; k = "<<k<<std::endl; 
      exit(EXIT_FAILURE);
   }
   std::cout<<"Done. Cleaning memory"<<std::endl;
      delete [] photoelectron_wavefunction;
      delete [] mo_val;
      delete [] kx;
      delete [] ky;
      delete [] kz;
      delete [] thet;
      delete [] phi;
      delete [] mo_ft;
      delete [] pot_vec;

   std::cout<<"Leaving photoelectron density routine"<<std::endl;
   return;
}
//##########################################################################
//
//##########################################################################
//END OF THE WAVE FUNCTION OBJECT
