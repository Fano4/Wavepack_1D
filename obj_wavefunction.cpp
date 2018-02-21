
//WAVE FUNCTION OBJECT
//constructor of the wavefunction
//##########################################################################
//
//##########################################################################
wavefunction::wavefunction(int gsize_x,int tgsize_x,int n_states_neut,int n_states_cat,int n_states_cont)
{

   this->m_gsize_x=gsize_x;
   this->m_tgsize_x=tgsize_x;
   this->m_n_states_neut=n_states_neut;
   this->m_n_states_cat=n_states_cat;
   this->m_n_states_cont=n_states_cont;
   this->m_neut_part=new std::complex<double>[tgsize_x*n_states_neut];
   if(n_states_cat!=0 && n_states_cont != 0)
      this->m_cat_part=new std::complex<double>[tgsize_x*n_states_cat*n_states_cont];
   this->m_dipole=new double[3];
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
   delete [] this->m_dipole;
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
   if(this->m_neut_part==NULL)
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

   }
   else
   return this->m_neut_part[state_index*this->m_tgsize_x+grid_index];
   //The direct dimension of the neutral part is the grid size. The second dimension is the state index.
}
//##########################################################################
//
//##########################################################################
//set the cation part of the wavefunction
std::complex<double> wavefunction::show_cat_psi(int grid_index,int state_cat_index,int state_cont_index)
{
   if(this->m_cat_part==NULL)
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
   else
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
   double min_pot(0);
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
         if(i==j)
         H_mat_gs[i*(this->m_tgsize_x)+i]=H->pot_neut(0,i);
         H_mat_gs[i*(this->m_tgsize_x)+j]+=H->kinetic_energy_matrix(i,j);
         cmatrix[i*(this->m_tgsize_x)+j]=std::complex<double>(H_mat_gs[i*(this->m_tgsize_x)+j],0);
         //std::cout<<H_mat_gs[i*(this->m_tgsize_x)+j]<<"   ";
      }
      //std::cout<<std::endl;
   }
/*
   std::cout<<"PES OF THE GS "<<std::endl;
   for(int i=0;i!=this->m_tgsize_x;i++)
   {
      std::cout<<H->pot_neut(0,i)<<std::endl;
   }
   */
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
   return ;
   double *vector=new double[3];
   int state_index;
   int grid_index;
   int state_index_cont;

   for(int i=0;i!=this->m_gsize_x*this->m_n_states_neut;i++)
   {
      H->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);
      for(int s=0;s!=this->m_n_states_neut;s++)
      {
//         vector[0]+=H->m_dmx_neut[state_index*this->m_n_states_neut+s][grid_index]*;
      }
      for(int j=0;j!=this->m_gsize_x*this->m_n_states_cat*this->m_n_states_cont;j++)
      {
         continue;
      }
   }
   for(int i=0;i!=this->m_gsize_x*this->m_n_states_cat*this->m_n_states_cont;i++)
   {
      for(int j=0;j!=this->m_gsize_x*this->m_n_states_neut;j++)
      {
         continue;
      }
      for(int j=0;j!=this->m_gsize_x*this->m_n_states_cat*this->m_n_states_cont;j++)
      {
         continue;
      }
   }
   //this->m_dipole[0]=vector[0];
   //this->m_dipole[1]=vector[1];
   //this->m_dipole[2]=vector[2];
}
//##########################################################################
//
//##########################################################################
void wavefunction::show_dipole(double* vector)
{
   vector[0]=this->m_dipole[0];
   vector[1]=this->m_dipole[1];
   vector[2]=this->m_dipole[2];
}
//##########################################################################
//
//##########################################################################
void wavefunction::set_wf(wavefunction* x,bool cat)
{
   std::complex<double>* x_neut = new std::complex<double>[this->m_n_states_neut*this->m_tgsize_x];
   if(cat)
   {
      std::complex<double>* x_cat = new std::complex<double>[this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x];
      x->wf_vec(x_neut,x_cat);
      cblas_zcopy(this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x,x_cat,1,this->m_cat_part,1);
      delete [] x_cat;
   }
   else
      x->wf_vec(x_neut,NULL);
   cblas_zcopy(this->m_n_states_neut*this->m_tgsize_x,x_neut,1,this->m_neut_part,1);

   delete [] x_neut;
}
//##########################################################################
//
//##########################################################################
void wavefunction::add_wf(std::complex<double>* a,wavefunction* y,bool cat)
{
   if( y != NULL)
   {
      std::complex<double>* x_neut = new std::complex<double>[this->m_n_states_neut*this->m_tgsize_x];
      if(cat)
      {
         std::complex<double>* x_cat = new std::complex<double>[this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x];
         y->wf_vec(x_neut,x_cat);
         cblas_zaxpy(this->m_n_states_cat*this->m_n_states_cont*this->m_tgsize_x,a,x_cat,1,this->m_cat_part,1);
         delete [] x_cat;
      }
      else
         y->wf_vec(x_neut,NULL);
      cblas_zaxpy(this->m_n_states_neut*this->m_tgsize_x,a,x_neut,1,this->m_neut_part,1);
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
   else if(i < (this->n_states_neut())*(this->gsize_x()))
   {
      grid_index_1 = int(i%this->gsize_x());
      state_index_1 = int((i-grid_index_1)/this->gsize_x());
      state_index_cont_1 = -1;
   }
   else
   {
      grid_index_1=int((i-this->n_states_neut()*this->gsize_x())%this->gsize_x());
      state_index_cont_1=int(((i-this->n_states_neut()*this->gsize_x()-grid_index_1))/(this->gsize_x()))%this->n_states_cont();
      state_index_1=int(((i-this->n_states_neut()*this->gsize_x()-grid_index_1-(this->gsize_x())*state_index_cont_1))/(this->gsize_x()*this->n_states_cont()));
   }
   std::cout<<grid_index_1<<","<<state_index_cont_1<<","<<state_index_1<<std::endl;

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
//END OF THE WAVE FUNCTION OBJECT
