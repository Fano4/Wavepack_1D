

//THIS FILE CONTAINS THE ROUTINES ASSOCIATED TO THE COMPUTATION RELATED TO THE HAMILTON MATRIX OBJECT.
//
//TABLE OF CONTENT:
//
//HAMILTON MATRIX OBJECT COMPUTATION


//HAMILTON MATRIX OBJECT COMPUTATION
//##########################################################################
//
// /        |               \
// |  N-N   |      N-C      |
// |        |               |
// |------------------------|
// |        |               |
// |        |               |
// |        |               |
// |  C-N   |      C-C      |
// |        |               |
// |        |               |
// \        |               /
//##########################################################################
std::complex<double> hamilton_matrix::hamilt_element(double time_index,int i, int j)
{
   int state_index_1(0);
   int grid_index_1(0);
   int state_index_cont_1(-1);
   int state_index_2(0);
   int grid_index_2(0);
   int state_index_cont_2(-1);
   double elec_field[3];
   double pot_vector[3];

   this->show_indexes(i,j,&state_index_1,&grid_index_1,&state_index_cont_1,&state_index_2,&grid_index_2,&state_index_cont_2);

   this-> electric_field(time_index,elec_field);
   this-> potential_vector(time_index,pot_vector);

   if(state_index_cont_1 == -1 && state_index_cont_2 == -1)//N-N section
   {
      if(state_index_1==state_index_2)//diagonal element in the electronic state basis
      {
         if(grid_index_1==grid_index_2)//diagonal in the grid basis
         {
            return std::complex<double>(this->kinetic_energy[grid_index_1*this->m_gsize_x+grid_index_2]+this->m_pot_neut[state_index_1*this->m_gsize_x+grid_index_1]-this->m_dmx_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[0]-this->m_dmy_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[1]-this->m_dmz_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[2],0); //diagonal term of kinetic energy + potential energy + permanent dipole interaction
         }
         else // non-diagonal term of kinetic energy
         {
            return std::complex<double>(this->kinetic_energy[grid_index_1*this->m_gsize_x+grid_index_2],0);
         }
      }
      else //coupling elements between electronic states of the neutral
      {
         if(grid_index_1==grid_index_2)//diagonal in the grid basis
         {
            return std::complex<double>(-this->m_dmx_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[0]-this->m_dmy_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[1]-this->m_dmz_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[2],0); //dipole coupling
         }
         else //NACME
         {
            if(state_index_1 < state_index_2)//upper triangle of the matrix
            {
               return std::complex<double>(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*this->derivative_matrix[grid_index_1*this->m_gsize_x+grid_index_2]/this->m_mass,0);//Non adiabatic coupling
            }
            else //lower triangle of the matrix
            {
               return std::complex<double>(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*this->derivative_matrix[grid_index_1*this->m_gsize_x+grid_index_2]/this->m_mass,0);//Non adiabatic coupling
            }
         }
      }
   }

   else if( (state_index_cont_1 ==-1 || state_index_cont_2==-1) && !(state_index_cont_2==-1 && state_index_cont_1 ==-1) )//N-C or C-N section => PICE
   {
      if(state_index_cont_1 ==-1 && grid_index_1==grid_index_2)//N-C
      {
         return -(this->m_PICE_x[state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+(this->translation_vector[state_index_cont_2][int(time_index)])][grid_index_1])*elec_field[0]-(this->m_PICE_y[state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+(this->translation_vector[state_index_cont_2][int(time_index)])][grid_index_1])*elec_field[1]-(this->m_PICE_z[state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+(this->translation_vector[state_index_cont_2][int(time_index)])][grid_index_1])*elec_field[2];//PICE dipole interaction
      }
      else if(grid_index_1==grid_index_2) // C-N
      {
         return -std::conj(this->m_PICE_x[state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+this->translation_vector[state_index_cont_1][int(time_index)]][grid_index_1])*elec_field[0]-std::conj(this->m_PICE_y[state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+this->translation_vector[state_index_cont_1][int(time_index)]][grid_index_1])*elec_field[1]-std::conj(this->m_PICE_z[state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+this->translation_vector[state_index_cont_1][int(time_index)]][grid_index_1])*elec_field[2];//PICE dipole interaction
      }
      else
         return 0;//non-diagonal in grid index
   }
   else//C-C section
   {
      if(state_index_cont_1==state_index_cont_2)
      {
         if(state_index_1==state_index_2)//diagonal element in the electronic state basis
         {
            if(grid_index_1==grid_index_2)//diagonal in the grid basis
            {
               return std::complex<double>(0.5*(pow(this->k_modulus[state_index_cont_1%this->m_n_k]*sin(this->k_orientation[0][state_index_cont_1-state_index_cont_1%this->m_n_k])*cos(this->k_orientation[1][state_index_cont_1-state_index_cont_1%this->m_n_k]),2)+pow(k_modulus[state_index_cont_1%this->m_n_k]*sin(this->k_orientation[0][state_index_cont_1-state_index_cont_1%this->m_n_k])*sin(this->k_orientation[1][state_index_cont_1-state_index_cont_1%this->m_n_k]),2)+this->k_modulus[state_index_cont_1%this->m_n_k]*cos(this->k_orientation[0][state_index_cont_1-state_index_cont_1%this->m_n_k]))+this->kinetic_energy[grid_index_1*this->m_gsize_x+grid_index_2]+this->m_pot_cat[state_index_1*this->m_gsize_x+grid_index_1]-this->m_dmx_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[0]-this->m_dmy_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[1]-this->m_dmz_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[2],0); //diagonal term of kinetic energy + potential energy + permanent dipole interaction
            }
            else // non-diagonal term of kinetic energy
            {
               return std::complex<double>(this->kinetic_energy[grid_index_1*this->m_gsize_x+grid_index_2],0);
            }
         }
         else //coupling elements between electronic states of the cation
         {
            if(grid_index_1==grid_index_2)//diagonal in the grid basis
            {
               return std::complex<double>(-this->m_dmx_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[0]-this->m_dmy_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[1]-this->m_dmz_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[2]); //dipole coupling
            }
            else
            {
               if(state_index_1 < state_index_2)//upper triangle of the matrix
               {
                     return 0;//std::complex<double>(0,-m_NAC[state_index_1*n_states_neut+state_index_2][grid_index_1]*derivative_matrix[grid_index_1*gsize_x+grid_index_2]/m_mass);//Non adiabatic coupling
               }
               else //lower triangle of the matrix
               {
                  return 0;//std::complex<double>(0,m_NAC[state_index_1*n_states_neut+state_index_2][grid_index_1]*derivative_matrix[grid_index_1*gsize_x+grid_index_2]/m_mass);//Non adiabatic coupling
               }
            }
         }
      }
      else
         return 0;
   }
}
//##########################################################################
//
//This function computes the energy and the norm and the total dipole of the wavefunction that is given for this hamiltonian. 
//It is for gaining time that the norm and dipole are computed here because it is usually reauired at the same time as energy.
//
//##########################################################################
double hamilton_matrix::energy(wavefunction* Psi,double time_index)
{
   using namespace std;
   complex<double> ec(complex<double>(0,0));
   double e(0);
   double cnorm(0);
   std::complex<double> *dipvec=new std::complex<double> [3];
   int index_1(0);
   int index_2(0);
   int state_index;
   int grid_index;
   int state_index_cont;
   int state_index2;
   int grid_index2;
   int state_index_cont2;


//#pragma omp parallel for
      for(int i=0;i!=(Psi->n_states_neut())*(Psi->gsize_x());i++)
      {
         int j(0);
         
         j=(i%Psi->gsize_x()-2)*bool(i%Psi->gsize_x()-2 >= 0);
         this->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);

         cnorm+=(conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index))).real();

         for(int s=0;s!=Psi->n_states_neut();s++)
         {
            while( (j%Psi->gsize_x()) <= i%Psi->gsize_x()+2)
            {
               this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
               ec+=(this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index2,state_index2));

               if(i%Psi->gsize_x() == j%Psi->gsize_x())
               {
                  dipvec[0]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmx_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
                  dipvec[1]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmy_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
                  dipvec[2]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmz_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
               }

               if(!(j%Psi->gsize_x() == Psi->gsize_x()-1))
                  j++;
               else
               {
                  j++;
                  break;
               }
            }
            j += Psi->gsize_x() - 5 + 2*bool(i%Psi->gsize_x() == 0) + bool (i%Psi->gsize_x() == 1) + 2*bool(i%Psi->gsize_x() == Psi->gsize_x()-1) + bool (i%Psi->gsize_x() ==Psi->gsize_x()-2);
         }
         j+=2;
         for(int s=0;s!=Psi->n_states_cont()*Psi->n_states_cat();s++)
         {
            this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
            ec+=(this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_cat_psi(grid_index2,state_index2,state_index_cont2));
            j += Psi->gsize_x();
               if(i%Psi->gsize_x() == j%Psi->gsize_x())
               {
                  dipvec[0]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_PICE_x[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
                  dipvec[1]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_PICE_y[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
                  dipvec[2]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_PICE_z[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
               }
         }

      }

//#pragma omp parallel for
      for(int i=(Psi->n_states_neut())*(Psi->gsize_x());i!=(Psi->n_states_neut())*(Psi->gsize_x())+(Psi->n_states_cat())*(Psi->n_states_cont())*(Psi->gsize_x());i++)
      {
         int j(0);
         j=i%Psi->gsize_x();

         this->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);

         cnorm+=(conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index,state_index_cont))).real();

         for(int s=0;s!=Psi->n_states_neut();s++)
         {
            this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
               if(i%Psi->gsize_x() == j%Psi->gsize_x())
               {
                  dipvec[0]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_neut_psi(grid_index,state_index2))*conj(this->m_PICE_x[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
                  dipvec[1]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_neut_psi(grid_index,state_index2))*conj(this->m_PICE_y[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
                  dipvec[2]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_neut_psi(grid_index,state_index2))*conj(this->m_PICE_z[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
               }
            ec+=(this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_neut_psi(grid_index2,state_index2));
            j += Psi->gsize_x();
         }
         for(int s=0;s!=Psi->n_states_cat()*Psi->n_states_cont();s++)
         {
            while(j%Psi->gsize_x() <= (i+2)%Psi->gsize_x())
            {
               this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
               ec+=(this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index2,state_index2,state_index_cont2));
               if(i%Psi->gsize_x() == j%Psi->gsize_x())
               {
                  dipvec[0]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmx_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
                  dipvec[1]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmy_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
                  dipvec[2]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmz_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
               }

               if(!(j%Psi->gsize_x() == Psi->gsize_x()-1))
                  j++;
               else
               {
                  j++;
                  break;
               }
            }
            j += Psi->gsize_x() - 5 + 2*bool(i%Psi->gsize_x() == 0) + bool (i%Psi->gsize_x() == 1) + 2*bool(i%Psi->gsize_x() == Psi->gsize_x()-1) + bool (i%Psi->gsize_x() ==Psi->gsize_x()-2);
         }
      }
/*
   for(int g=0;g!=this->m_gsize_x;g++)
   {
      //N-N and C-C terms
      for(int h=0;h!=this->m_gsize_x;h++)
      {
         //N-N terms
         for(int m=0;m!=this->m_n_states_neut;m++)
         {
            if( g == h )
               cnorm+=(conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(g,m))).real();

            for(int n=0;n!=this->m_n_states_neut;n++)
            {
               index_1=m*this->m_gsize_x+g;
               index_2=n*this->m_gsize_x+h;

               if( g == h )
               {
                  dipvec[0]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(g,n))*(this->m_dmx_neut[m*this->m_n_states_neut+n][g]);
                  dipvec[1]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(g,n))*(this->m_dmy_neut[m*this->m_n_states_neut+n][g]);
                  dipvec[2]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(g,n))*(this->m_dmz_neut[m*this->m_n_states_neut+n][g]);
             // *     if(m == n)
                  {
                     dipvec[0]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(g,n))*nucl_dip_x;
                     dipvec[1]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(g,n))*nucl_dip_y;
                     dipvec[2]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(g,n))*nucl_dip_z;
                  }* /
               }
               ec+=(this->hamilt_element(time_index,index_1,index_2))*conj(Psi->show_neut_psi(g,m))*(Psi->show_neut_psi(h,n));
            }
         }
         //C-C terms
         for(int m=0;m!=this->m_n_states_cat;m++)
         {
            for(int n=0;n!=this->m_n_states_cat;n++)
            {
               for(int j=0;j!=this->m_n_states_cont;j++)
               {
                  index_1=this->m_n_states_neut*this->m_gsize_x+m*this->m_n_states_cont*this->m_gsize_x+j*this->m_gsize_x+g;
                  index_2=this->m_n_states_neut*this->m_gsize_x+n*this->m_n_states_cont*this->m_gsize_x+j*this->m_gsize_x+h;
                  if( g == h && m == n)
                  cnorm+=(conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(g,m,j))).real();
                  if( g == h )
                  {
                     dipvec[0]+=conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(g,n,j))*(this->m_dmx_cat[m*this->m_n_states_neut+n][g]);
                     dipvec[1]+=conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(g,n,j))*(this->m_dmy_cat[m*this->m_n_states_neut+n][g]);
                     dipvec[2]+=conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(g,n,j))*(this->m_dmz_cat[m*this->m_n_states_neut+n][g]);
                    / * if(m == n)
                     {
                        dipvec[0]+=conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(g,n,j))*nucl_dip_x;
                        dipvec[1]+=conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(g,n,j))*nucl_dip_y;
                        dipvec[2]+=conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(g,n,j))*nucl_dip_z;
                     }* /
                  }

                  ec+=(this->hamilt_element(time_index,index_1,index_2))*conj(Psi->show_cat_psi(g,m,j))*(Psi->show_cat_psi(h,n,j));
               }
            }
         }
      }
      //C-N and N-C terms
      for(int m=0;m!=this->m_n_states_neut;m++)
      {
         for(int n=0;n!=this->m_n_states_cat;n++)
         {
            for(int j=0;j!=this->m_n_states_cont;j++)
            {
               index_1=m*this->m_gsize_x+g;
               index_2=this->m_n_states_neut*this->m_gsize_x+n*this->m_n_states_cont*this->m_gsize_x+j*this->m_gsize_x+this->m_h;
               ec+=(this->hamilt_element(time_index,index_1,index_2))*conj(Psi->show_neut_psi(g,m))*(Psi->show_cat_psi(g,n,j));
               ec+=(this->hamilt_element(time_index,index_2,index_1))*conj(Psi->show_cat_psi(g,n,j))*(Psi->show_neut_psi(g,m));
               dipvec[0]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_cat_psi(g,n,j))*conj(this->m_PICE_x[m*this->m_n_states_neut+n][g]);
               dipvec[1]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_cat_psi(g,n,j))*conj(this->m_PICE_y[m*this->m_n_states_neut+n][g]);
               dipvec[2]+=conj(Psi->show_neut_psi(g,m))*(Psi->show_cat_psi(g,n,j))*conj(this->m_PICE_z[m*this->m_n_states_neut+n][g]);
               dipvec[0]+=(Psi->show_neut_psi(g,m))*conj(Psi->show_cat_psi(g,n,j))*(this->m_PICE_x[m*this->m_n_states_neut+n][g]);
               dipvec[1]+=(Psi->show_neut_psi(g,m))*conj(Psi->show_cat_psi(g,n,j))*(this->m_PICE_y[m*this->m_n_states_neut+n][g]);
               dipvec[2]+=(Psi->show_neut_psi(g,m))*conj(Psi->show_cat_psi(g,n,j))*(this->m_PICE_z[m*this->m_n_states_neut+n][g]);
            }
         }
      }
   }
*/
   double* rdipvec = new double [3];

   rdipvec[0]=dipvec[0].real();
   rdipvec[1]=dipvec[1].real();
   rdipvec[2]=dipvec[2].real();

   delete [] dipvec;

   e=ec.real();
   Psi->set_norm(cnorm);
   Psi->set_dipole(rdipvec);
   delete [] rdipvec;

   return e;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::h()
{
   return this->m_h;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::efield_thresh()
{
   return this->m_efield_thresh;
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::show_indexes(int index1,int index2,int *state_index_1,int *grid_index_1,int *state_index_cont_1,int *state_index_2,int *grid_index_2,int *state_index_cont_2)
{
   if(index1==0)
   {
      *state_index_1 = 0;
      *grid_index_1 = 0;
      *state_index_cont_1 = -1;
   }
   else if(index1 < (this->m_n_states_neut)*(this->m_gsize_x))
   {
      *grid_index_1 = int(index1%this->m_gsize_x);
      *state_index_1 = int((index1-*grid_index_1)/this->m_gsize_x);
      *state_index_cont_1 = -1;
   }
   else
   {
      *state_index_1=int((index1-this->m_n_states_neut*this->m_gsize_x)/(this->m_gsize_x*this->m_n_states_cont));
      *state_index_cont_1=int(((index1-this->m_n_states_neut*this->m_gsize_x)-*state_index_1*(this->m_gsize_x*this->m_n_states_cont))/this->m_gsize_x);
      *grid_index_1=int(((index1-this->m_n_states_neut*this->m_gsize_x)-*state_index_1*(this->m_gsize_x*this->m_n_states_cont))%this->m_gsize_x);
   }
   if(index2 == 0)
   {
      *state_index_2=0;
      *grid_index_2=0;
      *state_index_cont_2=-1;
   }
   else if(index2 < this->m_n_states_neut*this->m_gsize_x)
   {
      *state_index_2=int(index2/this->m_gsize_x);
      *grid_index_2=int(index2%this->m_gsize_x);
      *state_index_cont_2=-1;
   }
   else
   {
      *state_index_2=int((index2-this->m_n_states_neut*this->m_gsize_x)/(this->m_gsize_x*this->m_n_states_cont));
      *state_index_cont_2=int(((index2-this->m_n_states_neut*this->m_gsize_x)-*state_index_2*(this->m_gsize_x*this->m_n_states_cont))/this->m_gsize_x);
      *grid_index_2=int(((index2-this->m_n_states_neut*this->m_gsize_x)-*state_index_2*(this->m_gsize_x*this->m_n_states_cont))%this->m_gsize_x);
   }
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::xmin()
{
   return this->m_xmin;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::xmax()
{
   return this->m_xmax;
}
//##########################################################################
//
//##########################################################################
//END OF HAMILTON MATRIX OBJECT COMPUTATION
