

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
//   std::cout<<"kmax = "<<this->k_modulus[this->m_n_k-1]<<" kmin = "<<this->k_modulus[0]<<" nk = "<<this->m_n_k<<" => dk = "<<dk<<std::endl;
//   exit(EXIT_SUCCESS);

   this->show_indexes(i,j,&state_index_1,&grid_index_1,&state_index_cont_1,&state_index_2,&grid_index_2,&state_index_cont_2);


   this-> electric_field(time_index,elec_field);
   this-> potential_vector(time_index,pot_vector);

   if(state_index_cont_1 == -1 && state_index_cont_2 == -1)//N-N section
   {
         
      if(state_index_1==state_index_2)//diagonal element in the electronic state basis
      {
         if(grid_index_1==grid_index_2)//diagonal in the grid basis
         {
            return std::complex<double>(
                  this->kinetic_energy[this->m_tgsize_x*grid_index_1+grid_index_2]
                  +this->m_pot_neut[state_index_1*(this->m_tgsize_x)+grid_index_1]
                  -this->m_dmx_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[0]
                  -this->m_dmy_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[1]
                  -this->m_dmz_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[2]
                  ,0); //diagonal term of kinetic energy + potential energy + permanent dipole interaction
         }
         else // non-diagonal term of kinetic energy
         {
            return std::complex<double>(this->kinetic_energy[this->m_tgsize_x*grid_index_1+grid_index_2],0);
         }
      }
      else //coupling elements between electronic states of the neutral
      {
         if(grid_index_1==grid_index_2)//diagonal in the grid basis
         {
            //return 0;
            if(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]!=0 && state_index_1 < state_index_2)
                  return std::complex<double>(
                     -this->m_dmx_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[0]
                     -this->m_dmy_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[1]
                     -this->m_dmz_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[2]

            //         -0.5*(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]
            //         -(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][(grid_index_1-2)*bool(grid_index_1-2>=0)]
            //         -8*this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][(grid_index_1-1)*bool(grid_index_1-1>=0)]
            //         +8*this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][(grid_index_1+1)*bool(grid_index_1+1<this->m_tgsize_x)+this->m_tgsize_x*bool(!(grid_index_1+1<this->m_tgsize_x))]
            //         -this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][(grid_index_1+2)*bool(grid_index_1+2<this->m_tgsize_x)+this->m_tgsize_x*bool(!(grid_index_1+2<this->m_tgsize_x))]
            //         )/(12*(this->m_xmax-this->m_xmin)/this->m_gsize_x)
            //         )/((this->m_mass)*(fabs(this->m_pot_neut[state_index_1*this->m_tgsize_x+grid_index_1]-this->m_pot_neut[state_index_2*this->m_tgsize_x+grid_index_1])))
                    ,0); //dipole coupling
            else
               return std::complex<double>(
                  -this->m_dmx_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[0]
                  -this->m_dmy_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[1]
                  -this->m_dmz_neut[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*elec_field[2]
                  ,0);

         }
         else //NACME
         {
        /*    if(fabs(this->m_pot_neut[state_index_1*this->m_tgsize_x+grid_index_1]-this->m_pot_neut[state_index_2*this->m_tgsize_x+grid_index_1]) <= 1e-6)
            {
               //std::cout<<"DIVISION BY ZERO IN NAC !!! "<<state_index_1<<"="<<state_index_2<<" AT "<<grid_index_1<<std::endl;
               //std::cout<<"Energies of PES: "<<this->m_pot_neut[state_index_1*this->m_tgsize_x+grid_index_1]<<"  "<<this->m_pot_neut[state_index_2*this->m_tgsize_x+grid_index_1]<<std::endl;
               return std::complex<double>(0,0);
            }*/

//               std::cout<<"NACME : "<<this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]<<std::endl;
            if(state_index_1<state_index_2)
            {
               if(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]!=0)
                  return std::complex<double>(
                     -(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]
                     *this->derivative_matrix[grid_index_1*this->m_tgsize_x+grid_index_2]
                     )/((this->m_mass))//*(fabs(this->m_pot_neut[state_index_1*this->m_tgsize_x+grid_index_1]-this->m_pot_neut[state_index_2*this->m_tgsize_x+grid_index_1])))
                     ,0);//Non adiabatic coupling
               else
                  return std::complex<double>(0,0);
            }
            else
            {
               if(this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_2]!=0)
                  return std::complex<double>(
                     (this->m_NAC[state_index_2*this->m_n_states_neut+state_index_1][grid_index_2]
                     *this->derivative_matrix[grid_index_1*this->m_tgsize_x+grid_index_2]
                     )/((this->m_mass))//*(fabs(this->m_pot_neut[state_index_1*this->m_tgsize_x+grid_index_2]-this->m_pot_neut[state_index_2*this->m_tgsize_x+grid_index_2])))
                     ,0);//Non adiabatic coupling
               else
                  return std::complex<double>(0,0);
            }
         }
      }
   }


   else if( (state_index_cont_1 ==-1 || state_index_cont_2==-1) && !(state_index_cont_2==-1 && state_index_cont_1 ==-1) )//N-C or C-N section => PICE
   {
      if(state_index_cont_1 ==-1 && grid_index_1==grid_index_2)//N-C
      {
      //   std::cout<<"In hamiltonian, asking for element "<<state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+state_index_cont_2<<std::endl;
         return std::complex<double>((
                  -std::conj(this->m_PICE_sto_x[this->pice_time_mapping[int(floor(time_index))]][state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+state_index_cont_2][grid_index_1])*elec_field[0]
                  -std::conj(this->m_PICE_sto_y[this->pice_time_mapping[int(floor(time_index))]][state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+state_index_cont_2][grid_index_1])*elec_field[1]
                  -std::conj(this->m_PICE_sto_z[this->pice_time_mapping[int(floor(time_index))]][state_index_1*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_2*(this->m_n_states_cont)+state_index_cont_2][grid_index_1])*elec_field[2])
               *sqrt(this->dk(state_index_cont_2)));//PICE dipole interaction
      }
      else if(grid_index_1==grid_index_2) // C-N
      {
         return std::complex<double>((
                  -(this->m_PICE_sto_x[this->pice_time_mapping[int(floor(time_index))]][state_index_2*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_1*(this->m_n_states_cont)+state_index_cont_1][grid_index_1])*elec_field[0]
                  -(this->m_PICE_sto_y[this->pice_time_mapping[int(floor(time_index))]][state_index_2*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_1*(this->m_n_states_cont)+state_index_cont_1][grid_index_1])*elec_field[1]
                  -(this->m_PICE_sto_z[this->pice_time_mapping[int(floor(time_index))]][state_index_2*(this->m_n_states_cat)*(this->m_n_states_cont)+state_index_1*(this->m_n_states_cont)+state_index_cont_1][grid_index_1])*elec_field[2])
               *sqrt(this->dk(state_index_cont_1)));//PICE dipole interaction
      }
      else
         return 0;//non-diagonal in grid index
   }
   else//C-C section
   {
      if(state_index_cont_1==state_index_cont_2 && state_index_cont_1!=-1)
      {
         if(state_index_1==state_index_2)//diagonal element in the electronic state basis
         {
            if(grid_index_1==grid_index_2)//diagonal in the grid basis
            {
               return (std::complex<double>(
                        0.5*(
                           pow(this->k_modulus[(state_index_cont_1-state_index_cont_1%(this->m_n_states_cont/this->m_n_k))/(this->m_n_states_cont/this->m_n_k)]*sin(this->k_orientation[0][state_index_cont_1%(this->m_n_states_cont/this->m_n_k)])*cos(this->k_orientation[1][state_index_cont_1%(this->m_n_states_cont/this->m_n_k)])+pot_vector[0],2)
                           +pow(k_modulus[(state_index_cont_1-state_index_cont_1%(this->m_n_states_cont/this->m_n_k))/(this->m_n_states_cont/this->m_n_k)]*sin(this->k_orientation[0][state_index_cont_1%(this->m_n_states_cont/this->m_n_k)])*sin(this->k_orientation[1][state_index_cont_1%(this->m_n_states_cont/this->m_n_k)])+pot_vector[1],2)
                           +pow(this->k_modulus[(state_index_cont_1-state_index_cont_1%(this->m_n_states_cont/this->m_n_k))/(this->m_n_states_cont/this->m_n_k)]*cos(this->k_orientation[0][state_index_cont_1%(this->m_n_states_cont/this->m_n_k)])+pot_vector[2],2))
                        +this->kinetic_energy[grid_index_1*this->m_tgsize_x+grid_index_2]
                        +this->m_pot_cat[state_index_1*this->m_tgsize_x+grid_index_1]
                        -this->m_dmx_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[0]
                        -this->m_dmy_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[1]
                        -this->m_dmz_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[2],0)); //diagonal term of kinetic energy + potential energy + permanent dipole interaction
            }
            else // non-diagonal term of kinetic energy
            {
               return std::complex<double>(this->kinetic_energy[grid_index_1*this->m_tgsize_x+grid_index_2],0);//*this->dk(state_index_cont_2);
            }
         }
         else //coupling elements between electronic states of the cation
         {
            if(grid_index_1==grid_index_2)//diagonal in the grid basis
            {
               return std::complex<double>(
                     -this->m_dmx_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[0]
                     -this->m_dmy_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[1]
                     -this->m_dmz_cat[state_index_1*this->m_n_states_cat+state_index_2][grid_index_1]*elec_field[2]);
//                  *this->dk(state_index_cont_2); //dipole coupling
            }
            else
            {
               if(state_index_1 < state_index_2)//upper triangle of the matrix
               {
                  return 0;
                     //return std::complex<double>(0,m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*derivative_matrix[grid_index_1*this->m_tgsize_x+grid_index_2]/(this->m_mass*fabs(this->m_pot_neut[state_index_1*this->m_tgsize_x+grid_index_1]-this->m_pot_neut[state_index_2*this->m_tgsize_x+grid_index_2])));//*this->dk(state_index_cont_2));//Non adiabatic coupling
               }
               else //lower triangle of the matrix
               {
                  return 0;
                  //return std::complex<double>(0,m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*derivative_matrix[grid_index_1*this->m_tgsize_x+grid_index_2]/(this->m_mass*fabs(this->m_pot_neut[state_index_1*this->m_tgsize_x+grid_index_1]-this->m_pot_neut[state_index_2*this->m_tgsize_x+grid_index_2])));//*this->dk(state_index_cont_2));//Non adiabatic coupling
               }
            }
         }
      }
   }
   std::cout<<"FATAL ERROR !! UNASSIGNED MATRIX ELEMENT IN HAMILTON_MATRIX::HAMILT_ELEMENT"<<std::endl; 
   return -1;
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
   double rec(0);
   double imec(0);
   double e(0);
   double cnorm(0);
   double tempnorm(0);
   std::complex<double> *dipvec=new std::complex<double> [3];
   int index_1(0);
   int index_2(0);
   int state_index;
   int grid_index;
   int state_index_cont;
   int state_index2;
   int grid_index2;
   int state_index_cont2;
   int j(0);
   std::complex<double> imunit (0,1);
   bool cat=1;

   if(Psi->n_states_cat() == 0)
      cat=0;

   wavefunction* dPsi=new wavefunction(Psi->gsize_x(),Psi->tgsize_x(),Psi->n_states_neut(),Psi->n_states_cat(),Psi->n_states_cont());

   std::cout<<"Energy computation...";
   t_deriv(Psi,this,dPsi,time_index);
   dPsi->add_wf(&imunit,NULL,cat);
   ec=dPsi->dot_prod(Psi,this);
   std::cout<<"Done !"<<std::endl;

   //Psi->norm(this);

   delete dPsi;
   delete [] dipvec;
   return real(ec);


//This block reads the elements of the sparse hamiltonian matrix in the basis of electronic states.

//std::cout<<"N block... "<<(Psi->n_states_neut())*(Psi->gsize_x())<<" elements."<<std::endl;
/*#pragma omp parallel for reduction(+:rec,imec,tempnorm) private(state_index,grid_index,state_index_cont,state_index2,grid_index2,state_index_cont2,j)
   for(int i=0;i!=(Psi->n_states_neut())*(Psi->gsize_x());i++)//read every line in the Neutral part of the matrix
   {
      j=0;

      this->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);
      /////////NORM COMPUTATION
      tempnorm+=(conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index))).real();
      /////////NORM COMPUTATION

      *
       * N-N block of the hamilton matrix
       *
      for(int s=0;s!=Psi->n_states_neut();s++)
      {
         if(s==state_index)//Diagonal block in electronic state => pentadiagonal in grid index
         {
            j=(i%Psi->gsize_x()-2)*bool(i%Psi->gsize_x()-2 >= 0)+s*Psi->gsize_x();//initial grid index of column in diagonal block
            while( (j%Psi->gsize_x()) <= i%Psi->gsize_x()+2)
            {
               this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);

               rec+=real((this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index2,state_index2)));//energy
               imec+=imag((this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index2,state_index2)));//energy

               ///////////DIPOLE COMPUTATION
               if(grid_index == grid_index2)
               {
                  dipvec[0]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmx_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
                  dipvec[1]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmy_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
                  dipvec[2]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmz_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
               }
               ///////////DIPOLE COMPUTATION
               
               if(!(grid_index2 == Psi->gsize_x()-1))//if you do not reach the edge of the grid, continue
                  j++;
               else //end of the while loop
               {
                  j++;
                  break;
               }
            }
         }
         else//Non-diagonal block in electronic state => diagonal in grid index
         {
            j=grid_index+s*Psi->gsize_x();//grid index of column in diagonal block
            this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
            rec+=real((this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index2,state_index2)));//energy
            imec+=imag((this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index2,state_index2)));//energy

            ///////////DIPOLE COMPUTATION
            dipvec[0]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmx_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
            dipvec[1]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmy_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
            dipvec[2]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_neut_psi(grid_index,state_index2))*(this->m_dmz_neut[state_index*this->m_n_states_neut+state_index2][grid_index]);
            ///////////DIPOLE COMPUTATION
         }
      }
      *
       * End of  N-N block of the hamilton matrix
       *
      *
       *  N-C block of the hamilton matrix
       *
      for(int s=0;s!=this->m_n_states_cat*this->m_n_states_cont;s++)
      {
         j=this->m_n_states_neut*this->m_gsize_x+s*this->m_gsize_x+grid_index;
         this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
         rec+=real((this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*Psi->show_cat_psi(grid_index,state_index2,state_index_cont2));
         imec+=imag((this->hamilt_element(time_index,i,j))*conj(Psi->show_neut_psi(grid_index,state_index))*Psi->show_cat_psi(grid_index,state_index2,state_index_cont2));
      }
      *
       *  N-C block of the hamilton matrix
       *
   }
   ec+=std::complex<double>(rec,imec);
   cnorm+=tempnorm;
   tempnorm=0;
   rec=0;
   imec=0;
//std::cout<<"C block... "<<(Psi->n_states_cat())*(Psi->n_states_cont())*(Psi->gsize_x())<<" elements."<<std::endl;

#pragma omp parallel for reduction(+:rec,imec,tempnorm) private(state_index,grid_index,state_index_cont,state_index2,grid_index2,state_index_cont2,j)
   for(int i=(Psi->n_states_neut())*(Psi->gsize_x());i!=(Psi->n_states_neut())*(Psi->gsize_x())+(Psi->n_states_cat())*(Psi->n_states_cont())*(Psi->gsize_x());i++)//read every line in the Cation part of the matrix
   {
      j=0;

      this->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);
      /////////NORM COMPUTATION
      tempnorm+=(conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index,state_index_cont))).real();
      /////////NORM COMPUTATION
      *
       * C-N block of the hamilton matrix
       *
      for(int s=0;s!=this->m_n_states_neut;s++)
      {
         j=s*this->m_gsize_x+grid_index;
         this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
         rec+=real((this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*Psi->show_neut_psi(grid_index,state_index2));
         imec+=imag((this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*Psi->show_neut_psi(grid_index,state_index2));
      }
      *
       *  End of C-N block of the hamilton matrix
       *
      *
       * C-C block of the hamilton matrix
       *
      for(int s=0;s!=(Psi->n_states_cat());s++)
      {
            if(s==state_index)//Diagonal block in electronic state => pentadiagonal in grid index
            {
               j=(Psi->n_states_neut())*(Psi->gsize_x())+(i%Psi->gsize_x()-2)*bool(i%Psi->gsize_x()-2 >= 0)+s*Psi->n_states_cont()*Psi->gsize_x();//initial grid index of column in diagonal block
               while( (j%Psi->gsize_x()) <= i%Psi->gsize_x()+2)
               {
                  this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);

                  rec+=real((this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index2,state_index2,state_index_cont2)));//energy
                  imec+=imag((this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index2,state_index2,state_index_cont2)));//energy

                  ///////////DIPOLE COMPUTATION
                  if(grid_index == grid_index2)
                  {
                     dipvec[0]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmx_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
                     dipvec[1]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmy_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
                     dipvec[2]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmz_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
                  }
                  ///////////DIPOLE COMPUTATION
               
                  if(!(grid_index2 == Psi->gsize_x()-1))//if you do not reach the edge of the grid, continue
                     j++;
                  else //end of the loop
                  {
                     j++;
                     break;
                  }
               }
            }
            else//Non-diagonal block in electronic state => diagonal in grid index
            {
               j=(Psi->n_states_neut())*(Psi->gsize_x())+i%Psi->gsize_x()+s*Psi->gsize_x()*Psi->n_states_cont()+state_index_cont*Psi->gsize_x();//grid index of column in diagonal block
               this->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
               rec+=real((this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index2,state_index2,state_index_cont2)));//energy
               imec+=imag((this->hamilt_element(time_index,i,j))*conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index2,state_index2,state_index_cont2)));//energy

               ///////////DIPOLE COMPUTATION
               dipvec[0]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmx_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
               dipvec[1]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmy_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
               dipvec[2]+=conj(Psi->show_cat_psi(grid_index,state_index,state_index_cont))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_dmz_cat[state_index*this->m_n_states_cat+state_index2][grid_index]);
               ///////////DIPOLE COMPUTATION
            }
         }
      *
       *  End of C-C block of the hamilton matrix
       *
   }
   ec+=std::complex<double>(rec,imec);
   cnorm+=tempnorm;
   tempnorm=0;
*      for(int i=0;i!=(Psi->n_states_neut())*(Psi->gsize_x());i++)//read every line in the Neutral part of the matrix
      {
         int j(0);
         
         j=(i%Psi->gsize_x()-2)*bool(i%Psi->gsize_x()-2 >= 0);//grid index of the column 
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
            //std::cout<<"probe "<<i<<","<<s<<std::endl;//DEBUGGING
            j += Psi->gsize_x();
               if(i%Psi->gsize_x() == j%Psi->gsize_x())
               {
                  dipvec[0]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_PICE_x[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
                  dipvec[1]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_PICE_y[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
                  dipvec[2]+=conj(Psi->show_neut_psi(grid_index,state_index))*(Psi->show_cat_psi(grid_index,state_index2,state_index_cont2))*(this->m_PICE_z[i*this->m_n_states_cat*this->m_n_states_cont+j][grid_index]);
               }
         }
         //std::cout<<"probe "<<i<<std::endl;//DEBUGGING

      }

      //std::cout<<"probe_neutral"<<std::endl;//DEBUGGING
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
*/

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
   //double* rdipvec = new double [3];

   //rdipvec[0]=dipvec[0].real();
   //rdipvec[1]=dipvec[1].real();
   //rdipvec[2]=dipvec[2].real();

   //delete [] dipvec;

   //e=ec.real();
   //Psi->set_norm(cnorm);
   //Psi->set_dipole(rdipvec);
   //delete [] rdipvec;

   //return e;
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
   *state_index_1 = 0;
   *grid_index_1 = 0;
   *state_index_cont_1 = -1;
   *state_index_2 = 0;
   *grid_index_2 = 0;
   *state_index_cont_2 = -1;

   if(index1==0)
   {
      *state_index_1 = 0;
      *grid_index_1 = 0;
      *state_index_cont_1 = -1;
   }
   else if(index1 < (this->m_n_states_neut)*(this->m_tgsize_x))
   {
      *grid_index_1 = int(index1%this->m_tgsize_x);
      *state_index_1 = int((index1-*grid_index_1)/this->m_tgsize_x);
      *state_index_cont_1 = -1;
   }
   else
   {
      *state_index_1=int((index1-this->m_n_states_neut*this->m_tgsize_x)/(this->m_tgsize_x*this->m_n_states_cont));
      *state_index_cont_1=int(((index1-this->m_n_states_neut*this->m_tgsize_x)-*state_index_1*(this->m_tgsize_x*this->m_n_states_cont))/this->m_tgsize_x);
      *grid_index_1=int(((index1-this->m_n_states_neut*this->m_tgsize_x)-*state_index_1*(this->m_tgsize_x*this->m_n_states_cont))%this->m_tgsize_x);
   }
   if(index2 == 0)
   {
      *state_index_2=0;
      *grid_index_2=0;
      *state_index_cont_2=-1;
   }
   else if(index2 < (this->m_n_states_neut)*(this->m_tgsize_x))
   {
      *grid_index_2 = int(index2%this->m_tgsize_x);
      *state_index_2 = int((index2-*grid_index_2)/this->m_tgsize_x);
      *state_index_cont_2=-1;
   }
   else
   {
      *state_index_2=int((index2-this->m_n_states_neut*this->m_tgsize_x)/(this->m_tgsize_x*this->m_n_states_cont));
      *state_index_cont_2=int(((index2-this->m_n_states_neut*this->m_tgsize_x)-*state_index_2*(this->m_tgsize_x*this->m_n_states_cont))/this->m_tgsize_x);
      *grid_index_2=int(((index2-this->m_n_states_neut*this->m_tgsize_x)-*state_index_2*(this->m_tgsize_x*this->m_n_states_cont))%this->m_tgsize_x);
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
int hamilton_matrix::n_k()
{
   return this->m_n_k;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::dk(int i)
{
   if(i < 0 || i >= this->m_n_states_cont)
   {
      std::cout<<"ERROR : TRYING TO ACEES INDEX OUT OF RANGE IN DK VECTOR"<<std::endl;
      exit(EXIT_FAILURE);
   }
   if(this->m_n_k == 0 )
      return 0;
   else
      return this->m_dk_vec[i];
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::dk_vec(double *array)
{
   if(this->m_n_k == 0 )
      return;
   else
      cblas_dcopy(this->m_n_states_cont,this->m_dk_vec,1,array,1);
   
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::show_nac(int state_index_1,int state_index_2,int grid_index_1,int grid_index_2,int state_index_cont_1)
{
   if(state_index_cont_1==0)
   {
      if(state_index_1<state_index_2)
         return this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_1]*this->derivative_matrix[grid_index_1*this->m_tgsize_x+grid_index_2];
      else
         return this->m_NAC[state_index_1*this->m_n_states_neut+state_index_2][grid_index_2]*this->derivative_matrix[grid_index_2*this->m_tgsize_x+grid_index_1];
   }
   else
      return 0;
}
//##########################################################################
//
//##########################################################################
/*double hamilton_matrix::grid_k_cube_spacing() 
{
   return (this->m_kxmax-this->m_kxmin)/this->m_n_kx;
}*/
//##########################################################################
//
//##########################################################################
double hamilton_matrix::pot_vec_thresh() const
{
   return this->m_pot_vec_thresh;
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_pot_vec_mod(double value)
{
   this->m_pot_vec_mod=value;
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::set_pot_vec_tm_mod(double value)
{
   this->m_pot_vec_tm_mod=value;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::pot_vec_mod()
{
   return this->m_pot_vec_mod;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::pot_vec_tm_mod()
{
   return this->m_pot_vec_tm_mod;
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::plot_integrated_cross_section(std::string file_address,int neut_state,int cat_state,int pos_index)
{

   using namespace std;
   std::complex<double> cs(0);
   int dgsize(this->m_tgsize_x-this->m_small_gsize_x);

   ofstream cs_out;

   cs_out.open(file_address.c_str());
   if(!cs_out.is_open())
   {
      std::cout<<"ERROR, IMPOSSIBLE TO OPEN "<<file_address.c_str()<<std::endl;
   }

   if(pos_index==-1)
   {
      for(int x=0;x!=this->m_tgsize_x;x++)
      {
         for(int k=0;k!=this->m_n_k;k++)
         {
            cs=0;
            for(int i=0;i!=this->m_n_angles;i++)
            {
//            std::cout<<neut_state<<","<<k<<","<<i<<std::endl<<this->m_PICE_x[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]<<","<<this->m_PICE_y[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]<<","<<this->m_PICE_z[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]<<std::endl;
               cs+=pow(abs(this->m_PICE_sto_z[0][neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]),2)*pow(this->k_modulus[k],2)*4*acos(-1)/this->m_n_angles;
          //  cs=this->m_PICE_z[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+1][x];
            }
            cs_out<<(this->m_xmin + (x-dgsize)*(this->m_xmax-this->m_xmin)/this->m_small_gsize_x)*0.529<<"   "<<pow(this->k_modulus[0]+k*(this->k_modulus[this->m_n_k-1]-this->k_modulus[0])/this->m_n_k,2)*27.211/2<<"   "<<real(cs)<<"   "<<imag(cs)<<std::endl;
         }cs_out<<std::endl;
      }
      cs_out.close();
   }
   else
   {
      int x=pos_index;
         for(int k=0;k!=this->m_n_k;k++)
         {
            cs=0;
            for(int i=0;i!=this->m_n_angles;i++)
            {
//            std::cout<<neut_state<<","<<k<<","<<i<<std::endl<<this->m_PICE_x[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]<<","<<this->m_PICE_y[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]<<","<<this->m_PICE_z[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]<<std::endl;
               cs+=pow(abs(this->m_PICE_sto_z[0][neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+i][x]),2)*pow(this->k_modulus[k],2)*4*acos(-1)/this->m_n_angles;
          //  cs=this->m_PICE_z[neut_state*this->m_n_states_cat*this->m_n_states_cont+cat_state*this->m_n_states_cont+k*this->m_n_angles+1][x];
            }
            cs_out<<(this->m_xmin + (x-dgsize)*(this->m_xmax-this->m_xmin)/this->m_small_gsize_x)*0.529<<"   "<<pow(this->k_modulus[0]+k*(this->k_modulus[this->m_n_k-1]-this->k_modulus[0])/this->m_n_k,2)*27.211/2<<"   "<<real(cs)<<"   "<<imag(cs)<<std::endl;
         }cs_out<<std::endl;
         cs_out.close();
   }
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::PI_rate(int time_index,double** ionization_rate,wavefunction* Psi)
{


   double* elec_field=new double[3];

   this->electric_field(time_index,elec_field);

   for(int i=0;i!=this->m_n_states_neut;i++)
   {
      for(int j=0;j!=this->m_n_states_cat;j++)
      {
         ionization_rate[i][j]=0;

         for(int k=0;k!=this->m_n_states_cont;k++)
         {
            for(int r=0;r!=this->m_tgsize_x;r++)
            {
               ionization_rate[i][j]-=2*std::imag(
                              (
                                 (this->m_PICE_sto_x[this->pice_time_mapping[int(floor(time_index))]][i*(this->m_n_states_cat)*(this->m_n_states_cont)+j*(this->m_n_states_cont)+k][r]*elec_field[0]
                                 +this->m_PICE_sto_y[this->pice_time_mapping[int(floor(time_index))]][i*(this->m_n_states_cat)*(this->m_n_states_cont)+j*(this->m_n_states_cont)+k][r]*elec_field[1]
                                 +this->m_PICE_sto_z[this->pice_time_mapping[int(floor(time_index))]][i*(this->m_n_states_cat)*(this->m_n_states_cont)+j*(this->m_n_states_cont)+k][r]*elec_field[2])
                              )
               *(Psi->show_cat_psi(r,j,k))*std::conj(Psi->show_neut_psi(r,i)));

            }
         }
      }
   }
      delete [] elec_field;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::mass()
{
   return this->m_mass;
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::show_derivative_matrix(int grid_index_1,int grid_index_2)
{
   return this->derivative_matrix[grid_index_1*this->m_tgsize_x+grid_index_2];
}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::change_basis_dipole(wavefunction **change_basis,double *dipole_new_basis)
{
   double *S=new double[this->m_n_states_neut*this->m_n_states_neut*this->m_tgsize_x*this->m_tgsize_x];
   double *A=new double[this->m_n_states_neut*this->m_n_states_neut*this->m_tgsize_x*this->m_tgsize_x];
   double *C=new double[this->m_n_states_neut*this->m_n_states_neut*this->m_tgsize_x*this->m_tgsize_x];

   std::cout<<"probe dipole"<<std::endl;
   for(int m=0;m!=this->m_n_states_neut;m++)
   {
      for(int n=0;n!=this->m_n_states_neut;n++)
      {
         for(int i=0;i!=this->m_tgsize_x;i++)
         {
            for(int j=0;j!=this->m_tgsize_x;j++)
            {
               S[m*this->m_n_states_neut*this->m_tgsize_x*this->m_tgsize_x+i*this->m_n_states_neut*this->m_tgsize_x+n*this->m_tgsize_x+j]=real(change_basis[m*this->m_tgsize_x+i]->show_neut_psi(j,n));
            }
               A[m*this->m_n_states_neut*this->m_tgsize_x*this->m_tgsize_x+i*this->m_n_states_neut*this->m_tgsize_x+n*this->m_tgsize_x+i]=this->show_dm_neut(m,n,i,2);
         }
      }
   }
   std::cout<<"probe dipole"<<std::endl;
   //A*D
   cblas_dgemm (CblasRowMajor, CblasTrans,CblasNoTrans, this->m_n_states_neut*this->m_tgsize_x, this->m_n_states_neut*this->m_tgsize_x, this->m_n_states_neut*this->m_tgsize_x,1, A, this->m_n_states_neut*this->m_tgsize_x, S, this->m_n_states_neut*this->m_tgsize_x, 0,C,this->m_n_states_neut*this->m_tgsize_x);
   std::cout<<"probe dipole"<<std::endl;
   //DT*(A*D)
   cblas_dgemm (CblasRowMajor, CblasNoTrans,CblasNoTrans, this->m_n_states_neut*this->m_tgsize_x, this->m_n_states_neut*this->m_tgsize_x, this->m_n_states_neut*this->m_tgsize_x,1, S, this->m_n_states_neut*this->m_tgsize_x, C, this->m_n_states_neut*this->m_tgsize_x, 0,dipole_new_basis,this->m_n_states_neut*this->m_tgsize_x);

   std::cout<<"probe dipole"<<std::endl;

   delete [] S;
   delete [] A;


}
//##########################################################################
//
//##########################################################################
void hamilton_matrix::eigenstate(int grid_index,int state_index,int state_index_cont,wavefunction* psi )
{
   if(state_index_cont==-1)
   {
      psi->set_wf(this->m_eigenstate[state_index*this->m_tgsize_x+grid_index],1);
   }
   else
   {
      psi->set_wf(this->m_eigenstate[this->m_n_states_neut*this->m_tgsize_x+state_index*this->m_n_states_cont*this->m_tgsize_x+state_index_cont*this->m_tgsize_x+grid_index],1);
   }
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::eigenvalue_neut(int state_index,int grid_index)
{
   return this->m_eigenvalue_neut[state_index*this->m_tgsize_x+grid_index];
}
//##########################################################################
//
//##########################################################################
double hamilton_matrix::eigenvalue_cat(int state_index,int state_index_cont,int grid_index)
{
   return this->m_eigenvalue_cat[state_index*this->m_tgsize_x*this->m_n_states_cont+state_index_cont*this->m_tgsize_x+grid_index];
}
//##########################################################################
//
//##########################################################################
//END OF HAMILTON MATRIX OBJECT COMPUTATION
