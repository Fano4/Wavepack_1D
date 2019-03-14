#include "Computation.hpp"

//##########################################################################
//
//Runge Kutta routine is the ODE integrator. it is implemented to the sixth order, with an older implementation to fourth order. It takes both the wave function and the 
//Hamiltonian objects as arguments, along with the time-index for time-dependent Hamiltonian.
//
//##########################################################################
bool Runge_kutta(wavefunction *Psi0,hamilton_matrix *H, int time_index)
{

    wavefunction* k1=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k2=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k3=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k4=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k5=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k6=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k7=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());

    wavefunction* temp=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* temp2=new wavefunction(Psi0->gsize_x(),Psi0->tgsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());

    int state_index;
    int grid_index;
    int state_index_cont;
    std::complex<double> ctemp;

    //const double a21(.5),a32(.5),a43(1.),ddc1(1./6.),ddc2(1./3.),ddc3(1./3.),ddc4(1./6.);//RK coeff 4th order
    const double a21(1);
    const double a31(3./8.),a32(1./8.);
    const double a41(8./27.),a42(2./27.),a43(8./27.);
    const double  a51(3.*(3.*sqrt(21.)-7.)/392.),a52(-8.*(7.- sqrt(21.))/392.),a53(48.*(7.-sqrt(21.))/392.),a54(-3.*(21.-sqrt(21.))/392.);
    const double a61(-5.*(231.+51.*sqrt(21.))/1960.),a62(-40.*(7.+sqrt(21.))/1960.),a63(-320.*sqrt(21.)/1960.),a64(3.*(21.+121.*sqrt(21.))/1960.),a65(392.*(6.+sqrt(21.))/1960.);
    const double a71(15.*(22.+7*sqrt(21.))/180.),a72(120./180.),a73(40.*(7.*sqrt(21.)-5.)/180.),a74(-63.*(3.*sqrt(21.)-2.)/180.),a75(-14.*(49.+9.*sqrt(21))/180.),a76(70.*(7.-sqrt(21.))/180.);
    const double ddc1(9./180.),ddc2(0.),ddc3(64./180.),ddc4(0.),ddc5(49./180.),ddc6(49./180.),ddc7(9./180.);//SIXTH ORDER RK COEFF

    double vector[3];
    H->electric_field(time_index,vector);
    double efield_magnitude(sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]));
    bool cat=bool(efield_magnitude >= H->efield_thresh());


    // k1 = dPsi(t0) / dt | t = t0

    t_deriv(Psi0,H,k1,time_index);
    temp->set_wf(Psi0,cat);
    ctemp=std::complex<double>(a21*H->h(),0);
    temp->add_wf(&ctemp,k1,cat);

    // k2 = d[ a21*k1 + Psi(t0) ]/dt | t = t0 + h

    t_deriv(temp,H,k2,time_index+1);
    temp->set_wf(Psi0,cat);
    ctemp=std::complex<double>(a31*H->h(),0);
    temp->add_wf(&ctemp,k1,cat);
    ctemp=std::complex<double>(a32*H->h(),0);
    temp->add_wf(&ctemp,k2,cat);

    // k3 = d[ a31*k1 + a32*k2 + Psi(t0) ]/dt | t = t0 + .5*h

    t_deriv(temp,H,k3,time_index+.5);
    temp->set_wf(Psi0,cat);
    ctemp=std::complex<double>(a41*H->h(),0);
    temp->add_wf(&ctemp,k1,cat);
    ctemp=std::complex<double>(a42*H->h(),0);
    temp->add_wf(&ctemp,k2,cat);
    ctemp=std::complex<double>(a43*H->h(),0);
    temp->add_wf(&ctemp,k3,cat);

    // k4 = d[ a41*k1 + a42*k2 + a43*k3 + Psi(t0) ]/dt | t = t0 + (2/3)*h

    t_deriv(temp,H,k4,time_index+2./3.);
    temp->set_wf(Psi0,cat);
    ctemp=std::complex<double>(a51*H->h(),0);
    temp->add_wf(&ctemp,k1,cat);
    ctemp=std::complex<double>(a52*H->h(),0);
    temp->add_wf(&ctemp,k2,cat);
    ctemp=std::complex<double>(a53*H->h(),0);
    temp->add_wf(&ctemp,k3,cat);
    ctemp=std::complex<double>(a54*H->h(),0);
    temp->add_wf(&ctemp,k4,cat);

    // k5 = d[ a51*k1 + a52*k2 + a53*k3 + a54*k4 + Psi(t0) ]/dt | t = t0 + ( 7-sqrt(21)/14 )*h

    t_deriv(temp,H,k5,time_index+(7.-sqrt(21))/14);
    temp->set_wf(Psi0,cat);
    ctemp=std::complex<double>(a61*H->h(),0);
    temp->add_wf(&ctemp,k1,cat);
    ctemp=std::complex<double>(a62*H->h(),0);
    temp->add_wf(&ctemp,k2,cat);
    ctemp=std::complex<double>(a63*H->h(),0);
    temp->add_wf(&ctemp,k3,cat);
    ctemp=std::complex<double>(a64*H->h(),0);
    temp->add_wf(&ctemp,k4,cat);
    ctemp=std::complex<double>(a65*H->h(),0);
    temp->add_wf(&ctemp,k5,cat);

    // k6 = d[ a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5 + Psi(t0) ]/dt | t = t0 + ( 7+sqrt(21)/14 )*h

    t_deriv(temp,H,k6,time_index+(7.+sqrt(21))/14);
    temp->set_wf(Psi0,cat);
    ctemp=std::complex<double>(a71*H->h(),0);
    temp->add_wf(&ctemp,k1,cat);
    ctemp=std::complex<double>(a72*H->h(),0);
    temp->add_wf(&ctemp,k2,cat);
    ctemp=std::complex<double>(a73*H->h(),0);
    temp->add_wf(&ctemp,k3,cat);
    ctemp=std::complex<double>(a74*H->h(),0);
    temp->add_wf(&ctemp,k4,cat);
    ctemp=std::complex<double>(a75*H->h(),0);
    temp->add_wf(&ctemp,k5,cat);
    ctemp=std::complex<double>(a76*H->h(),0);
    temp->add_wf(&ctemp,k6,cat);

    // k7 = d[ a71*k1 + a72*k2 + a73*k3 + a74*k4 + a75*k5 + a76*k6 + Psi(t0) ]/dt | t = t0 + h

    t_deriv(temp,H,k7,time_index+1);
    temp->set_wf(Psi0,cat);

    // Psi(t) = ddc1 * h * k1 + ddc2 * h * k2 + ddc3 * h * k3 + ddc4 * h * k4 + ddc5 * h * k5 + ddc6 * h * k6 + ddc7 * h * k7

    ctemp=std::complex<double>(ddc1*H->h(),0);
    temp->add_wf(&ctemp,k1,cat);
    ctemp=std::complex<double>(ddc2*H->h(),0);
    temp->add_wf(&ctemp,k2,cat);
    ctemp=std::complex<double>(ddc3*H->h(),0);
    temp->add_wf(&ctemp,k3,cat);
    ctemp=std::complex<double>(ddc4*H->h(),0);
    temp->add_wf(&ctemp,k4,cat);
    ctemp=std::complex<double>(ddc5*H->h(),0);
    temp->add_wf(&ctemp,k5,cat);
    ctemp=std::complex<double>(ddc6*H->h(),0);
    temp->add_wf(&ctemp,k6,cat);
    ctemp=std::complex<double>(ddc7*H->h(),0);
    temp->add_wf(&ctemp,k7,cat);
    Psi0->set_wf(temp,cat);

    delete k1;
    delete k2;
    delete k3;
    delete k4;
    delete k5;
    delete k6;
    delete k7;
    delete temp;
    delete temp2;

    return 0;
}
//##########################################################################
//
//##########################################################################
bool adam_bashforth_moulton(wavefunction *dPsim4,wavefunction *dPsim3, wavefunction *dPsim2,wavefunction *dPsim1,hamilton_matrix *H,wavefunction *Psim1,int time_index)
{
    wavefunction *temp=new wavefunction(Psim1->gsize_x(),Psim1->tgsize_x(),Psim1->n_states_neut(),Psim1->n_states_cat(),Psim1->n_states_cont());
    wavefunction *dtemp=new wavefunction(Psim1->gsize_x(),Psim1->tgsize_x(),Psim1->n_states_neut(),Psim1->n_states_cat(),Psim1->n_states_cont());

    double vector[3];
    H->electric_field(time_index,vector);
    double efield_magnitude(sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]));
    bool cat=bool(efield_magnitude >= H->efield_thresh());
    std::complex<double> ctemp(std::complex<double>(0,0));

    const double am1(55./24.), am2(-59./24.), am3(37./24.),am4(-3./8.);
    const double b0(251./720.), bm1(646./720.), bm2(-264./720.), bm3(106./720.),bm4(-19./720.);

    temp->set_wf(Psim1,cat);
    ctemp=std::complex<double>(am1*H->h(),0);
    temp->add_wf(&ctemp,dPsim1,cat);
    ctemp=std::complex<double>(am2*H->h(),0);
    temp->add_wf(&ctemp,dPsim2,cat);
    ctemp=std::complex<double>(am3*H->h(),0);
    temp->add_wf(&ctemp,dPsim3,cat);
    ctemp=std::complex<double>(am4*H->h(),0);
    temp->add_wf(&ctemp,dPsim4,cat);

/*
    t_deriv(temp,H,dtemp,time_index);

    temp->set_wf(Psim1,cat);
    ctemp=std::complex<double>(b0*H->h(),0);
    temp->add_wf(&ctemp,dtemp,cat);
    ctemp=std::complex<double>(bm1*H->h(),0);
    temp->add_wf(&ctemp,dPsim1,cat);
    ctemp=std::complex<double>(bm2*H->h(),0);
    temp->add_wf(&ctemp,dPsim2,cat);
    ctemp=std::complex<double>(bm3*H->h(),0);
    temp->add_wf(&ctemp,dPsim3,cat);
    ctemp=std::complex<double>(bm4*H->h(),0);
    temp->add_wf(&ctemp,dPsim4,cat);
*/

    dPsim4->set_wf(dPsim3,cat);
    dPsim3->set_wf(dPsim2,cat);
    dPsim2->set_wf(dPsim1,cat);
    t_deriv(temp,H,dPsim1,time_index);

    Psim1->set_wf(temp,cat);


    delete temp;
    delete dtemp;
    return 0;
}
//##########################################################################
//
//The t_deriv routine computes the time-derivative of the wavefunction Psi at the instant given by time_index using the TDSE, from the Hamiltonian H. 
//The derivative is written in the wavefunction dPsi.
//This routine is usually called either by the Runge-Kutta routine or for energy computation.
//
//##########################################################################

bool t_deriv(wavefunction *Psi,hamilton_matrix *H,wavefunction *dPsi,double time_index)
{
   using namespace std;

   double* vector=new double[3];
   H->electric_field(time_index,vector);
   double efield_magnitude(sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]));
   bool condition(0);
   int state_index;
   int grid_index;
   int state_index_cont;
   int state_index2;
   int grid_index2;
   int state_index_cont2;
   double rsum;
   double imsum;
   int j(0);

#pragma omp parallel for reduction(+:rsum,imsum) private(state_index,grid_index,state_index_cont,state_index2,grid_index2,state_index_cont2,j)
   for(int i=0;i<(Psi->n_states_neut())*(Psi->tgsize_x());i++)//read every line in the Neutral part of the matrix
   {
      std::complex<double> ctemp;
      j=0;
      rsum=0;
      imsum=0;
      H->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);
      /*
       * N-N block of the hamilton matrix
       */
      for(int s=0;s<Psi->n_states_neut();s++)
      {
            j=(grid_index-2)*bool(grid_index-2 >= 0)+s*Psi->tgsize_x();//initial grid index of column in diagonal block
            while( (j%Psi->tgsize_x()) <= grid_index+2 && (j%Psi->tgsize_x()) < Psi->tgsize_x()-1)
            {
               H->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);

               ctemp=(H->hamilt_element(time_index,i,j))*(Psi->show_neut_psi(grid_index2,state_index2));
               rsum+=real(ctemp);//energy
               imsum+=imag(ctemp);//energy
               j++;
            }
      }
      /*
       * End of  N-N block of the hamilton matrix
       */
      if(efield_magnitude >= H->efield_thresh())
      {
         /*
          *  N-C block of the hamilton matrix
          */
         for(int s=0;s<Psi->n_states_cat()*Psi->n_states_cont();s++)
         {
            j=Psi->n_states_neut()*Psi->tgsize_x()+s*Psi->tgsize_x()+grid_index;
            H->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
            ctemp=(H->hamilt_element(time_index,i,j))*Psi->show_cat_psi(grid_index,state_index2,state_index_cont2);
            rsum+=real(ctemp);
            imsum+=imag(ctemp);
         }
         /*
          * End of N-C block of the hamilton matrix
          */
      }
      std::complex<double> sum(0,0);
      sum=std::complex<double>(rsum,imsum);
      dPsi->set_neut_psi(state_index,grid_index,std::complex<double>(0,-1)*sum);
   }
   if(efield_magnitude >= H->efield_thresh())
   {
#pragma omp parallel for reduction(+:rsum,imsum) private(state_index,grid_index,state_index_cont,state_index2,grid_index2,state_index_cont2,j)
      for(int i=(Psi->n_states_neut())*(Psi->tgsize_x());i<(Psi->n_states_neut())*(Psi->tgsize_x())+(Psi->n_states_cat())*(Psi->n_states_cont())*(Psi->tgsize_x());i++)//read every line in the Cation part of the matrix
      {
         rsum=0;
         imsum=0;
         j=0;
         std::complex<double> ctemp;

         H->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);
         /*
          * C-N block of the hamilton matrix
          */
         for(int s=0;s<Psi->n_states_neut();s++)
         {
            j=s*Psi->tgsize_x()+grid_index;
            ctemp=(H->hamilt_element(time_index,i,j))*Psi->show_neut_psi(grid_index,s);
            rsum+=real(ctemp);//state_index2));
            imsum+=imag(ctemp);//state_index2));
         }
         /*
          *  End of C-N block of the hamilton matrix
          */

         /*
          * C-C block of the hamilton matrix
          */
         for(int s=0;s<(Psi->n_states_cat());s++)
         {
            j=(Psi->n_states_neut())*(Psi->tgsize_x())+(i%Psi->tgsize_x()-2)*bool(i%Psi->tgsize_x()-2 >= 0)+s*Psi->n_states_cont()*Psi->tgsize_x()+state_index_cont*Psi->tgsize_x();//initial grid index of column in diagonal block

            while( (j%Psi->tgsize_x()) <= grid_index+2 && (j%Psi->tgsize_x()) < Psi->tgsize_x()-1)
            {
               H->show_indexes(i,j,&state_index,&grid_index,&state_index_cont,&state_index2,&grid_index2,&state_index_cont2);
               ctemp=(H->hamilt_element(time_index,i,j))*Psi->show_cat_psi(grid_index2,state_index2,state_index_cont2);
               rsum+=real(ctemp);//energy
               imsum+=imag(ctemp);//energy
               j++;

            }
         }
         /*
          *  End of C-C block of the hamilton matrix
          */
         std::complex<double> sum(0,0);
         sum=std::complex<double>(rsum,imsum);
         dPsi->set_cat_psi(state_index,state_index_cont,grid_index,std::complex<double>(0,-1)*sum);
      }
   }
   delete [] vector;
    return 0;
}
//##########################################################################
//
//The propagate routine is a routine called by the main function to propagate the wave function. It is this routine that calls Runge-Kutta and makes the efield and 
//pot vec evolve. it also checks that the thresholds for taking the effect of efield or pot vec of ionization are crossed or not.
//When it is needed, this routine calls for the reco;putation of the PICE, under translation by the pot_vec.
//
//##########################################################################
void propagate(wavefunction *Psi, hamilton_matrix *H,int* time_index,int num_of_loop)
{
   double *pot_vec=new double [3];
   double vector[3];
   int temp_time_index(*time_index);
   double efield_magnitude(0);
   bool analytic=1;


   for(int i=0;i!=num_of_loop;i++)
   {
      H->electric_field(*time_index,vector);
      efield_magnitude=(sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]));

      if(efield_magnitudei >= H->efield_thresh())
      {
         analytic=0;
      }
      temp_time_index=temp_time_index+1;
   }

   if(!analytic)
   {
      for(int i=0;i!=num_of_loop;i++)
      {
         H->potential_vector(*time_index,pot_vec);
         H->electric_field(*time_index,vector);
         efield_magnitude=(sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]));
         Runge_kutta(Psi,H,*time_index);
         *time_index=*time_index+1;
      }
   }
   else
   {
      Psi->projection_eigenstates(H,1);
      Psi->analytic_propagation(H,num_of_loop)
      Psi->projection_eigenstates(H,-1);
   }
   delete [] pot_vec;
}

