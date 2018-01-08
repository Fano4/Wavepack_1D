#include "Computation.hpp"

//##########################################################################
//
//##########################################################################
bool Runge_kutta_notdH(wavefunction *Psi0,hamilton_matrix *H, int time_index)
{
    wavefunction* k1=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k2=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k3=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k4=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k5=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k6=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k7=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());

    wavefunction* temp=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* temp2=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());

    wavefunction** dpsi_mat=new wavefunction *[(Psi0->n_states_neut())*(Psi0->gsize_x())+(Psi0->n_states_cat())*(Psi0->n_states_cont())*(Psi0->gsize_x())]; 
       for(int i=0;i!=(Psi0->n_states_neut())*(Psi0->gsize_x())+(Psi0->n_states_cat())*(Psi0->n_states_cont())*(Psi0->gsize_x());i++)
       {
          dpsi_mat[i]=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
       }


    int state_index;
    int grid_index;
    int state_index_cont;
    std::complex<double> ctemp;

    const double a21(1);
    const double a31(3./8.),a32(1./8.);
    const double a41(8./27.),a42(2./27.),a43(8./27.);
    const double  a51(3.*(3.*sqrt(21.)-7.)/392.),a52(-8.*(7.- sqrt(21.))/392.),a53(48.*(7.-sqrt(21.))/392.),a54(-3.*(21.-sqrt(21.))/392.);
    const double a61(-5.*(231.+51.*sqrt(21.))/1960.),a62(-40.*(7.+sqrt(21.))/1960.),a63(-320.*sqrt(21.)/1960.),a64(3.*(21.+121.*sqrt(21.))/1960.),a65(392.*(6.+sqrt(21.))/1960.);
    const double a71(15.*(22.+7*sqrt(21.))/180.),a72(120./180.),a73(40.*(7.*sqrt(21.)-5.)/180.),a74(-63.*(3.*sqrt(21.)-2.)/180.),a75(-14.*(49.+9.*sqrt(21))/180.),a76(70.*(7.-sqrt(21.))/180.);
    const double ddc1(9./180.),ddc2(0.),ddc3(64./180.),ddc4(0.),ddc5(49./180.),ddc6(49./180.),ddc7(9./180.);//SIXTH ORDER RK COEFF

    t_deriv_matrix(H,dpsi_mat,time_index);

    k1->matrix_prod(dpsi_mat,Psi0);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a21*H->h(),0);
    temp->add_wf(&ctemp,k1);

    k2->matrix_prod(dpsi_mat,temp);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a31*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a32*H->h(),0);
    temp->add_wf(&ctemp,k2);

    k3->matrix_prod(dpsi_mat,temp);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a41*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a42*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a43*H->h(),0);
    temp->add_wf(&ctemp,k3);

    k4->matrix_prod(dpsi_mat,temp);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a51*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a52*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a53*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(a54*H->h(),0);
    temp->add_wf(&ctemp,k4);

    k5->matrix_prod(dpsi_mat,temp);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a61*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a62*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a63*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(a64*H->h(),0);
    temp->add_wf(&ctemp,k4);
    ctemp=std::complex<double>(a65*H->h(),0);
    temp->add_wf(&ctemp,k5);

    k6->matrix_prod(dpsi_mat,temp);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a71*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a72*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a73*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(a74*H->h(),0);
    temp->add_wf(&ctemp,k4);
    ctemp=std::complex<double>(a75*H->h(),0);
    temp->add_wf(&ctemp,k5);
    ctemp=std::complex<double>(a76*H->h(),0);
    temp->add_wf(&ctemp,k6);


    k7->matrix_prod(dpsi_mat,temp);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(ddc1*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(ddc2*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(ddc3*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(ddc4*H->h(),0);
    temp->add_wf(&ctemp,k4);
    ctemp=std::complex<double>(ddc5*H->h(),0);
    temp->add_wf(&ctemp,k5);
    ctemp=std::complex<double>(ddc6*H->h(),0);
    temp->add_wf(&ctemp,k6);
    ctemp=std::complex<double>(ddc7*H->h(),0);
    temp->add_wf(&ctemp,k7);
    Psi0->set_wf(temp);

    delete k1;
    delete k2;
    delete k3;
    delete k4;
    delete k5;
    delete k6;
    delete k7;
    delete temp;
    delete temp2;
    delete [] dpsi_mat;
}
//##########################################################################
//
//##########################################################################
bool t_deriv_matrix(hamilton_matrix *H,wavefunction **dPsi_mat,double time_index)
{
   double* vector=new double[3];
   H->electric_field(time_index,vector);
   double efield_magnitude(sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]));
   bool condition(0);
   std::complex<double> ctemp(0,-1);

#pragma omp parallel for
      for(int i=0;i!=(dPsi_mat[0]->n_states_neut())*(dPsi_mat[0]->gsize_x());i++)
      {//LOOP OVER THE LINES OF THE H MATRIX
         int j(0);
         j=(i%dPsi_mat[0]->gsize_x()-2)*bool(i%dPsi_mat[0]->gsize_x()-2 >= 0);
      //WE USE FINITE DIFFERENCE METHOD TO FOURTH ORDER. ONLY THE 5 DIAGONALS ARE NON-ZERO. WE BEGIN AT i-2
         //std::cout<<"probe1 : "<<i<<" , "<<j<<std::endl;
         wavefunction* temp=new wavefunction (dPsi_mat[0]->gsize_x(),dPsi_mat[0]->n_states_neut(),dPsi_mat[0]->n_states_cat(),dPsi_mat[0]->n_states_cont());

         for(int s=0;s!=dPsi_mat[0]->n_states_neut();s++)
         {//LOOP OVER THE ELECTRONIC STATES OF THE NEUTRAL, COLUMN OF H
            while( (j%dPsi_mat[0]->gsize_x()) <= i%dPsi_mat[0]->gsize_x()+2)
            {//LOOP OVER THE GRID POINTS AROUND THE DIAGONAL
               //std::cout<<" j = "<<j<<" : "<<j%Psi->gsize_x()<<"/"<<(i%Psi->gsize_x()+2)<<" => "<<bool( (j%Psi->gsize_x()) <= (i%Psi->gsize_x()+2))<<std::endl;
               temp->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
               if(!(j%dPsi_mat[0]->gsize_x() == dPsi_mat[0]->gsize_x()-1))
                  j++;
               else
               {
                  j++;
                  break;
               }
            }
            j += dPsi_mat[0]->gsize_x() - 5 + 2*bool(i%dPsi_mat[0]->gsize_x() == 0) + bool (i%dPsi_mat[0]->gsize_x() == 1) + 2*bool(i%dPsi_mat[0]->gsize_x() == dPsi_mat[0]->gsize_x()-1) + bool (i%dPsi_mat[0]->gsize_x() == dPsi_mat[0]->gsize_x()-2);
            //std::cout<<"probe2 j = "<<j<<std::endl;
         }
         j+=2;
         if(efield_magnitude > H->efield_thresh())
         {
            for(int s=0;s!=dPsi_mat[0]->n_states_cat()*dPsi_mat[0]->n_states_cont();s++)
            {
               temp->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
               j += dPsi_mat[0]->gsize_x();
            }
         }
         dPsi_mat[i]->add_wf(&ctemp,temp);

         delete temp;
      }

#pragma omp parallel for
      for(int i=(dPsi_mat[0]->n_states_neut())*(dPsi_mat[0]->gsize_x());i!=(dPsi_mat[0]->n_states_neut())*(dPsi_mat[0]->gsize_x())+(dPsi_mat[0]->n_states_cat())*(dPsi_mat[0]->n_states_cont())*(dPsi_mat[0]->gsize_x());i++)
      {
         int j(0);
         wavefunction* temp=new wavefunction (dPsi_mat[0]->gsize_x(),dPsi_mat[0]->n_states_neut(),dPsi_mat[0]->n_states_cat(),dPsi_mat[0]->n_states_cont());
         j=i%dPsi_mat[0]->gsize_x();
         if(efield_magnitude > H->efield_thresh())
         {
            for(int s=0;s!=dPsi_mat[0]->n_states_neut();s++)
            {
               temp->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
               j += dPsi_mat[0]->gsize_x();
            }
         }
         for(int s=0;s!=dPsi_mat[0]->n_states_cat()*dPsi_mat[0]->n_states_cont();s++)
         {
            while(j%dPsi_mat[0]->gsize_x() <= (i+2)%dPsi_mat[0]->gsize_x())
            {
               temp->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
               if(!(j%dPsi_mat[0]->gsize_x() == dPsi_mat[0]->gsize_x()-1))
                  j++;
               else
               {
                  j++;
                  break;
               }
            }
            j += dPsi_mat[0]->gsize_x() - 5 + 2*bool(i%dPsi_mat[0]->gsize_x() == 0) + bool (i%dPsi_mat[0]->gsize_x() == 1) + 2*bool(i%dPsi_mat[0]->gsize_x() == dPsi_mat[0]->gsize_x()-1) + bool (i%dPsi_mat[0]->gsize_x() ==dPsi_mat[0]->gsize_x()-2);
         }
         dPsi_mat[i]->add_wf(&ctemp,temp);
         delete temp;
      }
}
//##########################################################################
//
//##########################################################################
bool Runge_kutta(wavefunction *Psi0,hamilton_matrix *H, int time_index)
{
/*    std::complex<double> **k1=new std::complex<double> * [2];
    k1[0]=new std::complex<double>[Psi0->n_states_neut()*Psi0->gsize_x()];
    k1[1]=new std::complex<double>[Psi0->n_states_cat()*Psi0->n_states_cont()*Psi0->gsize_x()];
    std::complex<double> **k2=new std::complex<double> * [2];
    k2[0]=new std::complex<double>[Psi0->n_states_neut()*Psi0->gsize_x()];
    k2[1]=new std::complex<double>[Psi0->n_states_cat()*Psi0->n_states_cont()*Psi0->gsize_x()];
    std::complex<double> **k3=new std::complex<double> * [2];
    k3[0]=new std::complex<double>[Psi0->n_states_neut()*Psi0->gsize_x()];
    k3[1]=new std::complex<double>[Psi0->n_states_cat()*Psi0->n_states_cont()*Psi0->gsize_x()];
    std::complex<double> **k4=new std::complex<double> * [2];
    k4[0]=new std::complex<double>[Psi0->n_states_neut()*Psi0->gsize_x()];
    k4[1]=new std::complex<double>[Psi0->n_states_cat()*Psi0->n_states_cont()*Psi0->gsize_x()];
    std::complex<double> **k5=new std::complex<double> * [2];
    k5[0]=new std::complex<double>[Psi0->n_states_neut()*Psi0->gsize_x()];
    k5[1]=new std::complex<double>[Psi0->n_states_cat()*Psi0->n_states_cont()*Psi0->gsize_x()];
    std::complex<double> **k6=new std::complex<double> * [2];
    k6[0]=new std::complex<double>[Psi0->n_states_neut()*Psi0->gsize_x()];
    k6[1]=new std::complex<double>[Psi0->n_states_cat()*Psi0->n_states_cont()*Psi0->gsize_x()];
    std::complex<double> **k7=new std::complex<double> * [2];
    k7[0]=new std::complex<double>[Psi0->n_states_neut()*Psi0->gsize_x()];
    k7[1]=new std::complex<double>[Psi0->n_states_cat()*Psi0->n_states_cont()*Psi0->gsize_x()];
*/

    wavefunction* k1=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k2=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k3=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k4=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k5=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k6=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* k7=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());

    wavefunction* temp=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());
    wavefunction* temp2=new wavefunction(Psi0->gsize_x(),Psi0->n_states_neut(),Psi0->n_states_cat(),Psi0->n_states_cont());

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

//    std::cout<<Psi0->show_neut_psi(21,0).real()<<std::endl;

    t_deriv(Psi0,H,k1,time_index);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a21*H->h(),0);
    temp->add_wf(&ctemp,k1);
//#pragma omp parallel for
/*
    for(int i=0;i!=(Psi0->n_states_neut()+Psi0->n_states_cat()*Psi0->n_states_cont())*Psi0->gsize_x();i++)
    {
       if(i < Psi0->n_states_neut()*Psi0->gsize_x())
       {
          state_index=int(i/Psi0->gsize_x());
          grid_index=int(i%Psi0->gsize_x());
          state_index_cont=-1;
          k1[0][i]=temp.show_neut_psi(grid_index,state_index);
          temp.set_neut_psi(state_index,grid_index,Psi0->show_neut_psi(grid_index,state_index)+a21*H->h()*k1[0][i]);
       }
       else
       {
          state_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())/(Psi0->gsize_x()*Psi0->n_states_cont()));
          state_index_cont=int(((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont()))/Psi0->gsize_x());
          grid_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont())%Psi0->gsize_x());
          k1[1][i]=temp.show_cat_psi(state_index,state_index_cont,grid_index);
          temp.set_cat_psi(state_index,state_index_cont,grid_index,Psi0->show_cat_psi(state_index,state_index_cont,grid_index)+a21*H->h()*k1[1][i]);
       }

    }*/
    t_deriv(temp,H,k2,time_index+1);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a31*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a32*H->h(),0);
    temp->add_wf(&ctemp,k2);
//#pragma omp parallel for
/*
    for(int i=0;i!=(Psi0->n_states_neut()+Psi0->n_states_cat()*Psi0->n_states_cont())*Psi0->gsize_x();i++)
    {
       if(i < Psi0->n_states_neut()*Psi0->gsize_x())
       {
          state_index=int(i/Psi0->gsize_x());
          grid_index=int(i%Psi0->gsize_x());
          state_index_cont=-1;
          k2[0][i]=temp2.show_neut_psi(grid_index,state_index);
          temp.set_neut_psi(state_index,grid_index,Psi0->show_neut_psi(grid_index,state_index)+a31*H->h()*k1[0][i]+a32*H->h()*k2[0][i]);
       }
       else
       {
          state_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())/(Psi0->gsize_x()*Psi0->n_states_cont()));
          state_index_cont=int(((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont()))/Psi0->gsize_x());
          grid_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont())%Psi0->gsize_x());
          k2[1][i]=temp2.show_cat_psi(state_index,state_index_cont,grid_index);
          temp.set_cat_psi(state_index,state_index_cont,grid_index,Psi0->show_cat_psi(state_index,state_index_cont,grid_index)+a31*H->h()*k1[1][i]+a32*H->h()*k2[1][i]);
       }
    }*/
    t_deriv(temp,H,k3,time_index+.5);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a41*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a42*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a43*H->h(),0);
    temp->add_wf(&ctemp,k3);
    /*
//#pragma omp parallel for
    for(int i=0;i!=(Psi0->n_states_neut()+Psi0->n_states_cat()*Psi0->n_states_cont())*Psi0->gsize_x();i++)
    {
       if(i < Psi0->n_states_neut()*Psi0->gsize_x())
       {
          state_index=int(i/Psi0->gsize_x());
          grid_index=int(i%Psi0->gsize_x());
          state_index_cont=-1;
          k3[0][i]=temp2.show_neut_psi(grid_index,state_index);
          temp.set_neut_psi(state_index,grid_index,Psi0->show_neut_psi(grid_index,state_index)+a41*H->h()*k1[0][i]+a42*H->h()*k2[0][i]+a43*H->h()*k3[0][i]);
       }
       else
       {
          state_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())/(Psi0->gsize_x()*Psi0->n_states_cont()));
          state_index_cont=int(((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont()))/Psi0->gsize_x());
          grid_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont())%Psi0->gsize_x());
          k3[1][i]=temp2.show_cat_psi(state_index,state_index_cont,grid_index);
          temp.set_cat_psi(state_index,state_index_cont,grid_index,Psi0->show_cat_psi(state_index,state_index_cont,grid_index)+a41*H->h()*k1[1][i]+a42*H->h()*k2[1][i]+a43*H->h()*k3[1][i]);
       }
    }*/
    t_deriv(temp,H,k4,time_index+2./3.);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a51*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a52*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a53*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(a54*H->h(),0);
    temp->add_wf(&ctemp,k4);
//#pragma omp parallel for
/*
    for(int i=0;i!=(Psi0->n_states_neut()+Psi0->n_states_cat()*Psi0->n_states_cont())*Psi0->gsize_x();i++)
    {
       if(i < Psi0->n_states_neut()*Psi0->gsize_x())
       {
          state_index=int(i/Psi0->gsize_x());
          grid_index=int(i%Psi0->gsize_x());
          state_index_cont=-1;
          k4[0][i]=temp2.show_neut_psi(grid_index,state_index);
          temp.set_neut_psi(state_index,grid_index,Psi0->show_neut_psi(grid_index,state_index)+a51*H->h()*k1[0][i]+a52*H->h()*k2[0][i]+a53*H->h()*k3[0][i]+a54*H->h()*k4[0][i]);
       }
       else
       {
          state_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())/(Psi0->gsize_x()*Psi0->n_states_cont()));
          state_index_cont=int(((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont()))/Psi0->gsize_x());
          grid_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont())%Psi0->gsize_x());
          k4[1][i]=temp2.show_cat_psi(state_index,state_index_cont,grid_index);
          temp.set_cat_psi(state_index,state_index_cont,grid_index,Psi0->show_cat_psi(state_index,state_index_cont,grid_index)+a51*H->h()*k1[1][i]+a52*H->h()*k2[1][i]+a53*H->h()*k3[1][i]+a54*H->h()*k4[1][i]);
       }
    }*/
    t_deriv(temp,H,k5,time_index+(7.-sqrt(21))/14);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a61*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a62*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a63*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(a64*H->h(),0);
    temp->add_wf(&ctemp,k4);
    ctemp=std::complex<double>(a65*H->h(),0);
    temp->add_wf(&ctemp,k5);
//#pragma omp parallel for
/*
    for(int i=0;i!=(Psi0->n_states_neut()+Psi0->n_states_cat()*Psi0->n_states_cont())*Psi0->gsize_x();i++)
    {
       if(i < Psi0->n_states_neut()*Psi0->gsize_x())
       {
          state_index=int(i/Psi0->gsize_x());
          grid_index=int(i%Psi0->gsize_x());
          state_index_cont=-1;
          k5[0][i]=temp2.show_neut_psi(grid_index,state_index);
          temp.set_neut_psi(state_index,grid_index,Psi0->show_neut_psi(grid_index,state_index)+a61*H->h()*k1[0][i]+a62*H->h()*k2[0][i]+a63*H->h()*k3[0][i]+a64*H->h()*k4[0][i]+a65*H->h()*k5[0][i]);
       }
       else
       {
          state_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())/(Psi0->gsize_x()*Psi0->n_states_cont()));
          state_index_cont=int(((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont()))/Psi0->gsize_x());
          grid_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont())%Psi0->gsize_x());
          k5[1][i]=temp2.show_cat_psi(state_index,state_index_cont,grid_index);
          temp.set_cat_psi(state_index,state_index_cont,grid_index,Psi0->show_cat_psi(state_index,state_index_cont,grid_index)+a61*H->h()*k1[1][i]+a62*H->h()*k2[1][i]+a63*H->h()*k3[1][i]+a64*H->h()*k4[1][i]+a65*H->h()*k5[1][i]);
       }
    }
    */
    t_deriv(temp,H,k6,time_index+(7.+sqrt(21))/14);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(a71*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(a72*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(a73*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(a74*H->h(),0);
    temp->add_wf(&ctemp,k4);
    ctemp=std::complex<double>(a75*H->h(),0);
    temp->add_wf(&ctemp,k5);
    ctemp=std::complex<double>(a76*H->h(),0);
    temp->add_wf(&ctemp,k6);
//#pragma omp parallel for
/*
    for(int i=0;i!=(Psi0->n_states_neut()+Psi0->n_states_cat()*Psi0->n_states_cont())*Psi0->gsize_x();i++)
    {
       if(i < Psi0->n_states_neut()*Psi0->gsize_x())
       {
          state_index=int(i/Psi0->gsize_x());
          grid_index=int(i%Psi0->gsize_x());
          state_index_cont=-1;
          k6[0][i]=temp2.show_neut_psi(grid_index,state_index);
          temp.set_neut_psi(state_index,grid_index,Psi0->show_neut_psi(grid_index,state_index)+a71*H->h()*k1[0][i]+a72*H->h()*k2[0][i]+a73*H->h()*k3[0][i]+a74*H->h()*k4[0][i]+a75*H->h()*k5[0][i]+a76*H->h()*k6[0][i]);
       }
       else
       {
          state_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())/(Psi0->gsize_x()*Psi0->n_states_cont()));
          state_index_cont=int(((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont()))/Psi0->gsize_x());
          grid_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont())%Psi0->gsize_x());
          k6[1][i]=temp2.show_cat_psi(state_index,state_index_cont,grid_index);
          temp.set_cat_psi(state_index,state_index_cont,grid_index,Psi0->show_cat_psi(state_index,state_index_cont,grid_index)+a71*H->h()*k1[1][i]+a72*H->h()*k2[1][i]+a73*H->h()*k3[1][i]+a74*H->h()*k4[1][i]+a75*H->h()*k5[1][i]+a76*H->h()*k6[1][i]);
       }
    }
    */
    t_deriv(temp,H,k7,time_index+1);
    temp->set_wf(Psi0);
    ctemp=std::complex<double>(ddc1*H->h(),0);
    temp->add_wf(&ctemp,k1);
    ctemp=std::complex<double>(ddc2*H->h(),0);
    temp->add_wf(&ctemp,k2);
    ctemp=std::complex<double>(ddc3*H->h(),0);
    temp->add_wf(&ctemp,k3);
    ctemp=std::complex<double>(ddc4*H->h(),0);
    temp->add_wf(&ctemp,k4);
    ctemp=std::complex<double>(ddc5*H->h(),0);
    temp->add_wf(&ctemp,k5);
    ctemp=std::complex<double>(ddc6*H->h(),0);
    temp->add_wf(&ctemp,k6);
    ctemp=std::complex<double>(ddc7*H->h(),0);
    temp->add_wf(&ctemp,k7);
    Psi0->set_wf(temp);

//    std::cout<<"Then "<<Psi0->show_neut_psi(21,0).real()<<std::endl;
//#pragma omp parallel for
/*
    for(int i=0;i!=(Psi0->n_states_neut()+Psi0->n_states_cat()*Psi0->n_states_cont())*Psi0->gsize_x();i++)
    {
       if(i < Psi0->n_states_neut()*Psi0->gsize_x())
       {
          state_index=int(i/Psi0->gsize_x());
          grid_index=int(i%Psi0->gsize_x());
          state_index_cont=-1;
          k7[0][i]=temp2.show_neut_psi(grid_index,state_index);
          Psi0->set_neut_psi(state_index,grid_index,Psi0->show_neut_psi(grid_index,state_index)+ddc1*H->h()*k1[0][i]+ddc2*H->h()*k2[0][i]+ddc3*H->h()*k3[0][i]+ddc4*H->h()*k4[0][i]+ddc5*H->h()*k5[0][i]+ddc6*H->h()*k6[0][i]+ddc7*H->h()*k7[0][i]);
       }
       else
       {
          state_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())/(Psi0->gsize_x()*Psi0->n_states_cont()));
          state_index_cont=int(((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont()))/Psi0->gsize_x());
          grid_index=int((i-Psi0->n_states_neut()*Psi0->gsize_x())-state_index*(Psi0->gsize_x()*Psi0->n_states_cont())%Psi0->gsize_x());
          k7[1][i]=temp2.show_cat_psi(state_index,state_index_cont,grid_index);
          Psi0->set_cat_psi(state_index,state_index_cont,grid_index,Psi0->show_cat_psi(state_index,state_index_cont,grid_index)+ddc1*H->h()*k1[1][i]+ddc2*H->h()*k2[1][i]+ddc3*H->h()*k3[1][i]+ddc4*H->h()*k4[1][i]+ddc5*H->h()*k5[1][i]+ddc6*H->h()*k6[1][i]+ddc7*H->h()*k7[1][i]);
       }
    }
    */
/*
     delete [] k1[0];
     delete [] k1[1];
     delete [] k1;
     delete [] k2[0];
     delete [] k2[1];
     delete [] k2;
     delete [] k3[0];
     delete [] k3[1];
     delete [] k3;
     delete [] k4[0];
     delete [] k4[1];
     delete [] k4;
     delete [] k5[0];
     delete [] k5[1];
     delete [] k5;
     delete [] k6[0];
     delete [] k6[1];
     delete [] k6;
     delete [] k7[0];
     delete [] k7[1];
     delete [] k7;
     */
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

bool t_deriv(wavefunction *Psi,hamilton_matrix *H,wavefunction *dPsi,double time_index)
{

   double* vector=new double[3];
   H->electric_field(time_index,vector);
   double efield_magnitude(sqrt(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]));
   bool condition(0);
   //std::cout<<"entering time deriv routine "<<std::endl;

   //((Psi->n_states_neut())*(Psi->gsize_x())+(Psi->n_states_cat())*(Psi->n_states_cont())*(Psi->gsize_x()))
#pragma omp parallel for
      for(int i=0;i!=(Psi->n_states_neut())*(Psi->gsize_x());i++)
      {
         int j(0);
         int state_index;
         int grid_index;
         int state_index_cont;
         wavefunction* Hvec=new wavefunction (Psi->gsize_x(),Psi->n_states_neut(),Psi->n_states_cat(),Psi->n_states_cont());

         H->show_indexes(i,i,&state_index,&grid_index,&state_index_cont,&state_index,&grid_index,&state_index_cont);

//#pragma omp parallel for
         for(int s=0;s!=Psi->n_states_neut();s++)
         {
               j=(s*(Psi->gsize_x())+(grid_index-2))*bool(grid_index-2 >= 0);

               while( (j%Psi->gsize_x()) <= grid_index + 2 && j%Psi->gsize_x() != Psi->gsize_x()-1 ) 
               {
                  Hvec->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
                  j++;
               }
               //j += Psi->gsize_x() - 5 + 2*bool(i%Psi->gsize_x() == 0) + bool (i%Psi->gsize_x() == 1) + 2*bool(i%Psi->gsize_x() == Psi->gsize_x()-1) + bool (i%Psi->gsize_x() ==Psi->gsize_x()-2);
         }
/*         j+=2;
         if(efield_magnitude >= H->efield_thresh())
         {
            for(int s=0;s!=Psi->n_states_cat()*Psi->n_states_cont();s++)
            {
               Hvec->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
               j += Psi->gsize_x();
            }
         }*/
         //std::cout<<"Then "<<Hvec->show_neut_psi(21,0).real()<<std::endl;
         dPsi->set_psi_elwise(i,std::complex<double>(0,-1)*Psi->dot_prod(Hvec));

         delete Hvec;
      }
/*         for(int m=0;m!=Psi->gsize_x();m++)
         {
            std::cout<<(H->xmin()+m*(H->xmax()-H->xmin())/(Psi->gsize_x()))*0.529<<"  ,"<<H->pot_neut(1,m)<<" ,"<<norm(Psi->show_neut_psi(m,1))<<","<<norm(dPsi->show_neut_psi(m,1))<<std::endl;
         }std::cout<<"##########################################################"<<std::endl;
*/
      for(int i=(Psi->n_states_neut())*(Psi->gsize_x());i!=(Psi->n_states_neut())*(Psi->gsize_x())+(Psi->n_states_cat())*(Psi->n_states_cont())*(Psi->gsize_x());i++)
      {
         int j(0);
         wavefunction* Hvec=new wavefunction (Psi->gsize_x(),Psi->n_states_neut(),Psi->n_states_cat(),Psi->n_states_cont());
         j=i%Psi->gsize_x();
         if(efield_magnitude > H->efield_thresh())
         {
            for(int s=0;s!=Psi->n_states_neut();s++)
            {
               Hvec->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
               j += Psi->gsize_x();
            }
         }
         for(int s=0;s!=Psi->n_states_cat()*Psi->n_states_cont();s++)
         {
            while(j%Psi->gsize_x() <= (i+2)%Psi->gsize_x())
            {
               Hvec->set_psi_elwise(j,H->hamilt_element(time_index,i,j));
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
         dPsi->set_psi_elwise(i,std::complex<double>(0,-1)*Psi->dot_prod(Hvec));
         delete Hvec;
      }
    //std::cout<<"leaving time deriv routine "<<std::endl;
    return 0;
}
//##########################################################################
//
//##########################################################################
void propagate(wavefunction *Psi, hamilton_matrix *H,int* time_index,int num_of_loop)
{
   for(int i=0;i!=num_of_loop;i++)
   {
      //std::cout<<"loop "<<i<<std::endl;
      Runge_kutta(Psi,H,*time_index);
      //Runge_kutta_notdH(Psi,H,*time_index);
      *time_index=*time_index+1;
   }
}

