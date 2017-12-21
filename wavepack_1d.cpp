#include "wavepack_1d.hpp"

int main( int argc, char * argv [])
{
    using namespace std;

    omp_set_num_threads(16);
    //PATH LINKING OTHER FILES
    string neutral_pes("/data1/home/stephan/LiH_512_points/LiH_neut_");
    string cation_pes("/data1/home/stephan/LiH_512_points/LiH_cat_");
    string neutral_dipole("/data1/home/stephan/LiH_512_points/LiH_");
    //string ionization_coupling_file="LiH_PICE_";//LiH_PICE_0_0.txt
    //string out_file="Output.log";
    //string read_file="gnu-out.txt";
    //string elec_wavepack_file="LiH_elec_wvpck.txt";
    //PARAMETERS OF THE SIMULATION
    int gsize_x(512);
    int n_states_neut(15);
    int n_states_cat(0);
    int n_angles(0);
    int n_k(0);
    double xmin(0.8/0.529);//!!! THESE VALUES ARE IN ATOMIC UNITS AND NOT IN ANGSTROM
    double xmax(21.6/0.529);
    double mass(1836*(1.007825*6.015122795/(1.007825+6.015122795)));
    double total_time(100/0.02418884);
    double h(0.001/0.02418884);
    int n_times(int(total_time/h));
    int time_index(0);
    double dipole[3];
    double efield_thresh(1e-5);

    wavefunction* Psi= new wavefunction(gsize_x, n_states_neut,n_states_cat,n_angles*n_k);
    hamilton_matrix* H=new hamilton_matrix(gsize_x,n_states_neut,n_states_cat,n_k,n_angles,xmin,xmax,mass,n_times,h,efield_thresh);

    H->set_pot_neut(neutral_pes.c_str());
    H->set_pot_cat(cation_pes.c_str());
    H->set_dm_neut(neutral_dipole.c_str());

    Psi->initialize(H);

    for(int i=0;i!=gsize_x;i++)
    {
       std::cout<<H->pot_neut(0,i)<<","<<Psi->show_neut_psi(i,0).real()<<","<<Psi->show_neut_psi(i,0).imag()<<std::endl;
    }

    std::cout<<"initial energy of the system "<<setprecision(15)<<H->energy(Psi,0)<<std::endl<<"initial norm of the system ";
    std::cout<<setprecision(15)<<Psi->norm()<<std::endl;

    wavefunction* dPsi= new wavefunction(gsize_x, n_states_neut,n_states_cat,n_angles*n_k);

    while(time_index*h <= total_time)
    {
       propagate(Psi,H,&time_index,10);
       std::cout<<"Time "<<time_index*h*0.02418884<<std::endl;
       std::cout<<"Energy of the system "<<setprecision(15)<<H->energy(Psi,0)<<std::endl;
       Psi->show_dipole(dipole);
       std::cout<<"Total dipole"<<setprecision(15)<<dipole[0]<<","<<dipole[1]<<","<<dipole[2]<<std::endl;
       std::cout<<"Norm of the system "<<setprecision(15)<<Psi->norm()<<std::endl;
    }

    delete Psi;
    delete dPsi;
    delete H;


   return 0;
}
