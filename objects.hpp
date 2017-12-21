#include <cstdlib>
#include <cmath>
#include <complex>

class hamilton_matrix;

class wavefunction {
   public:
      //constructor
      wavefunction(int gsize_x,int n_states_neut,int n_states_cat,int n_states_cont);
      ~wavefunction();
      void set_neut_psi(int state_index,int grid_index,std::complex<double> value);
      void set_cat_psi(int state_cat_index,int states_cont_index,int grid_index,std::complex<double> value);
      std::complex<double> show_neut_psi(int grid_index,int state_index);
      std::complex<double> show_cat_psi(int grid_index,int state_cat_index,int states_cont_index);
      int gsize_x();
      int n_states_neut();
      int n_states_cat();
      int n_states_cont();
      void initialize(hamilton_matrix* H);
      void set_norm(double value);
      double norm();
      void set_dipole(double *vector);
      void show_dipole(double *vector);
      void add_wf(std::complex<double>* a,wavefunction* y=NULL);
      void set_wf(wavefunction* x);
      void set_psi_elwise(int i,std::complex<double> val);
      void wf_vec(std::complex<double>* neut_vec,std::complex<double>* cat_vec);
      std::complex<double> dot_prod(wavefunction* Bra);
      void matrix_prod(wavefunction* mat,wavefunction* ket);
   private:
      int m_gsize_x;
      int m_n_states_neut;
      int m_n_states_cat;
      int m_n_states_cont;
      std::complex<double> *m_neut_part;//neutral part of the wavefunction vector
      std::complex<double> *m_cat_part;//cation and continuum part of the wavefunction vector
      double m_norm;
      double m_ion_pop;
      double *m_dipole;
};

class hamilton_matrix {
   private:
      //GRID AND STATES PARAMETERS
      int m_gsize_x;
      int m_n_states_neut;
      int m_n_states_cat;
      int m_n_states_cont;
      int m_n_k;
      int m_n_angles;
      double m_xmin;
      double m_xmax;
      double m_mass;
      //ELECTRIC FIELD THRESHOLD
      double m_efield_thresh;
      //TIME VARIABLE SETTINGS
      int m_n_times;
      double m_h;
      //ELECTRONIC STRUCTURE ARRAYS
      double *m_pot_neut;
      double *m_pot_cat;
      double **m_dmx_neut;
      double **m_dmy_neut;
      double **m_dmz_neut;
      double **m_dmx_cat;
      double **m_dmy_cat;
      double **m_dmz_cat;
      std::complex<double> **m_PICE_x;
      std::complex<double> **m_PICE_y;
      std::complex<double> **m_PICE_z;
      double **k_orientation;
      double *k_modulus;
      double **m_NAC;
      //OTHER HAMILTONIAN MATRIX ARRAYS
      int **translation_vector;
      double min_distance;
      double *kinetic_energy;
      double *derivative_matrix;

   public:
      hamilton_matrix(int gsize_x,int n_states_neut,int n_states_cat,int n_k,int n_angles,double xmin,double xmax,double mass,int n_times,double h,double efield_thresh);
      ~hamilton_matrix();
      void set_pot_neut(std::string file_address);
      void set_pot_cat(std::string file_address);
      void set_dm_neut(std::string file_address);
      void set_dm_cat(std::string file_address);
      void set_PICE(std::string file_address);
      void set_NAC(std::string file_address);
      std::complex<double> hamilt_element(double time_index,int i,int j);
      double h();
      double efield_thresh();
      double electric_field(double time_index,double* vector);
      double potential_vector(double time_index,double* vector);
      double kinetic_energy_matrix(int i,int j);
      double pot_neut(int state_index,int grid_index);
      double energy(wavefunction* Psi,double time_index);
      void rescale_pot(double min_pot);
      void show_indexes(int index1,int index2,int *state_index_1,int *grid_index_1,int *state_index_cont_1,int *state_index_2,int *grid_index_2,int *state_index_cont_2);
};
#include "obj_wavefunction.cpp"
#include "obj_hamilt_mat_init.cpp"
#include "obj_hamilt_mat_comp.cpp"
#include "electric_field_func.cpp"
