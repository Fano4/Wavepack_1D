#include <cstdlib>
#include <cmath>
#include <complex>
#include <ctime>


class hamilton_matrix;

class wavefunction {
   public:
      //constructor
      wavefunction(int gsize_x,int tgsize_x,int n_states_neut,int n_states_cat,int n_states_cont);
      ~wavefunction();
      void set_neut_psi(int state_index,int grid_index,std::complex<double> value);
      void set_cat_psi(int state_cat_index,int states_cont_index,int grid_index,std::complex<double> value);
      std::complex<double> show_neut_psi(int grid_index,int state_index);
      std::complex<double> show_cat_psi(int grid_index,int state_cat_index,int states_cont_index);
      int gsize_x();
      int tgsize_x();
      int n_states_neut();
      int n_states_cat();
      int n_states_cont();
      void initialize(hamilton_matrix* H);
      void set_norm(double value);
      double norm(hamilton_matrix *H);
      void set_dipole(hamilton_matrix *H);
      void show_dipole(double *vector,bool species);
      void add_wf(std::complex<double>* a,wavefunction* y=NULL,bool cat=1);
      void set_wf(wavefunction* x,bool cat=1);
      void set_psi_elwise(int i,std::complex<double> val);
      void wf_vec(std::complex<double>* neut_vec,std::complex<double>* cat_vec);
      std::complex<double> dot_prod(wavefunction* Bra, hamilton_matrix *H);
      void matrix_prod(wavefunction** mat,wavefunction* ket,hamilton_matrix *H);
      double state_pop(bool species,int state,hamilton_matrix* H=NULL);
   private:
      int m_gsize_x;
      int m_tgsize_x;
      int m_n_states_neut;
      int m_n_states_cat;
      int m_n_states_cont;
      std::complex<double> *m_neut_part;//neutral part of the wavefunction vector
      std::complex<double> *m_cat_part;//cation and continuum part of the wavefunction vector
      double m_norm;
      double m_ion_pop;
      double *m_dipole_neut;
      double *m_dipole_cat;
};

class hamilton_matrix {
   private:
      //TEMPORARY SIGN CORRECTION ARRAY!!!

      double** sign_corr;

      //GRID AND STATES PARAMETERS
      int m_gsize_x;
      int m_tgsize_x;
      int m_small_gsize_x;
      int m_n_states_neut;
      int m_n_states_cat;
      int m_n_states_cont;
      int m_n_k;
      int m_n_angles;
      double m_xmin;
      double m_xmax;
      double m_mass;
      double m_kmin;
      double m_kmax;
//      double m_kymin;
//      double m_kymax;
//      double m_kzmin;
//      double m_kzmax;
//      int m_n_kx;
//      int m_n_ky;
//      int m_n_kz;
      //ELECTRIC FIELD THRESHOLD
      double m_pot_vec_thresh;
      double m_efield_thresh;
      double m_pot_vec_mod;
      double m_pot_vec_tm_mod;
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
      double *kinetic_energy;
      double *derivative_matrix;
      double *position_array_script;
      double *m_dk_vec;
      pice_set *pice_data;

   public:
      hamilton_matrix(int gsize_x,int tgsize_x,int small_gsize_x,int n_states_neut,int n_states_cat,int n_k,int n_angles,double kmin,double kmax,double xmin,double xmax,double mass,int n_times,double h,double efield_thresh,double pot_vec_thresh,std::string pice_data_loc);
      ~hamilton_matrix();
      void set_pot_neut(std::string file_address);
      void set_pot_cat(std::string file_address);
      void set_dm_neut(std::string file_address);
      void set_dm_cat(std::string file_address);
      double show_dm_neut(int state_index_1,int state_index_2,int grid_index,int component);
      double show_dm_cat(int state_index_1,int state_index_2,int grid_index,int component);
      void set_PICE(double *pot_vec=NULL);
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
      double xmin();
      double xmax();
      int n_k();
      double dk(int i);
      void set_phase(std::string file_address);
      void print_dipole_neut();
      void dk_vec(double *array);
      double show_nac(int state_index_1,int state_index_2,int grid_index_1,int grid_index_2);
//      void spherical_extract_from_cube(int neut_state,int cat_state,int r_index,int component,double *Recube,double *Imcube,double *pot_vec);
      void sphere_dist_gen(bool randiso=1,int n_phi=0);
//      bool cube_reader(std::string MO_cube_loc,double *cube_array,bool extract_dimensions=0);
//      double grid_k_cube_spacing(); 
      void set_pot_vec_mod(double value);
      void set_pot_vec_tm_mod(double value);
      double pot_vec_mod();
      double pot_vec_tm_mod();
      double pot_vec_thresh() const;
      void plot_integrated_cross_section(std::string file_address,int neut_state,int cat_state,int pos_index=-1);
      double k_mod_val(int k_index) const;
      double k_spher_orient(bool component, int index) const;
};

