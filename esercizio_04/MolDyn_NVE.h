/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
int n_props;
int iv,ik,it,ie,iw,igofr;
double bin_size, nbins;
double stima_pot, stima_kin, stima_etot;
double stima_temp, stima_press, stima_gdir;
double err_pot, err_kin, err_etot;
double err_temp, err_press, err_gdir;
double current_temp;

// averages
const int m_props=1000;
double acc,att;
double walker[m_props];
double blk_av[m_props], blk_norm;
double glob_av[m_props],glob_av2[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblocks;
double delta;
bool restart, rescale;
double rescaletemp;

// functions
void Input(void);
void Welcome();
void Move(void);
void Measure(void);
double Force(int, int);
double Pbc(double);

// begin
void Prepare();
void Rescale();
void Restart();

// medie di blocco
void Reset(int);
void Accumulate(void);
void Averages(int);

// stats
double Error(double, double, int);

// conf
void ConfFinal(void);
void ConfPenult(void);
void ConfOut(std::string);
void ConfXYZ(int);




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
