#ifndef __MOLDYN__
#define __MOLDYN__

class MolDyn{

private:
	// parameters, observables
	const int m_props = 4;
	int n_props;
	int iv, ik, it, ie; // indices for each observable
	double stima_pot, stima_kin, stima_etot, stima_temp;

	// averages
	double acc,att;

	// configuration (--> changed to vectors)
	static const int m_part = 108;
	double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
	double vx[m_part],vy[m_part],vz[m_part];

	// thermodynamical state
	int npart;
	double energy, temp, vol, rho, box, rcut;

	// simulation
	int seed;
	double delta;
	bool restart, quick_rescale; // new


protected:

public:

	MolDyn();
	~MolDyn();

	int nstep, iprint;

	void Input();
	void Move();

	void ConfFinal();
	void ConfXYZ(int);
	void Measure();
	double Force(int, int);
	double Pbc(double);

	// new
	void ConfOut(std::string);
	void PrepareVelocities();
	void Restart();
	void Rescale();

	// blocked stats
	static const int nblocks = 100;
	double ave_epot[nblocks], ave_ekin[nblocks], ave_etot[nblocks], ave_temp[nblocks];
 	double av2_epot[nblocks], av2_ekin[nblocks], av2_etot[nblocks], av2_temp[nblocks];
 	double err_epot[nblocks], err_ekin[nblocks], err_etot[nblocks], err_temp[nblocks];
 	void BlockedStats();
};

#endif
