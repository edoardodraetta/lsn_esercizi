
class MolDyn{

public:

	MolDyn(std::string input = "input.dat");
	~MolDyn();

	void Simulate();
	void PrintStats(std::string datadir = "../../data");

private:

	void Welcome();
	void Initialize();
	void Restart();
	void Rescale();
	void PrepareVelocities();

	void ConfFinal();
	void ConfXYZ(int);
	void ConfOut(std::string);

	void Move();
	void Measure();
	double Force(int, int);
	double Pbc(double);

	// simulation
	std::string inputfile;
	int seed {1};
	double delta;
	int nstep, iprint;

	// parameters, observables
	const int m_props = 4;
	int n_props {4};
	int iv{0}, ik{1}, it{2}, ie{3}; // index for each observable
	double stima_pot, stima_kin, stima_etot, stima_temp;

	// statistics
	static const int nblocks = 100;
	double ave_etot[nblocks] = {}, av2_etot[nblocks] = {};
	double ave_ekin[nblocks] = {}, av2_ekin[nblocks] = {};
	double ave_epot[nblocks] = {}, av2_epot[nblocks] = {};
	double ave_temp[nblocks] = {}, av2_temp[nblocks] = {};

	// averages
	double acc,att;

	// thermodynamical state
	int npart;
	double energy, temp, vol, rho, box, rcut;

	// configuration
	static const int m_part = 108;
	double x[m_part],y[m_part],z[m_part];
	double xold[m_part],yold[m_part],zold[m_part];
	double vx[m_part],vy[m_part],vz[m_part];

	// options (new)
	bool restart, quick_rescale;
};
