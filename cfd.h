#ifndef CFDSIM_CFD_H_
#define CFDSIM_CFD_H_

class cfd
{
	public:
	cfd(float _lx, float _ly, int _nx, int _ny, int _GS, int _IOP, float _h, float _grav);
	~cfd();

	float* getColorMap();
	float* getSourceMap();
	float* getObstructionMap();
	float* getColorMapWithBoundary();
	float* getDivergenceSourceMap();

	void updateTimeStep(float _h);
	void updateGrav(float _grav);

	void getVelocity(int i, int j, float &u, float &v);
	void getDensity(int i, int j, float &den);
	void getColor(int i, int j, int channel, float &col);
	void getDivergence(int i, int j, float &div);
	void getPressure(int i, int j, float &pre);

	void swapColMap();
	void swapDenMap();
	void swapVelMap();

	void interpolateVelocity(float x, float y, float &u, float &v);
	void interpolateDensity(float x, float y, float &res);
	void interpolateColor(float x, float y, int channel, float &res);

	void advectDensity();
	void advectVelocity();
	void advectColor();
	void bouyancy();
	void addForces();
	void drawSourcesToDensity();
	void modifyVelocity();
	void computeDivergence();
	void computePressure();
	void modifyBoundary();
	void project();
	void updateFluid();
	void Initialize( float *data, int size, float value );

private:

	int nx, ny;
	int size;
	float lx, ly;
	int dx, dy;
	int GS, IOP;
	float h, grav;

	float *color_map_boundary;
	float *color_map;
	float *obstruction_map;
	float *source_map;
	float *density_map;
	float *velocity_map;
	float *divergence_map;
	float *divergence_source_map;
	float *pressure_map;

	float *swap_vel;
	float *swap_den;
	float *swap_col;

	float *temp_vel_map;
	float *temp_den_map;
	float *temp_col_map;
	float *temp_pre_map;
};

#endif // CFDSIM_CFD_H_