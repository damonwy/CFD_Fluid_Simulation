#ifndef CFDSIM_CFD_H_
#define CFDSIM_CFD_H_

enum Scheme { BFECC, SL, MM };

class cfd
{
	public:
	cfd(float _lx, float _ly, int _nx, int _ny, int _GS, int _IOP, float _h, float _grav, Scheme _advect);
	~cfd();
	bool has_surfacetension;
	bool has_vorticity;
	float T, V;

	float* getColorMap();
	float* getSourceMap();
	float* getObstructionMap();
	float* getColorMapWithBoundary();
	float* getDivergenceSourceMap();

	void getVelocity(int i, int j, float &u, float &v, float *velocityMap);
	void getDensity(int i, int j, float &den, float *densityMap);
	void getColor(int i, int j, int channel, float &col, float *colorMap);
	void getDivergence(int i, int j, float &div);
	void getPressure(int i, int j, float &pre);
	void getAreaRotation(int i, int j, float &w);

	void interpolateVelocity(float x, float y, float &u, float &v, float *velocityMap);
	void interpolateDensity(float x, float y, float &res, float *densityMap);
	void interpolateColor(float x, float y, int channel, float &res, float *colorMap);

	void advectDensity(float *&target_density_map, float *&source_density_map, float *&temp_density_map, bool isBackWards);
	void advectVelocity(float *&target_velocity_map, float *&source_velocity_map, float *&temp_velocity_map, bool isBackWards);
	void advectColor(float *&target_color_map, float *&source_color_map, float *&temp_color_map, bool isBackWards);

	void advectDensityScheme(Scheme s);
	void advectVelocityScheme(Scheme s);
	void advectColorScheme(Scheme s);


	void bouyancy();
	void surface_tension();
	void vorticity_confinement();
	void addForces();
	
	void computeDivergence();
	void computePressure();
	void modifyBoundary();
	void modifyVelocity();
	void project();
	void updateFluid();

	void drawSourcesToDensity();
	void updateTimeStep(float _h);
	void updateGrav(float _grav);
	void updateScheme(Scheme s);

	void Initialize( float *data, int size, float value );

private:
	
	int nx, ny;
	int size;
	float lx, ly;
	int dx, dy;
	int GS, IOP;
	float h, grav;
	Scheme advection_scheme;


	float *color_map_boundary;
	float *obstruction_map;
	float *source_map;

	float *density_map;
	float *density_map_b;
	float *density_map_f;
	float *density_map_err;
	float *density_map_bfe;
	float *temp_density_map_b;
	float *temp_density_map_f;


	float *velocity_map;
	float *velocity_map_f;
	float *velocity_map_b;
	float *velocity_map_err;
	float *velocity_map_bfe;
	float *temp_velocity_map_b;
	float *temp_velocity_map_f;

	float *color_map;
	float *color_map_f;
	float *color_map_b;
	float *color_map_err;
	float *color_map_bfe;
	float *temp_color_map_b;
	float *temp_color_map_f;

	float *divergence_map;
	float *divergence_source_map;
	float *pressure_map;

	float *w_map;

	float *temp_vel_map;
	float *temp_den_map;
	float *temp_col_map;
	float *temp_pre_map;
};

#endif // CFDSIM_CFD_H_