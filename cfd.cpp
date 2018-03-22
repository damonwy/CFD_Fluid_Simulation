#include <iostream>
#include "cfd.h"

cfd::cfd(float _lx, float _ly, int _nx, int _ny, int _GS, int _IOP, float _h, float _grav){
    nx = _nx;
    ny = _ny;
    lx = _lx;
    ly = _ly;
    GS = _GS;
    IOP = _IOP;
    h = _h;
    grav = _grav;

    dx = lx / (nx - 1);
    dy = ly / (ny - 1);
    size = nx * ny;

    color_map = new float[3*size];
    color_map_boundary = new float[3*size];

    source_map = new float[size];
    Initialize(source_map, size, 0.0);

    density_map = new float[size];
    Initialize(density_map, size, 0.0);

    obstruction_map = new float[size];
    Initialize(obstruction_map, size, 1.0);

    velocity_map = new float[2*size];
    Initialize(velocity_map, 2*size, 0.0);

    pressure_map = new float[size];
    Initialize(pressure_map, size, 0.0);

    divergence_map = new float[size];
    Initialize(divergence_map, size, 0.0);

    divergence_source_map = new float[size];
    Initialize(divergence_source_map, size, 0.0);

    temp_vel_map = new float[2*size];
    temp_den_map = new float[size];
    temp_col_map = new float[3*size];
    temp_pre_map = new float[size];
    
    swap_vel = new float[2*size];
    swap_den = new float[size];
    swap_col = new float[3*size];
}

void cfd::updateTimeStep(float _h){
    h = _h;
}

void cfd::updateGrav(float _grav){
    grav = _grav;
}

float* cfd::getColorMap(){
    return color_map;
}

float* cfd::getSourceMap(){
    return source_map;
}

float* cfd::getObstructionMap(){
    return obstruction_map;
}

float* cfd::getDivergenceSourceMap(){
    return divergence_source_map;
}


float* cfd::getColorMapWithBoundary(){
    for(int j=0; j<ny; j++){
      for(int i = 0; i < nx; i++){
          int index = i + nx * j;
          color_map_boundary[3*index+0] = color_map[3*index+0] * obstruction_map[index];
          color_map_boundary[3*index+1] = color_map[3*index+1] * obstruction_map[index];
          color_map_boundary[3*index+2] = color_map[3*index+2] * obstruction_map[index];
       }
     }
     return color_map_boundary;
}

void cfd::getVelocity(int i, int j, float &u, float &v){
   int index = i + nx*j;
   if(i >= 0 && i < nx && j >= 0 && j < ny){
      u = velocity_map[2*index+0]*obstruction_map[index];
      v = velocity_map[2*index+1]*obstruction_map[index];
   }else{
      u = 0.0;
      v = 0.0;
   }  
}

void cfd::getDensity(int i, int j, float &den){
   int index = i + nx*j;
   if(i >= 0 && i < nx && j >= 0 && j < ny){
      den = density_map[index]*obstruction_map[index];
   }else{
      den = 0.0;
   }   
}

void cfd::getColor(int i, int j, int channel, float &col){
   int index = i + nx*j;
   if(i >= 0 && i < nx && j >= 0 && j < ny){
      col = color_map[3*index+channel];
   }else{
      col = 0.0;
   } 
}

void cfd::getDivergence(int i, int j, float &div){
   int index = i + nx * j;
   if(i >= 0 && i < nx && j >= 0 && j < ny){
      div = divergence_map[index];
   }else{
      div = 0.0;
   }
}

void cfd::getPressure(int i, int j, float &pre){
   int index = i + nx * j;
   if(i >= 0 && i < nx && j >= 0 && j < ny){
      pre = pressure_map[index];
   }else{
      pre = 0.0;
   }
}

void cfd::interpolateVelocity(float x, float y, float &u, float &v){


   int i = int(x / dx);
   float wx = x - i * dx;
   int j = int(y / dy);
   float wy = y - j * dy;

   float u00, u10, u01, u11;
   float v00, v10, v01, v11;

      getVelocity(i+0, j+0, u00, v00);
      getVelocity(i+1, j+0, u10, v10);
      getVelocity(i+0, j+1, u01, v01);
      getVelocity(i+1, j+1, u11, v11);

      u =   u00 * (1 - wx) * (1 - wy) +
            u10 * wx       * (1 - wy) +
            u01 * (1 - wx) * wy       +
            u11 * wx       * wy;

      v =   v00 * (1 - wx) * (1 - wy) +
            v10 * wx       * (1 - wy) +
            v01 * (1 - wx) * wy       +
            v11 * wx       * wy; 
}

void cfd::interpolateDensity(float x, float y, float &res){

    int i = int(x / dx);
    int j = int(y / dy);
    float wx = x - i * dx;
    float wy = y - j * dy;

    float f00, f10, f01, f11;
    getDensity(i+0, j+0, f00);
    getDensity(i+1, j+0, f10);
    getDensity(i+0, j+1, f01);
    getDensity(i+1, j+1, f11);
    res = f00 * (1.0 - wx) * (1.0 - wy) +
          f10 * wx       * (1.0 - wy) +
          f01 * (1.0 - wx) * wy       +
          f11 * wx       * wy;   
}

void cfd::interpolateColor(float x, float y, int channel, float &res){

    int i = int(x / dx);
    int j = int(y / dy);
    float wx = x - i * dx;
    float wy = y - j * dy;

    float f00, f10, f01, f11;

    // Get the colors of four neighbors
    getColor(i+0, j+0, channel, f00);
    getColor(i+1, j+0, channel, f10);
    getColor(i+0, j+1, channel, f01);
    getColor(i+1, j+1, channel, f11);

    res = f00 * (1.0 - wx) * (1.0 - wy) +
            f10 * wx       * (1.0 - wy) +
            f01 * (1.0 - wx) * wy       +
            f11 * wx       * wy;
}

void cfd::advectDensity(){
   Initialize(temp_den_map, size, 0.0);
   int index;
   float u, v;
   float xp, yp, interpolatedDensity;

   for( int j=0;j<ny;j++ )
   {
       for(int i=0;i<nx;i++ )
       {
         index = i + nx*j;
         getVelocity(i,j,u,v);

         // back trace previous position (xp, yp)
         xp = i * dx - u * h;
         yp = j * dy - v * h;
         interpolateDensity(xp, yp, interpolatedDensity);
         temp_den_map[index] = interpolatedDensity; 
       }
   }
    swapDenMap();
}

void cfd::advectVelocity(){
   int index;
   float u, v;
   float up, vp; // interpolated velocity
   float xp, yp; // previous position 

   for( int j=0;j<ny;j++ )
   {
       for(int i=0;i<nx;i++ )
       {
         index = i + nx*j;
         getVelocity(i,j,u,v);

         // back trace previous position (xp, yp)
         xp = i * dx - u * h;
         yp = j * dy - v * h;
 
         interpolateVelocity(xp, yp, up, vp);
          
         temp_vel_map[2*index+0] = up;
         temp_vel_map[2*index+1] = vp;
       }
   }
   swapVelMap();
}

void cfd::advectColor(){

   int index;
   float u, v;
   float xp, yp;
   float interpolatedR, interpolatedG, interpolatedB;

   for( int j=0;j<ny;j++ )
   {
       for(int i=0;i<nx;i++ )
       {    
         index = i + nx*j;
         getVelocity(i,j,u,v);
         // back trace previous position (xp, yp)
         xp = i * dx - u * h;
         yp = j * dy - v * h;       

         interpolateColor(xp, yp, 0, interpolatedR);
         interpolateColor(xp, yp, 1, interpolatedG);
         interpolateColor(xp, yp, 2, interpolatedB);

         temp_col_map[3*index+0] = interpolatedR;
         temp_col_map[3*index+1] = interpolatedG;
         temp_col_map[3*index+2] = interpolatedB;

       }
   }
   swapColMap();
}

void cfd::bouyancy(){
   int index;
   float force;

   for( int j=0;j<ny;j++ )
   {
       for(int i=0;i<nx;i++ )
       {
         index = i + nx*j;
         force = density_map[index] * grav;
         velocity_map[2*index + 1] += force * h;    
       }
   }
}

void cfd::addForces(){
   bouyancy();
}

void cfd::drawSourcesToDensity(){

   for( int j=0;j<ny;j++ )
   {
       for(int i=0;i<nx;i++ )
       {
          int index = i + nx*j;
          density_map[index] += source_map[index];       
       }
   }
   Initialize(source_map, size, 0.0); // clear source_map
}

void cfd::modifyVelocity(){
   for(int j = 0; j < ny; j++){
      for(int i = 0; i < nx; i++){
         int index = i + nx * j;
         float p0, p1, p2, p3;
         getPressure(i+1, j, p0);
         getPressure(i-1, j, p1);
         getPressure(i,j+1,p2);
         getPressure(i,j-1,p3);
         
         temp_vel_map[2*index+0] = velocity_map[2*index+0]-(p0-p1)/(2.0*dx);
         temp_vel_map[2*index+1] = velocity_map[2*index+1]-(p2-p3)/(2.0*dy);
         
      }
   }

   swapVelMap();
}

void cfd::computeDivergence(){
   float u_up, v_up, u_down, v_down, u_left, v_left, u_right, v_right;
   int index;

   for(int j = 0; j < ny; j++){
      for(int i = 0; i < nx; i++){
         index = i + nx * j;
         getVelocity(i, j+1, u_up, v_up);
         getVelocity(i-1, j, u_left, v_left);
         getVelocity(i+1, j, u_right, v_right);
         getVelocity(i, j-1, u_down, v_down);
         divergence_map[index] = (u_right - u_left)/(2.0*dx) + (v_up - v_down)/(2.0*dy)+divergence_source_map[index];

      }
   }
}

void cfd::computePressure(){
      Initialize(pressure_map, size, 0.0);
      int index;
      float div, p0, p1, p2, p3;

      for(int t = 0; t<GS; t++){

         for( int j = 0; j < ny; j++){
            for(int i = 0; i < nx; i++){
               index = i + nx * j;
               
               getDivergence(i, j, div);
               getPressure(i+1,j,p0);
               getPressure(i-1,j,p1);
               getPressure(i,j+1,p2);
               getPressure(i,j-1,p3);

               pressure_map[index] = 0.25 * (p0 + p1 + p2 + p3)-0.25*div*dx*dx;
            }
         } 
      }     
}

void cfd::modifyBoundary(){
   for(int j=0; j<ny; j++){
      for(int i = 0; i < nx; i++){
         int index = i + nx * j;

         if (i == 0 || i == nx - 1){
            velocity_map[2*index+0] = 0.0;
         }

         if (j == 0 || j == ny - 1){
            velocity_map[2*index+1] = 0.0;
         } 
      }
   }
}

void cfd::project(){
   for(int i = 0; i < IOP; i++){
         computeDivergence();
         modifyBoundary();
         computePressure();
         modifyVelocity();
      }
      Initialize(divergence_source_map, size, 0.0);
      modifyBoundary();
}

void cfd::updateFluid(){
   advectDensity();
   advectVelocity();
   drawSourcesToDensity();
   addForces();
   project();
   advectColor();
}

void cfd::Initialize( float *data, int size, float value )
{
#pragma omp parallel for
   for(int i=0;i<size;i++ ) { data[i] = value; }
}

void cfd::swapColMap(){
    swap_col = temp_col_map;
    temp_col_map = color_map;
    color_map = swap_col;
}

void cfd::swapDenMap(){
    swap_den = temp_den_map;
    temp_den_map = density_map;
    density_map = swap_den;
  
}

void cfd::swapVelMap(){
    swap_vel = temp_vel_map;
    temp_vel_map = velocity_map;
    velocity_map = swap_vel;
}

cfd::~cfd(){
    delete color_map_boundary;
    delete density_map;
    delete color_map;
    delete velocity_map;
    delete pressure_map;
    delete divergence_map;
    delete temp_vel_map;
    delete temp_col_map;
    delete temp_den_map;
    delete temp_pre_map;
    delete swap_vel;
    delete swap_col;
    delete swap_den;
}