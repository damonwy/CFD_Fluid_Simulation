#include <iostream>
#include <algorithm>
#include <cmath>
#include "cfd.h"

cfd::cfd(float _lx, float _ly, int _nx, int _ny, int _GS, int _IOP, float _h, float _grav, Scheme _advect){
    nx = _nx;
    ny = _ny;
    lx = _lx;
    ly = _ly;
    GS = _GS;
    IOP = _IOP;
    h = _h;
    grav = _grav;
    advection_scheme = _advect;
    has_vorticity = false;
    has_surfacetension = false;
    V = 1.0;
    T = 0.1;

    dx = lx / (nx - 1);
    dy = ly / (ny - 1);
    size = nx * ny;

    source_map = new float[size];
    Initialize(source_map, size, 0.0);

    density_map = new float[size];
    density_map_f = new float[size];
    density_map_b = new float[size];
    density_map_err = new float[size];
    density_map_bfe = new float[size];
    temp_density_map_b = new float[size];
    temp_density_map_f = new float[size];
    Initialize(density_map, size, 0.0);
    Initialize(density_map_f, size, 0.0);
    Initialize(density_map_b, size, 0.0);
    Initialize(temp_density_map_b, size, 0.0);
    Initialize(temp_density_map_f, size, 0.0);

    velocity_map = new float[2*size];
    velocity_map_b = new float[2*size];
    velocity_map_f = new float[2*size];
    velocity_map_err = new float[2*size];
    velocity_map_bfe = new float[2*size];
    temp_velocity_map_b = new float[2*size];
    temp_velocity_map_f = new float[2*size];
    Initialize(velocity_map, 2*size, 0.0);
    Initialize(velocity_map_f, 2*size, 0.0);
    Initialize(velocity_map_b, 2*size, 0.0);
    Initialize(temp_velocity_map_f, 2*size, 0.0);
    Initialize(temp_velocity_map_b, 2*size, 0.0);
    Initialize(velocity_map_bfe, 2*size, 0.0);


    color_map = new float[3*size];
    color_map_b = new float[3*size];
    color_map_f = new float[3*size];
    color_map_err = new float[3*size];
    color_map_bfe = new float[3*size];
    temp_color_map_b = new float[3*size];
    temp_color_map_f = new float[3*size];
    color_map_boundary = new float[3*size];
    Initialize(color_map, 3*size, 0.0);
    Initialize(color_map_f, 3*size, 0.0);
    Initialize(color_map_b, 3*size, 0.0);
    Initialize(temp_color_map_f, 3*size, 0.0);
    Initialize(temp_color_map_b, 3*size, 0.0);
    Initialize(color_map_bfe, 3*size, 0.0);

    w_map = new float[size];
    Initialize(w_map, size, 0.0);
    obstruction_map = new float[size];
    Initialize(obstruction_map, size, 1.0);

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
    Initialize(temp_vel_map, 2*size, 0.0);
    Initialize(temp_col_map, 3*size, 0.0);
    Initialize(temp_den_map, size, 0.0);
    Initialize(temp_pre_map, size, 0.0);
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

void cfd::getVelocity(int i, int j, float &u, float &v, float *velocityMap){
    int index = i + nx*j;
    if(i >= 0 && i < nx && j >= 0 && j < ny){
      u = velocityMap[2*index+0]*obstruction_map[index];
      v = velocityMap[2*index+1]*obstruction_map[index];
    }else{
      u = 0.0;
      v = 0.0;
    }  
}

void cfd::getAreaRotation(int i, int j, float &w){
    int index = i + nx * j;
    if(i >= 0 && i < nx && j >= 0 && j < ny){
        w = w_map[index];
    }else{
        w = 0.0;
    }
}

void cfd::getDensity(int i, int j, float &den, float* densityMap){

   int index = i + nx*j;

   if(i >= 0 && i < nx && j >= 0 && j < ny && den){

      den = densityMap[index]*obstruction_map[index];
       if(den < 0.0001){
         den = 0.0;
         densityMap[index] = 0.0;
       }

   }else{
      den = 0.0;
   }   
}

void cfd::getColor(int i, int j, int channel, float &col, float* colorMap){
   int index = i + nx*j;
   if(i >= 0 && i < nx && j >= 0 && j < ny){
      col = colorMap[3*index+channel];
   }else{
      col = 1.0;
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

void cfd::interpolateVelocity(float x, float y, float &u, float &v, float *velocityMap){

   int i = int(x / dx);
   int j = int(y / dy);
   float wx = x - i * dx;
   float wy = y - j * dy;

   float u00, u10, u01, u11;
   float v00, v10, v01, v11;

      getVelocity(i+0, j+0, u00, v00, velocityMap);
      getVelocity(i+1, j+0, u10, v10, velocityMap);
      getVelocity(i+0, j+1, u01, v01, velocityMap);
      getVelocity(i+1, j+1, u11, v11, velocityMap);

      u =   u00 * (1 - wx) * (1 - wy) +
            u10 * wx       * (1 - wy) +
            u01 * (1 - wx) * wy       +
            u11 * wx       * wy;

      v =   v00 * (1 - wx) * (1 - wy) +
            v10 * wx       * (1 - wy) +
            v01 * (1 - wx) * wy       +
            v11 * wx       * wy; 
}

void cfd::interpolateDensity(float x, float y, float &res, float *densityMap){

    int i = int(x / dx);
    int j = int(y / dy);
    float wx = x - i * dx;
    float wy = y - j * dy;

    float f00, f10, f01, f11;
    getDensity(i+0, j+0, f00, densityMap);
    getDensity(i+1, j+0, f10, densityMap);
    getDensity(i+0, j+1, f01, densityMap);
    getDensity(i+1, j+1, f11, densityMap);
    res = f00 * (1.0 - wx) * (1.0 - wy) +
          f10 * wx       * (1.0 - wy) +
          f01 * (1.0 - wx) * wy       +
          f11 * wx       * wy;   
}

void cfd::interpolateColor(float x, float y, int channel, float &res, float *colorMap){

    int i = int(x / dx);
    int j = int(y / dy);
    float wx = x - i * dx;
    float wy = y - j * dy;

    float f00, f10, f01, f11;

    // Get the colors of four neighbors
    getColor(i+0, j+0, channel, f00, colorMap);
    getColor(i+1, j+0, channel, f10, colorMap);
    getColor(i+0, j+1, channel, f01, colorMap);
    getColor(i+1, j+1, channel, f11, colorMap);

    res = f00 * (1.0 - wx) * (1.0 - wy) +
            f10 * wx       * (1.0 - wy) +
            f01 * (1.0 - wx) * wy       +
            f11 * wx       * wy;
}

void cfd::advectDensity(float *&target_density_map, float *&source_density_map, float *&temp_density_map, bool isBackWards){
   Initialize(temp_density_map, size, 0.0);
   int index;
   float u, v;
   float xp, yp, interpolatedDensity;

   if(isBackWards == false){
    // SL advection 
      for( int j=0;j<ny;j++ )
         {
             for(int i=0;i<nx;i++ )
             {
               index = i + nx*j;
               getVelocity(i,j,u,v, velocity_map);

               // back trace previous position (xp, yp)
               xp = i * dx - u * h;
               yp = j * dy - v * h;
               interpolateDensity(xp, yp, interpolatedDensity, source_density_map);
               temp_density_map[index] = interpolatedDensity; 
             }
         }
   }else{
    // SL advection backwards
      for( int j=0;j<ny;j++ )
         {
             for(int i=0;i<nx;i++ )
             {
               index = i + nx*j;
               getVelocity(i,j,u,v, velocity_map);

               // forward trace previous position (xp, yp)
               xp = i * dx + u * h;
               yp = j * dy + v * h;
               interpolateDensity(xp, yp, interpolatedDensity, source_density_map);
               temp_density_map[index] = interpolatedDensity; 
             }
         }
   }

   for( int j=0;j<ny;j++ )
       {
           for(int i=0;i<nx;i++ )
           {    
             index = i + nx * j;
             target_density_map[2*index+0] = temp_density_map[2*index+0];
             target_density_map[2*index+1] = temp_density_map[2*index+1];
             
           }
       }
   
   //std::swap(target_density_map, temp_density_map);
}

void cfd::advectDensityScheme(Scheme s){
   
    //pf = p0;
   int index;
   // 0. density at initial frame
    for(int j=0; j < ny; j++){
      for(int i = 0; i <  nx; i++){
        index = i + nx * j;
        density_map_f[index] = density_map[index];
      }
    } 
    // 1. SL advect
    advectDensity(density_map_f, density_map_f, temp_density_map_f, false);

    // 2. SL advection backwards
    advectDensity(density_map_b, density_map_f, temp_density_map_b, true);

    // 3. Compute error 
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        index = i+nx*j;
        density_map_err[index] = 0.5*(density_map[index] - density_map_b[index]);
      }
    }

    if(s == BFECC){
      // 4. Add error
   
        for(int j=0; j<ny; j++){
          for(int i=0; i<nx; i++){
            index = i+nx*j;
            density_map_bfe[index] = density_map[index] + density_map_err[index];
          }
        }

        // 5. Update density
        advectDensity(density_map, density_map_bfe, temp_den_map, false);
    }

    if(s == MM){
      // 4. Add error
      for(int j=0; j<ny; j++){
        for(int i = 0; i<nx; i++){
          index = i+nx*j;
          density_map[index] = density_map_f[index]+density_map_err[index];
        }
      }
    }
}

void cfd::advectVelocityScheme(Scheme s){
   
    //pf = p0;
   int index;
   // 0. velocity at initial frame
    for(int j=0; j < ny; j++){
      for(int i = 0; i <  nx; i++){
        index = i + nx * j;
        velocity_map_f[2*index+0] = velocity_map[2*index+0];
        velocity_map_f[2*index+1] = velocity_map[2*index+1];
      }
    } 
    // 1. SL advect
    advectVelocity(velocity_map_f, velocity_map_f, temp_velocity_map_f, false);

    // 2. SL advection backwards
    advectVelocity(velocity_map_b, velocity_map_f, temp_velocity_map_b, true);

    // 3. Compute error 
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        index = i+nx*j;
        velocity_map_err[2*index+0] = 0.5*(velocity_map[2*index+0] - velocity_map_b[2*index+0]);
        velocity_map_err[2*index+1] = 0.5*(velocity_map[2*index+1] - velocity_map_b[2*index+1]);

      }
    }

    if(s == BFECC){
      // 4. Add error
        for(int j=0; j<ny; j++){
          for(int i=0; i<nx; i++){
            index = i+nx*j;
            velocity_map_bfe[2*index+0] = velocity_map[2*index+0] + velocity_map_err[2*index+0];
            velocity_map_bfe[2*index+1] = velocity_map[2*index+1] + velocity_map_err[2*index+1];
          }
        }

        // 5. Update velocity
        advectVelocity(velocity_map, velocity_map_bfe, temp_vel_map, false);
    }

    if(s == MM){
      // 4. Add error
      for(int j=0; j<ny; j++){
        for(int i = 0; i<nx; i++){
          index = i+nx*j;
          velocity_map[2*index+0] = velocity_map_f[2*index+0] + velocity_map_err[2*index+0];
          velocity_map[2*index+1] = velocity_map_f[2*index+1] + velocity_map_err[2*index+1];
        }
      }
    } 
}


void cfd::advectVelocity(float *&target_velocity_map, float *&source_velocity_map, float *&temp_velocity_map, bool isBackWards){

   int index;
   float u, v;
   float up, vp; // interpolated velocity
   float xp, yp; // previous position 

   Initialize(temp_velocity_map, 2*size, 0.0);

   if(isBackWards == false){
    //SL advection
        for( int j=0;j<ny;j++ )
     {
         for(int i=0;i<nx;i++ )
         {
           index = i + nx*j;
           getVelocity(i,j,u,v, velocity_map);

           // back trace previous position (xp, yp)
           xp = i * dx - u * h;
           yp = j * dy - v * h;
   
           interpolateVelocity(xp, yp, up, vp, source_velocity_map);
            
           temp_velocity_map[2*index+0] = up;
           temp_velocity_map[2*index+1] = vp;
         }
     }
   }else{
    //SL advection backwards
      for( int j=0;j<ny;j++ )
     {
         for(int i=0;i<nx;i++ )
         {
           index = i + nx*j;
           getVelocity(i,j,u,v, velocity_map);

           // back trace previous position (xp, yp)
           xp = i * dx + u * h;
           yp = j * dy + v * h;
   
           interpolateVelocity(xp, yp, up, vp, source_velocity_map);
            
           temp_velocity_map[2*index+0] = up;
           temp_velocity_map[2*index+1] = vp;
         }
     }
   }

   for( int j=0;j<ny;j++ )
       {
           for(int i=0;i<nx;i++ )
           {    
             index = i + nx * j;
             target_velocity_map[2*index+0] = temp_velocity_map[2*index+0];
             target_velocity_map[2*index+1] = temp_velocity_map[2*index+1];
             
           }
       }
 
   //std::swap(target_velocity_map, temp_velocity_map);
}

void cfd::advectColor(float *&target_color_map, float *&source_color_map, float *&temp_color_map, bool isBackWards){
   Initialize(temp_color_map, 3*size, 0.0);
   int index;
   float u, v;
   float xp, yp;
   float interpolatedR, interpolatedG, interpolatedB;
    if(isBackWards == false){
        //SL advection
       for( int j=0;j<ny;j++ )
       {
           for(int i=0;i<nx;i++ )
           {    
             index = i + nx*j;
             getVelocity(i,j,u,v, velocity_map);
             // back trace previous position (xp, yp)
             xp = i * dx - u * h;
             yp = j * dy - v * h;       

             interpolateColor(xp, yp, 0, interpolatedR, source_color_map);
             interpolateColor(xp, yp, 1, interpolatedG, source_color_map);
             interpolateColor(xp, yp, 2, interpolatedB, source_color_map);

             temp_color_map[3*index+0] = interpolatedR;
             temp_color_map[3*index+1] = interpolatedG;
             temp_color_map[3*index+2] = interpolatedB;

           }
       }
     }else{
        for( int j=0;j<ny;j++ )
       {
           for(int i=0;i<nx;i++ )
           {    
             index = i + nx*j;
             getVelocity(i,j,u,v, velocity_map);
             // back trace previous position (xp, yp)
             xp = i * dx + u * h;
             yp = j * dy + v * h;       

             interpolateColor(xp, yp, 0, interpolatedR, source_color_map);
             interpolateColor(xp, yp, 1, interpolatedG, source_color_map);
             interpolateColor(xp, yp, 2, interpolatedB, source_color_map);

             temp_color_map[3*index+0] = interpolatedR;
             temp_color_map[3*index+1] = interpolatedG;
             temp_color_map[3*index+2] = interpolatedB;

           }
       }
        
     }

     for( int j=0;j<ny;j++ )
       {
           for(int i=0;i<nx;i++ )
           {    
             index = i + nx * j;
             target_color_map[3*index+0] = temp_color_map[3*index+0];
             target_color_map[3*index+1] = temp_color_map[3*index+1];
             target_color_map[3*index+2] = temp_color_map[3*index+2];

           }
       }
 // std::swap(target_color_map, temp_color_map);
}

void cfd::advectColorScheme(Scheme s){
   
    //pf = p0;
   int index;
   // 0. color at initial frame
    for(int j=0; j < ny; j++){
      for(int i = 0; i <  nx; i++){
        index = i + nx * j;
        color_map_f[3*index+0] = color_map[3*index+0];
        color_map_f[3*index+1] = color_map[3*index+1];
        color_map_f[3*index+2] = color_map[3*index+2];
      }
    } 
    // 1. SL advect
    advectColor(color_map_f, color_map_f, temp_color_map_f, false);

    // 2. SL advection backwards
    advectColor(color_map_b, color_map_f, temp_color_map_b, true);

    // 3. Compute error 
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        index = i+nx*j;
        color_map_err[3*index+0] = 0.5*(color_map[3*index+0] - color_map_b[3*index+0]);
        color_map_err[3*index+1] = 0.5*(color_map[3*index+1] - color_map_b[3*index+1]);
        color_map_err[3*index+2] = 0.5*(color_map[3*index+2] - color_map_b[3*index+2]);

      }
    }

    if(s == BFECC){
            // 4. Add error
          for(int j=0; j<ny; j++){
            for(int i=0; i<nx; i++){
              index = i+nx*j;
              color_map_bfe[3*index+0] = color_map[3*index+0] + color_map_err[3*index+0];
              color_map_bfe[3*index+1] = color_map[3*index+1] + color_map_err[3*index+1];
              color_map_bfe[3*index+2] = color_map[3*index+2] + color_map_err[3*index+2];
            }
          }

          // 5. Update velocity
          advectColor(color_map, color_map_bfe, temp_col_map, false);
     }

    if(s == MM){
        // 4. Add error
        for(int j=0; j<ny; j++){
          for(int i = 0; i<nx; i++){
            index = i+nx*j;
            color_map[3*index+0] = color_map_f[3*index+0] + color_map_err[3*index+0];
            color_map[3*index+1] = color_map_f[3*index+1] + color_map_err[3*index+1];
            color_map[3*index+2] = color_map_f[3*index+2] + color_map_err[3*index+2];
          }
        }
    }

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

void cfd::surface_tension(){
  int index;
  float pr, pl, pu, pd, p, pru, pld, prd, plu;
  float gradx, grady, normalizer, mxx, myy, mxy, k;

  for(int j=0; j<ny; j++){
    for(int i = 0; i<nx; i++){
      index = i + nx * j;

      getDensity(i+1, j+0, pr, density_map);
      getDensity(i-1, j+0, pl, density_map);
      getDensity(i+0, j+1, pu, density_map);
      getDensity(i+0, j-1, pd, density_map);
      getDensity(i+0, j+0, p,  density_map);
      getDensity(i+1, j+1, pru, density_map);
      getDensity(i-1, j-1, pld, density_map);
      getDensity(i+1, j-1, prd, density_map);
      getDensity(i-1, j+1, plu, density_map);

      gradx = (pr - pl)/(2*dx);
      grady = (pu - pd)/(2*dy);

      normalizer = sqrt(gradx*gradx + grady*grady);
      mxx = (pr - pl - 2 * p)/(dx * dx);
      myy = (pu - pd - 2 * p)/(dy * dy);
      mxy = (pru + pld - prd - plu)/(4*dx* dy);
      k = (mxx + myy - gradx * gradx * mxx - grady * grady * myy + 2 * gradx *grady * mxy);

      if(normalizer > 0.001){
        k /= normalizer;
        gradx /= normalizer;
        grady /= normalizer;
      }
      
      velocity_map[2*index + 0] += k * h * gradx * T;
      velocity_map[2*index + 1] += k * h * grady * T;    

    }
  }
}

void cfd::vorticity_confinement(){
    int index;
    float w, u_up, v_up,u_down, v_down, u_left, v_left, u_right, v_right;
    float ex, ey, normalizer;
    float e_right, e_left, e_up, e_down, e_ij, D;
    Initialize(w_map, size, 0.0);

    // Identify the area of rotation

    for(int j=0; j<ny; j++){
        for(int i = 0; i<nx; i++){
            index = i + nx * j;
            getVelocity(i+1, j+0, u_right,v_right, velocity_map);
            getVelocity(i-1, j+0, u_left,v_left, velocity_map);
            getVelocity(i+0, j+1, u_up,v_up, velocity_map);
            getVelocity(i+0, j-1, u_down,v_down, velocity_map);
            w = -(v_right - v_left)/(2*dx)+(u_up-u_down)/(2*dy);
            w_map[index]=std::abs(w);
        }
    }
    
    for(int j=0; j<ny; j++){
        for(int i = 0; i<nx; i++){
            index = i + nx * j;
            getAreaRotation(i+1, j+0, e_right);
            getAreaRotation(i-1, j+0, e_left);
            getAreaRotation(i+0, j+1, e_up);
            getAreaRotation(i+0, j-1, e_down);
            getAreaRotation(i+0, j+0, e_ij);
            getDivergence(i+0, j+0, D);

            ex = (e_right - e_left)/ (2 * dx);
            ey = (e_up - e_down)/(2*dy);
            normalizer = sqrt(ex * ex + ey * ey);

            if(normalizer > 0.0001){
                ex = ex/normalizer;
                ey = ey/normalizer;
            }
            
            velocity_map[2*index + 0] += e_ij * h * D * dx *(-ey) * V;
            velocity_map[2*index + 1] += e_ij * h * D * dx * ex * V;    
        }
    }
}

void cfd::updateScheme(Scheme s){
    advection_scheme = s;
}

void cfd::addForces(){
   bouyancy();
   if(has_surfacetension){
      surface_tension();
   }
   if(has_vorticity){
      vorticity_confinement();
   }
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

   for(int j = 0; j < ny; j++){
      for(int i = 0; i < nx; i++){
         int index = i + nx * j;
         velocity_map[2*index+0] = temp_vel_map[2*index+0];
         velocity_map[2*index+1] = temp_vel_map[2*index+1];
   }
 }
   //std::swap(velocity_map, temp_vel_map);
}

void cfd::computeDivergence(){
   float u_up, v_up, u_down, v_down, u_left, v_left, u_right, v_right;
   int index;

   for(int j = 0; j < ny; j++){
      for(int i = 0; i < nx; i++){
         index = i + nx * j;
         getVelocity(i, j+1, u_up, v_up, velocity_map);
         getVelocity(i-1, j, u_left, v_left, velocity_map);
         getVelocity(i+1, j, u_right, v_right, velocity_map);
         getVelocity(i, j-1, u_down, v_down, velocity_map);
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
         
         modifyBoundary();computeDivergence();
         computePressure();
         modifyVelocity();
      }
      Initialize(divergence_source_map, size, 0.0);
      //modifyBoundary();
}

void cfd::updateFluid(){
   if(advection_scheme == SL){
      advectDensity(density_map, density_map, temp_den_map, false);
      advectVelocity(velocity_map, velocity_map, temp_vel_map, false);
      drawSourcesToDensity();
      addForces();
      project();
      advectColor(color_map, color_map, temp_col_map, false);
   }

   if(advection_scheme != SL){
      advectDensityScheme(advection_scheme);
      advectVelocityScheme(advection_scheme);
      drawSourcesToDensity();
      addForces();
      project();
      advectColor(color_map, color_map, temp_col_map, false);
   }
   
}

void cfd::Initialize( float *data, int size, float value )
{
#pragma omp parallel for
   for(int i=0;i<size;i++ ) { data[i] = value; }
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

}