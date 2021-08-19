#pragma once


#include <stdio.h>

class Solver
{
public : 
    Solver();
    ~Solver();
        
    void init(int n);
    void update();
    void addSource();
    void reset();

    const float* getDensity();
    float getDensity(int x, int y);

    const float* getVelocityU();
    const float* getVelocityV();
    float getVelocityU(int x, int y);
    float getVelocityV(int x, int y);
 
 private:
 
    void add_source(float* value, float* value_prev);

    void wuSource();
    
    void boundary_condition(float* value, int flag);
    
    void lin_solve(int b, float* value, float* value_prev, float a, float c);
    
    void diffuse(int b, float* value, float* value_prev);

    void advect(int b, float* density, float* density_prev, float* velocity_u, float* velocity_v);

    void project(float* velocity_u, float* velocity_v, float* p, float* div);

    void get_density(float* density, float* density_prev, float* velocity_u, float* velocity_v);

    void get_velocity(float* velocity_u, float* velocity_v, float* velocity_u_prev, float* velocity_v_prev);

    int wuIndex(int i, int j);

private :

    int N;

   float timeStep;
   float diffuseCoef;
   float viscocity;
   float force;
   float source;

   float* velocityU;
   float* velocityV;
   //float* velocityW;

   float* previousVelocityU;
   float* previousVelocityV;
   //float* previousVelocityW;

   float* density;
   float* previousDensity;

   bool addDensity;
};