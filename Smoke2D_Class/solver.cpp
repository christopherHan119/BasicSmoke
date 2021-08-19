#pragma once

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "solver.h"

Solver::Solver()
{
}
 
Solver::~Solver()
{
    delete  velocityU;
    delete  velocityV;
       
    delete  previousVelocityU;
    delete  previousVelocityV;
       
    delete  density;
    delete  previousDensity;
}

void Solver::init(int n)
{
    N = n;

    timeStep = 0.1f;
    diffuseCoef = 0.0f;
    viscocity = 0.0f;
    force = 90.0f;
    source = 800.0f;
    addDensity = false;

    const int eachGridCount = N + 2;
    const int size = eachGridCount * eachGridCount;

    velocityU = (float*)malloc(size * sizeof(float));
    velocityV = (float*)malloc(size * sizeof(float));

    previousVelocityU = (float*)malloc(size * sizeof(float));
    previousVelocityV = (float*)malloc(size * sizeof(float));

    density = (float*)malloc(size * sizeof(float));
    previousDensity = (float*)malloc(size * sizeof(float));

    memset(velocityU, 0, sizeof(float) * size);
    memset(velocityV, 0, sizeof(float) * size);

    memset(previousVelocityU, 0, sizeof(float) * size);
    memset(previousVelocityV, 0, sizeof(float) * size);

    memset(density, 0, sizeof(float) * size);
    memset(previousDensity, 0, sizeof(float) * size);
}

void Solver::update()
{
    wuSource();
    get_density(density, previousDensity, velocityU, velocityV);
    get_velocity(velocityU, velocityV, previousVelocityU, previousVelocityV);
}

void Solver::addSource()
{
    addDensity = true;
}

void Solver::reset()
{
    int size = (N + 2) * (N + 2);

    memset(velocityU, 0, sizeof(float) * size);
    memset(velocityV, 0, sizeof(float) * size);

    memset(previousVelocityU, 0, sizeof(float) * size);
    memset(previousVelocityV, 0, sizeof(float) * size);

    memset(density, 0, sizeof(float) * size);
    memset(previousDensity, 0, sizeof(float) * size);
}

const float* Solver::getDensity()
{
    return density;
}

float Solver::getDensity(int x, int y)
{
    return density[wuIndex(x, y)];
}

const float* Solver::getVelocityU()
{
    return velocityU;
}

const float* Solver::getVelocityV()
{
    return velocityV;
}


float Solver::getVelocityU(int x, int y)
{
    return velocityU[wuIndex(x, y)];
}

float Solver::getVelocityV(int x, int y)
{
    return velocityV[wuIndex(x, y)];
}

void Solver::add_source(float* value, float* value_prev)
{
    int count;
    int size = (N + 2) * (N + 2);

    for (count = 0; count < size; count++)
    {
        value[count] = value[count] + timeStep * value_prev[count];
    }
}

void Solver::wuSource()
{
    float* density = previousDensity;
    float* u = previousVelocityU;
    float* v = previousVelocityV;

    int indexX, indexY;
    const int eachGridSize = N + 2;
    const int size = eachGridSize * eachGridSize;

    memset(u, 0, sizeof(float) * size);
    memset(v, 0, sizeof(float) * size);
    memset(density, 0, sizeof(float) * size);

    if (addDensity == true)
    {
        indexX = N / 2;
        indexY = N / 2;

        density[wuIndex(indexX, indexY)] = source;
        v[wuIndex(indexX, indexY)] = source;
        //u[wuIndex(indexX, indexY)] = source;
        addDensity = false;
    }

    // gravity source
    //for (int x = 1; x <= N; x++)
    //{
    //    for (int y = 1; y <= N; y++)
    //    {

    //        //if (density[wuIndex(x, y)] > 0)
    //        {
    //            v[wuIndex(x, y)] -= -1.f;
    //        }
    //    }
    //}
}

void Solver::boundary_condition(float* value, int flag)
{
    int i, j;

    for (i = 1; i <= N; i++) {
        value[wuIndex(0, i)] = flag == 1 ? -value[wuIndex(1, i)] : value[wuIndex(1, i)];
        value[wuIndex(N + 1, i)] = flag == 1 ? -value[wuIndex(N, i)] : value[wuIndex(N, i)];
        value[wuIndex(i, 0)] = flag == 2 ? -value[wuIndex(i, 1)] : value[wuIndex(i, 1)];
        value[wuIndex(i, N + 1)] = flag == 2 ? -value[wuIndex(i, N)] : value[wuIndex(i, N)];
    }

    value[wuIndex(0, 0)] = 0.5f * (value[wuIndex(1, 0)] + value[wuIndex(0, 1)]);
    value[wuIndex(0, N + 1)] = 0.5f * (value[wuIndex(1, N + 1)] + value[wuIndex(0, N)]);
    value[wuIndex(N + 1, 0)] = 0.5f * (value[wuIndex(N, 0)] + value[wuIndex(N + 1, 1)]);
    value[wuIndex(N + 1, N + 1)] = 0.5f * (value[wuIndex(N, N + 1)] + value[wuIndex(N + 1, N)]);
}

void Solver::lin_solve(int b, float* value, float* value_prev, float a, float c)
{
    int count;
    int count_x;
    int count_y;

    for (count = 0; count < 20; count++)
    {
        for (count_x = 1; count_x <= N; count_x++)
        {
            for (count_y = 1; count_y <= N; count_y++)
            {
                value[IX(count_x, count_y)] =
                    (value_prev[IX(count_x, count_y)] +
                        a * (value[IX(count_x - 1, count_y)] + value[IX(count_x + 1, count_y)] +
                            value[IX(count_x, count_y - 1)] + value[IX(count_x, count_y + 1)])) / c;
            }
        }

        boundary_condition(value, b);
    }
}

void Solver::diffuse(int b, float* value, float* value_prev)
{
    float alpha = timeStep * diffuseCoef * N * N * N;

    lin_solve(b, value, value_prev, alpha, 1 + 4 * alpha);
}

void Solver::advect(int b, float* d, float* d0, float* u, float* v)
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = timeStep * N;

    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
            if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; i0 = (int)x; i1 = i0 + 1;
            if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; j0 = (int)y; j1 = j0 + 1;

            s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;

            d[IX(i, j)] = 
                s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }

    boundary_condition(d, b);
}

void Solver::project(float* velocity_u, float* velocity_v, float* p, float* div)
{
    int count_x;
    int count_y;

    for (count_x = 1; count_x <= N; count_x++)
    {
        for (count_y = 1; count_y <= N; count_y++)
        {
            div[IX(count_x, count_y)] = (float)(-1.0 / 3.0 * (
                (velocity_u[IX(count_x + 1, count_y)] - velocity_u[IX(count_x - 1, count_y)]) / N +
                (velocity_v[IX(count_x, count_y + 1)] - velocity_v[IX(count_x, count_y - 1)]) / N ));
            p[IX(count_x, count_y)] = 0;
            
        }
    }

    boundary_condition(div, 0);
    boundary_condition(p, 0);

    lin_solve(0, p, div, 1, 4);

    for (count_x = 1; count_x <= N; count_x++)
    {
        for (count_y = 1; count_y <= N; count_y++)
        {
            velocity_u[IX(count_x, count_y)] -= 0.5f * N * (p[IX(count_x + 1, count_y)] - p[IX(count_x - 1, count_y)]);
            velocity_v[IX(count_x, count_y)] -= 0.5f * N * (p[IX(count_x, count_y + 1)] - p[IX(count_x, count_y - 1)]);

        }
    }

    boundary_condition(velocity_u, 1);
    boundary_condition(velocity_v, 2);
}

void Solver::get_density(float* density, float* density_prev, float* velocity_u, float* velocity_v)
{
    add_source(density, density_prev);

    SWAP(density_prev, density);
    diffuse(0, density, density_prev);

    SWAP(density_prev, density);
    advect(0, density, density_prev, velocity_u, velocity_v);
}

void Solver::get_velocity(float* velocity_u, float* velocity_v, float* velocity_u_prev, float* velocity_v_prev)
{
    add_source(velocity_u, velocity_u_prev);
    add_source(velocity_v, velocity_v_prev);

    SWAP(velocity_u_prev, velocity_u);
    SWAP(velocity_v_prev, velocity_v);

    diffuse(1, velocity_u, velocity_u_prev);
    diffuse(2, velocity_v, velocity_v_prev);

    project(velocity_u, velocity_v, velocity_u_prev, velocity_v_prev);


    SWAP(velocity_u_prev, velocity_u);
    SWAP(velocity_v_prev, velocity_v);

    advect(1, velocity_u, velocity_u_prev, velocity_u_prev, velocity_v_prev);
    advect(2, velocity_v, velocity_v_prev, velocity_u_prev, velocity_v_prev);

    project(velocity_u, velocity_v, velocity_u_prev, velocity_v_prev);
}

int Solver::wuIndex(int i, int j)
{
    const int eachGridSize = N + 2;
    int aa = i + eachGridSize * j;
    return aa;
}