#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"


void additional_forces_rad_r(struct reb_simulation* const r){
    // Drag-to-Keplerian force.
    double xi, yi, zi, ri;
    double vrx, vry, vrz, vri;
    double rt = 1.;
    double OK, at;
    struct reb_particle* const particles = r->particles;
    const int N = r->N;
    double tau = particles[1].hash;
    double tr = particles[2].hash;
    for (int i=0;i<N-1;i++){
        // position info
        xi = particles[i+1].x-particles[0].x;
        yi = particles[i+1].y-particles[0].y;
        zi = particles[i+1].z-particles[0].z;
        ri = sqrt(xi*xi+yi*yi+zi*zi);
        // find vr
        vrx = xi * (particles[i+1].vx - particles[0].vx) / ri;
        vry = yi * (particles[i+1].vy - particles[0].vy) / ri;
        vrz = zi * (particles[i+1].vz - particles[0].vz) / ri;
        vri = (vrx+vry+vrz);
        // find the traping force
        OK = sqrt( r->G * particles[0].m/ri/ri/ri);
        at = OK/tr * (ri - rt);
        // apply damping
        particles[i+1].ax -= vri * xi/ri / tau + at * xi/ri;
        particles[i+1].ay -= vri * yi/ri / tau + at * yi/ri;
        particles[i+1].az -= vri * zi/ri / tau + at * zi/ri;
        // heartbest
        // printf("%f, %f \n",vrx,particles[i+1].x);
    }
}



// Coordinate Note:
//
// r2 = x*x + y*y
// r-hat = [x,y,0]/r
// phi-hat = [-y,x,0]/r
// z-hat = [0,0,1]
//

void disk_drag(struct reb_simulation* const r){
    // Drag-to-Keplerian force.
    double xi, yi, zi, ri;
    double vxi, vyi, vzi;
    double vKx,vKy,vKz, vK;
    struct reb_particle* const pts = r->particles;
    const int N = r->N;
    double tau_d = pts[0].hash;
    for (int i=0;i<N-1;i++){
        // relative position info
        xi = pts[i+1].x-pts[0].x;
        yi = pts[i+1].y-pts[0].y;
        zi = pts[i+1].z-pts[0].z;
        ri = sqrt(xi*xi+yi*yi);
        // relative velocity info
        vxi = pts[i+1].vx-pts[0].vx;
        vyi = pts[i+1].vy-pts[0].vy;
        vzi = pts[i+1].vz-pts[0].vz;
        // find vK: phi-hat = [-y,x,0]/r
        vK = sqrt( r->G * pts[0].m/ri);
        vKx = vK * (-yi/ri);
        vKy = vK * (+xi/ri);
        vKz = 0;
        // apply damping
        pts[i+1].ax -= (vxi-vKx)/tau_d;
        pts[i+1].ay -= (vyi-vKy)/tau_d;
        pts[i+1].az -= (vzi-vKz)/tau_d;
    }
}


void disk_ecc_inc_damping(struct reb_simulation* const r){
    // Damp-to-circular force.
    double xi, yi, zi, ri;
    double vxi, vyi, vzi;
    double vrx, vry, vr;
    struct reb_particle* const pts = r->particles;
    const int N = r->N;
    double tau_e = pts[1].hash;
    double tau_i = pts[2].hash;
    for (int i=0;i<N-1;i++){
        // relative position info
        xi = pts[i+1].x-pts[0].x;
        yi = pts[i+1].y-pts[0].y;
        zi = pts[i+1].z-pts[0].z;
        ri = sqrt(xi*xi+yi*yi);
        // relative velocity info
        vxi = pts[i+1].vx-pts[0].vx;
        vyi = pts[i+1].vy-pts[0].vy;
        vzi = pts[i+1].vz-pts[0].vz;
        // find vr
        vr = (xi*vxi+yi*vyi)/ri;
        vrx = vr * (xi/ri);
        vry = vr * (yi/ri);
        // apply damping
        if (tau_e > 0.1){
            pts[i+1].ax -= vrx/tau_e*2.;
            pts[i+1].ay -= vry/tau_e*2.;
        }
        if (tau_i > 0.1){
            pts[i+1].az -= vzi/tau_i;
        }
    }
}



void disk_trap(struct reb_simulation* const r){
    // Migration-trap force.
    double xi, yi, zi, ri;
    double OK0;
    struct reb_particle* const pts = r->particles;
    const int N = r->N;
    double tau_t = pts[0].hash;
    double tau_i;
    double r0 = 1.;
    for (int i=0;i<N-1;i++){
        // relative position info
        xi = pts[i+1].x-pts[0].x;
        yi = pts[i+1].y-pts[0].y;
        zi = pts[i+1].z-pts[0].z;
        ri = sqrt(xi*xi+yi*yi);
        // find OK and tau
        OK0 = sqrt( r->G * pts[0].m/r0/r0/r0);
        tau_i = tau_t * (pts[1].m/pts[i+1].m);
        // apply damping
        if (tau_t > 0.1){
            pts[i+1].ax -= OK0 * (ri-r0) * (-yi/ri) / tau_i;
            pts[i+1].ay -= OK0 * (ri-r0) * (+xi/ri) / tau_i;
        }
    }
}



void disk_ecc_inc_trap(struct reb_simulation* const r){
    // Damp-to-circular force.
    double xi, yi, zi, ri;
    double vxi, vyi, vzi;
    double vrx, vry, vr;
    double OK0;
    struct reb_particle* const pts = r->particles;
    const int N = r->N;
    double tau_t = pts[0].hash;
    double tau_ti;
    double r0 = 1.;
    double tau_e = pts[1].hash;
    double tau_i = pts[2].hash;
    for (int i=0;i<N-1;i++){
        // relative position info
        xi = pts[i+1].x-pts[0].x;
        yi = pts[i+1].y-pts[0].y;
        zi = pts[i+1].z-pts[0].z;
        ri = sqrt(xi*xi+yi*yi);
        // relative velocity info
        vxi = pts[i+1].vx-pts[0].vx;
        vyi = pts[i+1].vy-pts[0].vy;
        vzi = pts[i+1].vz-pts[0].vz;
        // find vr
        vr = (xi*vxi+yi*vyi)/ri;
        vrx = vr * (xi/ri);
        vry = vr * (yi/ri);
        // apply damping
        if (tau_e > 0.1){
            pts[i+1].ax -= vrx/tau_e*2.;
            pts[i+1].ay -= vry/tau_e*2.;
        }
        if (tau_i > 0.1){
            pts[i+1].az -= vzi/tau_i;
        }
        // find OK and tau
        OK0 = sqrt( r->G * pts[0].m/r0/r0/r0);
        tau_ti = tau_t * (pts[1].m/pts[i+1].m);
        // apply damping
        if (tau_t > 0.1){
            pts[i+1].ax -= OK0 * (ri-r0) * (-yi/ri) / tau_ti;
            pts[i+1].ay -= OK0 * (ri-r0) * (+xi/ri) / tau_ti;
        }
    }
}




