import rebound
import numpy as np
import sys
import time
import os
import ctypes
import shutil
from itertools import combinations
from joblib import Parallel, delayed

###########################
##
##  Setup runs
##
###########################

def K2a(K,a1,M1,M2,Ms):
    '''
    Find a2 from a1 using orbital spacing in the unit of mutual Hill radius.
    Mutual Hill radius RH is defined as:
        RH :=  [(M1+M2)/(3Ms)]^(1/3) * [(a1+a2)/2]
    The outer planet:
        a2 = a1 + K*RH
    This function needs to solve for RH first by substitution.
    '''

    eta = ((M1+M2)/(3.*Ms))**(1./3.)    # mass ratio factor
    RH = (eta/(1.-0.5*K*eta))*a1        # Find RH by substitute a2 into RH
    a2 = a1 + K*RH

    return a2


def set_sim(Ms,rs,pt_init_list,R_CE,damp_info):
    '''
    Initialize a simluation with 1 star Ms  and multiple planets. 
    Note the default values of the parameters. Units are in AU-yr-Msun system.
    The ptn_init_list should have the formation:
    [ for each planet [m, a, K, e, r, inc, Omega, f, omega]]
    If K is positive, semi-major axis a will be overwritten.
    Hence, if using K, the planets must be listed in the order of increasing a.
    '''

    sim = rebound.Simulation()          # Create a REBOUND simulation
    sim.units = ('AU', 'yr', 'Msun')    # Units for planetary simulations

    # add the star
    sim.add(m=Ms,r=rs)

    # add the planets
    Npt = len(pt_init_list)
    for i in range(Npt):
        # read the initial condition input
        pt = pt_init_list[i]
        (m,a,e,r) = (pt[0],pt[1],pt[3],pt[4])
        (inc,Ome,f,ome) = (pt[5],pt[6],pt[7],pt[8])
        K = pt[2]
        # overwite a (by a1 and K) if a positive K is given
        if K > 0.:
            pt1 = pt_init_list[i-1]
            (m1,a1) = (pt1[0],pt1[1])
            a = K2a(K,a1,m1,m,Ms)
        # add this planet with the given initial condition
        sim.add(m=m,a=a,e=e,r=r,inc=inc,Omega=Ome,M=f,omega=ome)
    
    # move to the COM frame for better accuracy
    sim.move_to_com()

    # integration spec
    sim.integrator = "ias15" # Integrator handles close encounters.
    #sim.collision = "direct" 
    #sim.collision_resolve = "halt"
    sim.exit_min_distance = R_CE

    # disk force
    if disk_force_info[0] in not None:
        t_trap = disk_force_info[1]
        t_damp = disk_force_info[2]
        t_ecc = disk_force_info[3]
        t_inc = disk_force_info[4]
        if disk_force_info[0] == 'damp':
            sim.particles[1].hash(t_damp)
            sim.additional_forces = ctypes.CDLL('./disk_force.so').disk_damp_to_keplerian
        if disk_force_info[0] == 'trap':
            sim.particles[0].hash(t_trap)
            sim.additional_forces = ctypes.CDLL('./disk_force.so').disk_trap_zone
        if disk_force_info[0] = 'damp+trap':
            sim.particles[0].hash(t_trap)
            sim.particles[1].hash(t_damp)
            sim.additional_forces = ctypes.CDLL('./disk_force.so').disk_damp_and_trap
        sim.force_is_velocity_dependent = 1

    # all set
    return sim




###########################
##
##  Resolve collisions
##
###########################

def resolve_close(sim, R_check):
    '''(sim, cnt_ff, cnt_ism, cnt_ss) --> sim, cnt_ff, cnt_ism, cnt_ss
    Take a simulation that just reach the cloes-encounter boundary, zoom in and
    run to a point near the pericenter. Save the elements of the mutual binary orbits.
    '''

    # Get particles from simulation to apply our prescription
    sim.exit_min_distance = 0.
    #sim.automateSimulationArchive(sim_name,interval=dt_small,deletefile=False)

    # find the encounter pair

    d2_min = 20.
    pts = sim.particles

    t_in = sim.t

    for i1, i2 in combinations(range(sim.N),2):

        dpt = pts[i1] - pts[i2]
        d2 = dpt.x*dpt.x + dpt.y*dpt.y + dpt.z*dpt.z
        if d2<d2_min:
            d2_min = d2
            pt1 = pts[i1]
            pt2 = pts[i2]
            ipt1 = i1
            ipt2 = i2

    # remove the star crasher
    if pt1.m > 0.9:
        meat = pt2.m
        sim.remove(i2)
        sim.exit_min_distance = R_check
        return sim, sim.t, -1, pt1.m, meat, sim.t-t_in


    # integrate until until well separated

    d2_min = 20.
    Twait = 500.
    Nstep = int(Twait)*3000
    check_points = np.linspace(sim.t,sim.t+Twait,Nstep)

    for ck_pt in check_points:

        # shift to the BBH frame
        pcm = (pt1*pt1.m + pt2*pt2.m)/(pt1.m+pt2.m)
        for i in range(sim.N):
            pts[i].x = pts[i].x - pcm.x
            pts[i].y = pts[i].y - pcm.y
            pts[i].z = pts[i].z - pcm.z
            pts[i].vx = pts[i].vx - pcm.vx
            pts[i].vy = pts[i].vy - pcm.vy
            pts[i].vz = pts[i].vz - pcm.vz

        sim.integrate(ck_pt,exact_finish_time=1)

        for i1, i2 in combinations(range(sim.N),2):
            if (i1!=ipt1) or (i2!=ipt2):
                dpt = pts[i1]-pts[i2]
                d2 = dpt.x*dpt.x + dpt.y*dpt.y + dpt.z*dpt.z
                if d2 < R_check*R_check:
                    raise Exception('--------Two CEs at the same time!!!--------')

        oBBH = pt1.calculate_orbit(pt2)
        aBBH_neg = ( oBBH.a < 0) #and aBBH_neg

        dpt = pt1-pt2
        d2 = dpt.x*dpt.x + dpt.y*dpt.y + dpt.z*dpt.z

        if d2 < d2_min:
            d2_min = d2
            oBBH_rp = oBBH
            t_rp = sim.t

        if ( d2 > R_check*R_check and aBBH_neg ):
            sim.move_to_com()
            sim.exit_min_distance = R_check
            #sim.automateSimulationArchive(sim_name,interval=dt_large,deletefile=False)
            return sim, t_rp, pt1.m, pt2.m, oBBH_rp, sim.t-t_in

    # The code below should never be excuted.
    raise Exception('Resolving enconuter: cannot leave 4R. (%1.9f, %1.9f, %1.9f)'%(d2,oBBH.a,oBBH.e ))
    return None


###########################
##
##  Run and main
##
###########################

def run(sim,tout_list,run_id):
    '''(sim, list of time) --> tend, int, int, int
    Provided with a Rebound simulation and a list of output timestep,
    this function integrates the your system and writes output files.
    It also returns the final timestep and the outcome of integration.
    '''

    safe = 1        # Assume it is safe first, to be modified when necessary
    collision = 0   # Change it to one if merger happens
    ejection = 0    # Change it to one if ejection happens

    stop_flag = False
    ej_flag = False

    SC_list = np.zeros((0,2))   # 0: t,  1: mp
    CE_list = np.zeros((0,10))   # 0: t,  1: m1,  2: m2,  3: a,  4: e,  5: inc,  6: h,  7: v,  8: d, 9: dt
    Ej_list = np.zeros((0,2))   # 0: t,  1: mp

    Energy_list = np.zeros((len(tout_list),2+4*(sim.N-1))) # 0: t, 1: E, 2-5: m1, a1, e1, i1,  6-9: m2, a2, e2, i2

    for itt, tout in enumerate(tout_list):

        # Collisions will be resolved as mergers.
        # Integrate all the way the time of snapshot.
        if tout>sim.t:
            try:
                sim.integrate(tout,exact_finish_time=0)
            except rebound.Encounter as error:
                R_check = sim.exit_min_distance
                sim, CE_t, CE_m1, CE_m2, CE_res, CE_dt = resolve_close(sim, R_check)
                sim.exit_min_distance = R_check
                if CE_m1<0:
                    SC_list = np.concatenate(( SC_list , np.array([[CE_t,CE_res]]) ))
                else:
                    CE_add = np.array([[CE_t, CE_m1, CE_m2, CE_res.a, CE_res.e, CE_res.inc, CE_res.h, CE_res.v, CE_res.d, CE_dt]])
                    CE_list = np.concatenate(( CE_list, CE_add))

            save_add = np.zeros(2+4*(sim.N-1))
            save_add[0] = sim.t
            save_add[1] = sim.calculate_energy()
            for i in range(sim.N-1):
                pt = sim.particles[i+1]
                orbi = pt.calculate_orbit()
                save_add[2+i*4:6+i*4] = np.array([pt.m,orbi.a,orbi.e,orbi.inc])
            Energy_list[itt,0:len(save_add)] = save_add[:]

            # get ready to eject particles
            sim.move_to_com()
            d2max = 1000*1000
            # find all particles that need to be ejected
            d2list = np.zeros(sim.N)
            for i in range(sim.N):
                pt = sim.particles[i]
                d2list[i] =  pt.x*pt.x + pt.y*pt.y + pt.z*pt.z
            L2Ej = np.arange(sim.N)[d2list>d2max]
            # eject from the last particle
            for j in L2Ej[::-1]:
                ptej = sim.particles[j]
                np.savetxt('./i-data/Ej%05d-%02d.dat'%(run_id,j), np.array([ptej.x,ptej.y,ptej.vx,ptej.vy]) )
                Ej_list = np.concatenate(( Ej_list, [[sim.t, ptej.m]] ))
                sim.remove(j)
                sim.move_to_com()

        if sim.N<2.5:
            break

    #np.save('./i-data/Energy%05d.npy'%run_id,Energy_list)

    return sim, SC_list, CE_list, Ej_list



## Main program here
def set_and_run(run_id,K,a1):

    # time span of the simulation
    Torb = a1**1.5
    Norb = 1.e+5
    time_list = np.linspace(0.,1.,100001)
    time_list = time_list * Norb * Torb

    # stellar parameters
    Ms = 1.
    rs = 1e-8
    mu = 1e-5
    R_CE = a1*(2*mu/Ms/3.)**(1./3.)

    # planets
    np.random.seed()
    tp = 2.*np.pi 
    pt_init_list = [
            #  m   a    K    e     r      inc         Omega                     f                      omega
            [2*mu, 1., None, 0., rs*mu, R_CE/a1, tp*np.random.rand(), tp*np.pi*np.random.rand(), tp*np.pi*np.random.rand()],
            [mu,  1.3, 2,  1e-5, rs*mu, R_CE/a1, tp*np.random.rand(), tp*np.pi*np.random.rand(), tp*np.pi*np.random.rand()]
            ]

    # create a simulaiton
    sim_name = "sim%03d_J.bin" % (run_id+1)
    sim = set_sim(Ms,rs,pt_init_list,R_CE,damp_info)
    #sim.automateSimulationArchive(sim_name,interval=1e-2,deletefile=True)
    sim.simulationarchive_snapshot("./snap/sim%05d_t0.bin" % (run_id+1))

    # perform the simulation
    print('Simulation %03d starts...' % (run_id+1) )
    wall_time_start = time.time()
    sim_fin, SC_data, CE_data, Ej_data = run(sim,time_list,run_id)
    sim_fin.simulationarchive_snapshot("./snap/sim%05d_tf.bin" % (run_id+1))
    wall_time_end = time.time()
    print('Simulatoin %03d has finished. Wall time elapse: ' + str(wall_time_end-wall_time_start))

    # store the output
    NSC = len(SC_data[:,0])
    NCE = len(CE_data[:,0])
    NEj = len(Ej_data[:,0])
    #np.savetxt('./i-data/SC%05d.dat' %run_id, SC_data)
    #np.savetxt('./i-data/CE%05d.dat' %run_id, CE_data)
    #np.savetxt('./i-data/Ej%05d.dat' %run_id, Ej_data)
    #np.save('./i-data/SC%05d.npy' %run_id, SC_data)
    #np.save('./i-data/CE%05d.npy' %run_id, CE_data)
    #np.save('./i-data/Ej%05d.npy' %run_id, Ej_data)

    # all set
    return SC_data, CE_data, Ej_data, NSC, NCE, NEj



#########################
##
##    Parallel run
##
#########################

if __name__ == '__main__':

    # number of simulations tp run
    num_core = 50
    num_runs = 2000
    K = 20.e-1
    a1 = 1.

    # prepare the output folders
    shutil.rmtree('./i-data',ignore_errors=True)
    shutil.rmtree('./snap',ignore_errors=True)
    os.mkdir('./i-data')
    os.mkdir('./snap')

    # start the simulations in parallel
    #res_full = np.zeros((num_runs,4+5*3))
    res_list = Parallel(n_jobs=num_core)(delayed(set_and_run)(i,K,a1) for i in range(num_runs))

    # count the results
    (NSC,NEj,NCE) = (0,0,0)
    for i in range(num_runs):
        NSC = NSC + res_list[i][3]
        NCE = NCE + res_list[i][4]
        NEj = NEj + res_list[i][5]

    # store the reults
    CE_res = np.zeros((NCE,10))
    #SC_res = np.zeros((NSC,2))
    #Ej_res = np.zeros((Nej,2))
    (SCi,CEi,Eji) = (0,0,0)
    for i in range(num_runs):
        nCE = res_list[i][4]
        CE_res[CEi:CEi+nCE,:] = res_list[i][1]
        CEi = nCE + CEi
    np.savetxt('./CE_all.dat',CE_res)
    np.save('./CE_all.npy',CE_res)

    # print a summary message
    print('===========')
    print('--Results--')
    print('K0=%2.1f'%K)
    print('Close Enc: %9d (%12.2f%%)' %( NCE, (NCE*1.)/(num_runs*1.)*100) )
    print('Ejections: %9d (%12.2f%%)' %( NEj, (NEj*1.)/(num_runs*1.)*100) )
    print('Star feed: %9d (%12.2f%%)' %( NSC, (NSC*1.)/(num_runs*1.)*100) )
    print('-----------')
    print('===========')


