import taichi as ti
import bc_init
import part_dist_init
import sim
import output
import numpy as np
import matplotlib.pyplot as plt

######################
### INITIALIZATION ###
######################


# Initialize Taichi GPU backend
ti.init(arch=ti.gpu)

# D2Q9 Arrays
simD = 2
simQ = 9
w = ti.static([4/9, 1/9, 1/9, 1/9, 1/9,
               1/36, 1/36, 1/36, 1/36])
c = ti.static([
    [ 0, 1, 0, -1,  0, 1, -1, -1,  1],
    [ 0, 0, 1,  0, -1, 1,  1, -1, -1]
])
cbar_indx = ti.static([0, 3, 4, 1, 2, 7, 8, 5, 6])


# Simulation Parameters
dx = 1
dt = 1
tau = 0.91*dt
output_fp = ".\\2D_fluidic_osc.mat"

# Boundary Conditions
img_fp = "./boundary_conditions/imgs/final_grids/2D_fluidic_oscillator.png"
Nx, Ny, bc_output = bc_init.read_bc_img(img_fp)
bc_output, bc_SWB_conditions, bc_IOPB_conditions, bc_IOUB_conditions = bc_init.identify_boundaries(Nx,Ny,bc_output,c,simQ,cbar_indx)


# Initialize Particle Distribution
f_mat, r_mat, u_mat = part_dist_init.init_sim(Nx,Ny,simD,simQ,w)


######################
##### SIMULATION #####
######################
r_out, u_out, u_mat_frames = sim.run(dx, dt, tau, simQ, w, c, cbar_indx, Nx, Ny, f_mat, r_mat, u_mat, bc_output, bc_SWB_conditions, bc_IOPB_conditions, bc_IOUB_conditions)


######################
####### OUTPUT #######
######################
output.plot_flow_quantities(output_fp, r_out, u_out)

output.generate_gif(u_mat_frames)

