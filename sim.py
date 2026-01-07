import taichi as ti
import matplotlib.pyplot as plt
import numpy as np


def run(dx, dt, tau, simQ, w, c, cbar_indx, Nx, Ny, f_mat, r_mat, u_mat, bc_output, bc_SWB, bc_IOPB, bc_IOUB):
    # Initialization
    inv_cs2 = 3 * pow(dt / dx, 2)
    dp = 0.125
    p_in = 0.6
    p_out = p_in - dp
    uw_in = ti.static([-0.01,0])
    uw_out = ti.static([0.01,0])
    rho_in = p_in * inv_cs2
    rho_out = p_out * inv_cs2
    omega = dt / tau
    omega_p = (1 - omega)
    fp_mat = ti.field(ti.f32, shape=(Ny,Nx,simQ))

    # Simulation
    @ti.kernel
    def simulation_loop():
        for y, x in r_mat:
            r = r_mat[y,x]
            u = u_mat[y,x]

            for q in ti.static(range(simQ)):
                cx = c[0][q]
                cy = c[1][q]
                cu_dot = u.x*cx+u.y*cy
                uu_dot = u.x*u.x+u.y*u.y

                # Calculate equilibrium
                feq = w[q] * r * (1 + inv_cs2*cu_dot
                                  + 0.5*inv_cs2*inv_cs2*cu_dot*cu_dot
                                  - 0.5*inv_cs2*uu_dot)

                # Calculate collision
                f_mat[y,x,q] = f_mat[y,x,q]*omega_p + feq*omega

                # Propagate
                # Standard cell
                if bc_output[y,x] == 0:
                    fp_mat[y+cy, x+cx, q] = f_mat[y, x, q]

                # Boundary Cells
                elif bc_output[y,x] == 11:
                    # Solid Walls
                    if bc_SWB[y,x,q] == 0:
                        fp_mat[y+cy, x+cx, q] = f_mat[y, x, q]
                    else:
                        fp_mat[y, x, cbar_indx[q]] = f_mat[y, x, q]

                elif bc_output[y,x] == 2:
                    # Periodic
                    xp = (x + cx + Nx) % Nx
                    yp = (y + cy + Ny) % Ny

                    fp_mat[yp, xp, q] = f_mat[y, x, q]

                elif bc_output[y,x] == 12:
                    # Solid Walls - Periodic Hybrid
                    xp = (x + cx + Nx) % Nx
                    yp = (y + cy + Ny) % Ny
                    if bc_SWB[y,x,q] == 0:
                        fp_mat[yp, xp, q] = f_mat[y, x, q]
                    else:
                        fp_mat[y, x, cbar_indx[q]] = f_mat[y, x, q]

                elif bc_output[y,x] == 33:
                    # Inlet/Outlet BCs
                    if x < Nx//2:
                        # Left boundary - pressure inlet
                        if bc_IOPB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            uw = u_mat[y, x + 1]
                            cuw_dot = uw.x * cx + uw.y * cy
                            uwuw_dot = uw.x * uw.x + uw.y * uw.y
                            fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2*rho_in*w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)
                    else:
                        # Right boundary - pressure outlet
                        if bc_IOPB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            uw = u_mat[y, x - 1]
                            cuw_dot = uw.x * cx + uw.y * cy
                            uwuw_dot = uw.x * uw.x + uw.y * uw.y
                            fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2*rho_out*w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)

                elif bc_output[y,x] == 13:
                    # Wall - Inlet/Outlet hybrid BCs
                    if x < Nx // 2:
                        # Left boundary - pressure inlet
                        if bc_IOPB[y, x, q] == 0 and bc_SWB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            if bc_IOPB[y, x, q] != 0 or bc_SWB[y, x, q] != 0:
                                uw = u_mat[y, x+1]
                                cuw_dot = uw.x * cx + uw.y * cy
                                uwuw_dot = uw.x * uw.x + uw.y * uw.y
                                fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2 * rho_in * w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)
                    else:
                        # Right boundary - pressure outlet
                        if bc_IOPB[y, x, q] == 0 and bc_SWB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            if bc_IOPB[y, x, q] != 0 or bc_SWB[y, x, q] != 0:
                                uw = u_mat[y, x - 1]
                                cuw_dot = uw.x * cx + uw.y * cy
                                uwuw_dot = uw.x * uw.x + uw.y * uw.y
                                fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2 * rho_out * w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)

                elif bc_output[y,x] == 44:
                    # Inlet/Outlet BCs
                    if x < Nx//2:
                        # Left boundary - velocity inlet
                        if bc_IOUB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            cuw_dot = uw_in[0] * cx + uw_in[1] * cy
                            uwuw_dot = uw_in[0] * uw_in[0] + uw_in[1] * uw_in[1]
                            rho_uin = (f_mat[y, x, 0] + f_mat[y, x, 2] + f_mat[y, x, 4] + 2 * (f_mat[y, x, 3] + f_mat[y, x, 6] + f_mat[y, x, 7]))/(1-uw_in[0])
                            fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2*rho_uin*w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)
                    else:
                        # Right boundary - velocity outlet
                        if bc_IOUB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            cuw_dot = uw_out[0] * cx + uw_out[1] * cy
                            uwuw_dot = uw_out[0] * uw_out[0] + uw_out[1] * uw_out[1]
                            rho_uout = (f_mat[y, x, 0] + f_mat[y, x, 2] + f_mat[y, x, 4] + 2 * (f_mat[y, x, 1] + f_mat[y, x, 8] + f_mat[y, x, 5]))/(1-uw_out[0])
                            fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2*rho_uout*w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)

                elif bc_output[y,x] == 14:
                    # Wall - Inlet/Outlet hybrid BCs
                    if x < Nx // 2:
                        # Left boundary - velocity inlet
                        if bc_IOUB[y, x, q] == 0 and bc_SWB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            if bc_IOUB[y, x, q] != 0 or bc_SWB[y, x, q] != 0:
                                cuw_dot = uw_in[0] * cx + uw_in[1] * cy
                                uwuw_dot = uw_in[0] * uw_in[0] + uw_in[1] * uw_in[1]
                                rho_uin = r_mat[y,x+1]
                                fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2 * rho_uin * w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)
                    else:
                        # Right boundary - velocity outlet
                        if bc_IOUB[y, x, q] == 0 and bc_SWB[y, x, q] == 0:
                            fp_mat[y + cy, x + cx, q] = f_mat[y, x, q]
                        else:
                            if bc_IOUB[y, x, q] != 0 or bc_SWB[y, x, q] != 0:
                                cuw_dot = uw_out[0] * cx + uw_out[1] * cy
                                uwuw_dot = uw_out[0] * uw_out[0] + uw_out[1] * uw_out[1]
                                rho_uout = r_mat[y,x-1]
                                fp_mat[y, x, cbar_indx[q]] = -f_mat[y, x, q] + 2 * rho_uout * w[cbar_indx[q]] * (1 + inv_cs2*cuw_dot
                                  + 0.5*inv_cs2*inv_cs2*cuw_dot*cuw_dot
                                  - 0.5*inv_cs2*uwuw_dot)

        for y, x in r_mat:
            r_update = 0.0
            rux_update = 0.0
            ruy_update = 0.0
            for q in ti.static(range(simQ)):
                # Update distribution
                f_mat[y,x,q] = fp_mat[y,x,q]

                cx = c[0][q]
                cy = c[1][q]

                r_update += f_mat[y,x,q]
                rux_update += f_mat[y,x,q]*cx
                ruy_update += f_mat[y,x,q]*cy

            # Moment update
            r_mat[y,x] = r_update
            u_mat[y,x][0] = rux_update / r_update
            u_mat[y,x][1] = ruy_update / r_update

    u_mat_frames = np.zeros(shape=(Ny, Nx, 150),dtype=np.float32)
    n = 0
    steps = 60000
    for step in range(steps):
        simulation_loop()

        if step % 100 == 0:
            print(100*step/steps)

        if step > steps//2 and step % 300 == 0:
            u_mat_temp = u_mat.to_numpy()
            u_mat_frames[:,:,n] = u_mat_temp[:,:,0]
            n = n + 1


    return r_mat, u_mat, u_mat_frames