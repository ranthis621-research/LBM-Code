import taichi as ti


def init_sim(Nx, Ny, simD, simQ, w):
    f = ti.field(ti.f32, shape=(Ny,Nx,simQ))
    r_mat = ti.field(ti.f32, shape=(Ny,Nx))
    u_mat = ti.Vector.field(simD, ti.f32, shape=(Ny,Nx))

    @ti.kernel
    def generate_flow_properties():
        # Initialize particle distribution
        for y, x in r_mat:
            r_mat[y, x] = 1
            u_mat[y, x] = [0, 0]

            # Specify initial distribution using density distribution
            for q in ti.static(range(simQ)):
                f[y, x, q] = w[q] * r_mat[y, x]


    generate_flow_properties()
    return f, r_mat, u_mat