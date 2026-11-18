import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as ani

# Posició inicial
R = np.array([
    1.4710e11,       # r (m)
    -24.0,           # dr (m/s)
    np.deg2rad(113.05),  # theta (rad)
    -1.57e-8,        # dtheta (rad/s)
    0.0,             # phi (rad)
    -7.27e-5         # dphi (rad/s)
])

Omega = 7.292115e-5
Mt    = 5.9722e24
Ms    = 1.98847e30
G     = 6.67430e-11

year   = 3600*24*365


def fdr(R):
    r = R[0]
    theta = R[2]
    dtheta = R[3]
    dphi = R[5]
    return r*dtheta**2+r*np.sin(theta)**2*(dphi+Omega)**2-G*(Mt+Ms)/r**2

def fdtheta(R):
    r = R[0]
    dr = R[1]
    theta = R[2]
    dtheta = R[3]
    dphi = R[5]
    return np.sin(2*theta)/2*(dphi+Omega)**2-2*dtheta*dr/r

def fdphi(R):
    r = R[0]
    dr = R[1]
    theta = R[2]
    dtheta = R[3]
    dphi = R[5]
    return (-2*dr/r-np.sin(2*theta)/np.sin(theta)**2*dtheta)*(dphi+Omega)

def fr(R):
    return R[1]

def ftheta(R):
    return R[3]

def fphi(R):
    return R[5]

def F(R):
    return np.array([
        fr(R),
        fdr(R),
        ftheta(R),
        fdtheta(R),
        fphi(R),
        fdphi(R)
    ])


def step_RK4(R, t, F, h):
    K1 = F(R)
    R2 = R.copy() + h/2*K1
    K2 = F(R2)
    R3 = R.copy() + h/2*K2
    K3 = F(R3)
    R4 = R.copy() + h*K3
    K4 = F(R4)

    return R + h/6 * (K1 + 2*K2 + 2*K3 + K4)

def RK4(R=R, h=3600, t_final=year, F=F):
    R_hist = [R]
    for t in np.linspace(0, t_final, int(t_final/h)):

        R = step_RK4(R, t, F, h)

        R_hist.append(R)
    
    return R_hist

def sim(R_hist, name, h=3600, speed=3600*24, fps=10, dpi=100):
    R_hist = np.array(R_hist)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    r_max = np.max(R_hist[:, 0])
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    ax.set_zlim(-r_max, r_max)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    sun = ax.scatter([0], [0], [0])
    ax.scatter([0], [0], [0], color="blue")

    stride = max(1, int(round(speed/fps/h)))
    indexs = np.arange(0, len(R_hist), stride)

    def update(frame):
        R = R_hist[indexs[frame]]

        r = R[0]
        theta = R[2]
        phi = R[4]

        x = r*np.cos(phi)*np.sin(theta)
        y = r*np.sin(phi)*np.sin(theta)
        z = r*np.cos(theta)

        sun._offsets3d = ([x], [y], [z])

    animation = ani.FuncAnimation(fig=fig, func=update, frames=len(indexs), interval=1000/fps,
    blit=False)
    animation.save(f"multimedia/{name}.mp4", writer="ffmpeg", dpi=dpi, fps=fps)
    plt.close(fig)

        


if __name__ == "__main__":
    print("Running from main")

    R_hist = RK4()

    print("done")

    sim(R_hist, "test")

    