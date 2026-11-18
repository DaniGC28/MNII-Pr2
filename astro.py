import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as ani


# Constants

G  = 6.67430e-11     
Mt = 5.9722e24       
Ml = 7.342e22          
Ms = 1.98841e30       

hour = 3600
day  = hour*24
year = day*365

eps    = np.radians(23.44)
eix_terra = eix = np.array([0, np.sin(eps), np.cos(eps)])
omega  = 7.2921159e-5

# valors inicials
# Posicions respecte al Sol [m] i velocitats [m/s]
# Marc: eclíptica J2000 centrada al Sol, 2027-01-01 00:00 TDB
a  = np.array([-2.6073e10,  1.44770e11, -8.93e6])   # Sol -> Terra
da = np.array([-2.9790e4,  -5.397e3,     0.29])

b  = np.array([-2.5930e10,  1.45100e11,  2.280e7])  # Sol -> Lluna
db = np.array([-3.0792e4,  -4.978e3,     2.79])

R = np.array([a, da, b, db])

m0 = np.array([-0.1850, 0.9016, -0.3909])           # centre terra -> lat=0 lon=0
theta0 = np.radians(100.66)



def fa(R):
    return R[1]

def fb(R):
    return R[3]

def fda(R):
    a = R[0]
    b = R[2]
    La = np.linalg.norm(a)
    Lb = np.linalg.norm(b)
    Lab= np.linalg.norm(a-b)
    return -G*( (Ms+Mt)*a/La**3 + Ml*( (a-b)/Lab**3 + b/Lb**3 ) )

def fdb(R):
    a = R[0]
    b = R[2]
    La = np.linalg.norm(a)
    Lb = np.linalg.norm(b)
    Lab= np.linalg.norm(a-b)
    return -G*( (Ms+Ml)*b/Lb**3 + Mt*( (b-a)/Lab**3 + a/La**3 ) )

def F(R):
    return np.array([
        fa(R),
        fda(R),
        fb(R),
        fdb(R)
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

def meridia(t):
    th = theta0 + omega*t
    x_eq, y_eq = np.cos(th), np.sin(th)
    return np.array([x_eq,
                    y_eq*np.cos(eps),
                    -y_eq*np.sin(eps)])

def sim(R_hist, name, h=3600, speed=3600*24, fps=10, dpi=100):

    R_hist = np.array(R_hist)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    r_max = np.max(np.abs(R_hist[:, [0, 2]]))
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    ax.set_zlim(-r_max, r_max)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    sun   = ax.scatter([0], [0], [0], color="yellow")
    earth = ax.scatter([0], [0], [0], color="blue")
    #moon = ax.scatter([0], [0], [0], color="grey")

    L1 = 0.15 * r_max
    axis_line, = ax.plot([], [], [], color="green", lw=2)

    L2 = 0.15 * r_max
    facing_line, = ax.plot([], [], [], color="purple", lw=2)

    stride = max(1, int(round(speed/fps/h)))
    indexs = np.arange(0, len(R_hist), stride)


    def update(frame):
        R = R_hist[indexs[frame]]
        t = h*indexs[frame]

        a = R[0]
        b = R[2]

        xt = a[0]
        yt = a[1]
        zt = a[2]

        xl = b[0]
        yl = b[1]
        zl = b[2]

        earth._offsets3d = ([xt], [yt], [zt])
        #moon._offsets3d = ([xl], [yl], [zl])

        axis1 = a - L1*eix
        axis2 = a + L1*eix
        axis_line.set_data_3d([axis1[0], axis2[0]], [axis1[1], axis2[1]], [axis1[2], axis2[2]])


        m = meridia(t)

        facing1 = a
        facing2 = a + L2*m
        facing_line.set_data_3d([facing1[0], facing2[0]], [facing1[1], facing2[1]], [facing1[2], facing2[2]])

    animation = ani.FuncAnimation(fig=fig, func=update, frames=len(indexs), interval=1000/fps,
    blit=False)
    animation.save(f"multimedia/{name}.mp4", writer="ffmpeg", dpi=dpi, fps=fps)
    plt.close(fig)

        


if __name__ == "__main__":
    print("Running from main")

    R_hist = RK4()

    print("done")

    sim(R_hist, "test", speed=3600*10)

    