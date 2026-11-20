from astro import meridia
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

refugi_lat = np.radians(42.48121924621716)
refugi_lon = np.radians(1.7455967387680447)

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

def simMarcSol(R_hist, name, h=3600, speed=3600*24, fps=10, dpi=100):

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

def simMarcTerra(s_hist, l_hist, name, h=3600, speed=3600*24, fps=10, dpi=100, elev=25, azim=30):

    s_hist = np.array(s_hist)
    l_hist = np.array(l_hist)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    ax.set_box_aspect((1, 1, 1))
    ax.set_proj_type("ortho")
    ax.view_init(elev=elev, azim=azim)

    r_max = np.max(np.abs(l_hist))
    ax.set_xlim(-r_max, r_max)
    ax.set_ylim(-r_max, r_max)
    ax.set_zlim(-r_max, r_max)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    L = 0.4 * r_max
    r = 0.25 * r_max
    

    #sun   = ax.scatter([0], [0], [0], color="yellow")
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    earth  = ax.plot_surface(r*np.cos(u)*np.sin(v),
                    r*np.sin(u)*np.sin(v),
                    r*np.cos(v),
                    color="blue", alpha=0.8, linewidth=0)
    moon = ax.scatter([0], [0], [0], color="grey")
    refugi = ax.scatter([r*np.cos(refugi_lon)*np.cos(refugi_lat)], 
                        [r*np.sin(refugi_lon)*np.cos(refugi_lat)],
                        [r*np.sin(refugi_lat)], color="orange", edgecolor="k")

    axis_line, = ax.plot([], [], [], color="green", lw=2)

    sun_line, = ax.plot([], [], [], color="gold", lw=2)

    phi = np.linspace(0, 2*np.pi, 200)
    Rr = r * 1.02   # lleugerament més gran que la Terra perquè no quedi tapat

    # Equador (lat = 0): cercle al pla XY
    ax.plot(Rr*np.cos(phi), Rr*np.sin(phi), 0*phi, color="red", lw=1.5)

    # Meridià de Greenwich (lon = 0): cercle al pla XZ
    ax.plot(Rr*np.cos(phi), 0*phi, Rr*np.sin(phi), color="cyan", lw=1.5)

    stride = max(1, int(round(speed/fps/h)))
    indexs = np.arange(0, len(s_hist), stride)


    def update(frame):
        s = s_hist[indexs[frame]]
        l = l_hist[indexs[frame]]
        t = h*indexs[frame]

        xs = s[0]
        ys = s[1]
        zs = s[2]

        xl = l[0]
        yl = l[1]
        zl = l[2]

        moon._offsets3d = ([xl], [yl], [zl])

        axis1 = -L*np.array([0, 0, 1])
        axis2 = L*np.array([0, 0, 1])
        axis_line.set_data_3d([axis1[0], axis2[0]], [axis1[1], axis2[1]], [axis1[2], axis2[2]])

        sun1 = r*s/np.linalg.norm(s)
        sun2 = L*s/np.linalg.norm(s)*10
        sun_line.set_data_3d([sun1[0], sun2[0]], [sun1[1], sun2[1]], [sun1[2], sun2[2]])

    animation = ani.FuncAnimation(fig=fig, func=update, frames=len(indexs), interval=1000/fps,
    blit=False)
    animation.save(f"multimedia/{name}.mp4", writer="ffmpeg", dpi=dpi, fps=fps)
    plt.close(fig)