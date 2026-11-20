import astro
import simulacions

if __name__ == "__main__":
    print("running from main")

    H = 3600
 
    R_hist = astro.RK4(h=H)
    s_hist, l_hist = astro.canviMarcSolTerra(R_hist, h=H)

    print("done")

    simulacions.simMarcTerra(s_hist, l_hist, "marc terra", speed=3600*10, h=H)