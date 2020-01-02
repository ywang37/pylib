from matplotlib.colors import ListedColormap
import numpy as np
import os

dn = os.path.dirname(os.path.realpath(__file__))

WhGrYlRd_scheme = np.genfromtxt(dn+'/WhGrYlRd.txt', delimiter=' ')
WhGrYlRd_map = ListedColormap(WhGrYlRd_scheme/255.0)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    try:
        from viscm import viscm
        viscm(WhGrYlRd_map)
    except ImportError:
        print("viscm not found, falling back on simple display")
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',
                   cmap=WhGrYlRd_map)
    plt.show()
