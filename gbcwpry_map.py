from matplotlib.colors import LinearSegmentedColormap

cm_data = [[0.0, 1.0, 0.0],
           [0.0, 0.5, 0.5],
           [0.0, 0.0, 1.0],
           [0.0, 0.5, 1.0],
           [0.0, 1.0, 1.0],
           [1.0, 1.0, 1.0],
           [1.0, 0.0, 1.0],
           [1.0, 0.0, 0.5],
           [1.0, 0.0, 0.0],
           [1.0, 0.5, 0.0],
           [1.0, 1.0, 0.0]]

gbcwpry_map = LinearSegmentedColormap.from_list('gbcwpry', cm_data)
# For use of "viscm view"
test_cm = gbcwpry_map

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    try:
        from viscm import viscm
        viscm(gbcwpry_map)
    except ImportError:
        print("viscm not found, falling back on simple display")
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',
                   cmap=gbcwpry_map)
    plt.show()
