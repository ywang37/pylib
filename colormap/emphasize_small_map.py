from matplotlib.colors import LinearSegmentedColormap

cm_data = [[1.0, 1.0, 1.0],
           [0.0, 1.0, 1.0],
           [0.0, 0.5, 1.0],
           [0.0, 1.0, 0.0],
           [1.0, 1.0, 0.0],
           [1.0, 0.0, 0.0], 
           [0.85, 0.0, 0.35],
           [0.7, 0.0, 0.7],
           [0.63, 0.0, 0.63],
           [0.57, 0.0, 0.57],
           [0.5, 0.0, 0.5],
           [0.43, 0.0, 0.43],
           [0.37, 0.0, 0.37],
           [0.3, 0.0, 0.3]]

emphasize_small_map = LinearSegmentedColormap.from_list('emphasize_small', cm_data)
# For use of "viscm view"
test_cm = emphasize_small_map

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    try:
        from viscm import viscm
        viscm(emphasize_small_map)
    except ImportError:
        print("viscm not found, falling back on simple display")
        plt.imshow(np.linspace(0, 100, 256)[None, :], aspect='auto',
                   cmap=emphasize_small_map)
    plt.show()
