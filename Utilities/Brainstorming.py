import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import LinearLocator

def get_intensity(beta: float, Resolution: int = 30):

    theta = np.linspace(0,     np.pi, Resolution)
    phi   = np.linspace(0, 2 * np.pi, Resolution)

    phi, theta = np.meshgrid(phi,theta)

    Intensity = 1 / (1 - beta * np.cos(theta))**3 * (1 - np.sin(theta)**2 * np.cos(phi)**2 * (1 - beta**2) / (1 - beta * np.cos(theta))**2)

    Intensity = Intensity / np.max(Intensity)

    return Intensity, theta, phi

# Make data.

fig = plt.figure()

for index, beta in enumerate([0.4]):

    Intensity, theta, phi = get_intensity(beta = beta, Resolution = 100)

    X = Intensity * np.sin(theta) * np.cos(phi)
    Y = Intensity * np.sin(theta) * np.sin(phi)
    Z = Intensity * np.cos(theta)

    ax  = fig.add_subplot(1, 1 , index + 1, projection='3d') 

    colors = plt.cm.jet( Intensity / np.max(Intensity) )
    surf = ax.plot_surface(X, Y, Z, 
                        rstride = 1, 
                        cstride = 1, 
                        facecolors = colors,
                        linewidth = 0, 
                        antialiased = True, 
                        alpha = 0.5, 
                        zorder = 0.5)

    Acceleration_Vector = np.array([0, 0, 0, 1.5, 0, 0])
    Velocity_Vector     = np.array([0, 0, 0, 0, 0, 1.5])

    ax.quiver(*Acceleration_Vector, arrow_length_ratio = 0.1, colors = "k")
    ax.quiver(*Velocity_Vector, arrow_length_ratio = 0.1, colors = "k")

    ax.text(0, 0, 1.5,  r"$\,\,{\vec{\beta}}$", size=15, zorder=1, color='k') 
    ax.text(1.5, 0, 0,  r"$\,\dot{\vec{\beta}}$", size=15, zorder=1, color='k') 

    ax.set_box_aspect((1,1,1))
    ax.grid(False)
    ax.axis('on')

    ax.set_title(r'$\beta = $ {}'.format(beta))

    ax.set_xlabel(r'$X$')
    ax.set_ylabel(r'$Y$')
    ax.set_zlabel(r'$Z$')

    ax.set_xlim([-1.2,1.2])
    ax.set_ylim([-1.2,1.2])
    ax.set_zlim([-1.2,1.2])

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
plt.show()


plt.show()