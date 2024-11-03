import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

from r3 import coords, rs
from r3 import xpp, xpn, x0p, x0n, xnp, xnn, ypp, ypn, y0p, y0n, ynp, ynn, zpp, zpn, z0p, z0n, znp, znn

# prerequisites : sudo apt install ffmpeg

# face colors
# (x+, x-, y+, y-, z+, z-) : (blue, green, red, orange, yellow, white)
hex = {
    (0, 1): "#3F51B5",
    (0, -1): "#4CAF50",
    (1, 1): "#B51F1F",
    (1, -1): "#FF6F00",
    (2, 1): "#FFEB3B",
    (2, -1): "#FFFFFF"
}

colors = np.array([hex[(coords[i][-1], coords[i][coords[i][-1]])] for i in range(len(coords))])

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

def draw(ax, coord, color):
    # 2 * 2 grid
    grid = [None, None, None]
    n = coord[-1] # normal vector
    grid[n] = np.ones((2,2)) * (coord[n] + 0.5) if coord[n] > 0 else np.ones((2,2)) * (coord[n] - 0.5)
    i, j = sorted({0,1,2} - {n})
    grid[i], grid[j] = np.meshgrid(np.linspace(coord[i] - 0.5, coord[i] + 0.5, 2), np.linspace(coord[j] - 0.5, coord[j] + 0.5, 2))
    x, y, z = grid
    ax.plot_surface(x, y, z, color=color, edgecolor="black", shade=False)

def draw_cube(ax, coords):
    for coord, color in zip(coords, colors):
        draw(ax, coord, color)

def draw_axes(ax, arrow_offset=3, text_offset=3.5):
    # draw arrows
    x0, y0, z0 = np.zeros((3,3))
    x1, y1, z1 = np.array([[1,0,0], [0,1,0], [0,0,1]]) * arrow_offset
    ax.quiver(x0, y0, z0, x1, y1, z1, arrow_length_ratio=0.1, color="black")
    # annotate arrows
    for i, a in enumerate(["x", "y", "z"]):
        args = [0,0,0,0]
        args[i] = text_offset
        args[-1] = a
        ax.text(*args, color="black")

def draw_rotation(ax, r, radius=2.5, theta_bgn=0, theta_end=350):
    if r is None:
        return
    # init
    theta = np.linspace((2 * np.pi) * theta_bgn / 360, (2 * np.pi) * theta_end / 360 , theta_end - theta_bgn)
    if r.orient == -1:
        theta = -theta
    points = np.zeros((3, len(theta)))
    # assign value
    orient = [[1,2,0], [2,0,1], [0,1,2]]
    i, j, k = orient[r.axis]
    points[k] = np.ones(len(theta)) * r.level
    points[i] = radius * np.cos(theta)
    points[j] = radius * np.sin(theta)
    x, y, z = points
    # draw rotation curve
    ax.plot(x[:-1], y[:-1], z[:-1], color="red")
    # draw arrowhead
    ax.quiver(x[-2], y[-2], z[-2], x[-1]-x[-2], y[-1]-y[-2], z[-1]-z[-2], arrow_length_ratio=5, color="red")
    # label rotation
    ax.text2D(0.05, 0.95, str(r), color='red', fontsize='large', transform=ax.transAxes)

def animate(i, seq, coords_list, angle_per_frame, limit=2):
    angle = i * angle_per_frame
    seq += [None]
    ax = plt.gca()
    if angle % 360 == 0:
        ax.cla()
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_zlim(-limit, limit)
        plt.axis("off")
        r = seq[angle // 360]
        coords = coords_list[angle // 360]
        draw_cube(ax, coords)
        draw_axes(ax)
        draw_rotation(ax, r)

    azim = angle % 360
    elev = (angle * 4) % 360
    if 90 < elev <= 270:
        elev = 180 - elev
    elif 270 < elev < 360:
        elev -= 360
    elev /= 6
    ax.view_init(elev=elev, azim=azim)
    return fig,

if __name__ == "__main__":
    seq = [x0n, x0p, ypn, z0n] * 1
    coords_list = [coords]
    for r in seq:
        coords_list.append(r(coords_list[-1]))

    angle_per_frame = 6
    frames = (360 // angle_per_frame) * (len(seq) + 1)
    interval = 1
    fps = 12

    anim = animation.FuncAnimation(fig, animate, frames=frames, interval=interval, blit=True, fargs=(seq, coords_list, angle_per_frame))
    # mp4 format
    writer = animation.FFMpegWriter(fps=fps)
    anim.save("rubik.mp4", writer=writer)
