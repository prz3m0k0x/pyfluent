import numpy as np
from numpy import random as rnd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

POP_SIZE = 10
MAX_ITER = 30
w_init = 0.8
c1_init = 0.6
c2_init = 0.6
w_finish = 0.5
c1_finish = 0.9
c2_finish = 0.9

N_PARAMS = 4
N_RESPONSES = 3
objective_weights = np.array([0.5, 0.25, 0.25])

def function(particles):
    z1 = -(particles[0,:]**2 + 0.5*particles[1,:]**2)
    z2 = -(2*particles[0,:] + particles[1,:]**2)
    return np.vstack((z1, z2))



def populate(size):
  x1,x2 = -10, 3 #x1, x2 = right and left boundaries of our X axis
  pop = rnd.uniform(x1,x2, size) # size = amount of particles in population
  return pop

particles = np.vstack((populate(POP_SIZE), populate(POP_SIZE)))
best_positions = np.copy(particles)
velocity = np.zeros_like(particles) #velocity of each of the particle
obj_funct = function(particles)
gains = np.sum((np.multiply(obj_funct, objective_weights[:, np.newaxis])), axis=0)
positions = np.copy(particles) 
swarm_best_position = particles[:,np.argmax(gains)]
swarm_best_gain = np.max(gains) #highest gain

l = np.zeros((N_PARAMS, POP_SIZE, MAX_ITER))

algorithm_params = np.linspace((w_init, c1_init, c2_init), (w_finish, c1_finish, c2_finish), MAX_ITER)

x_lo, x_hi = -10, 3
pbest_positions = np.copy(particles)
pbest_gains     = gains.copy()
best_gain_history = []

for i in range(MAX_ITER):

    l[:, :, i] = np.copy(particles)

    w, c1, c2 = algorithm_params[i]
    r1 = rnd.uniform(0, 1, (N_PARAMS, POP_SIZE))
    r2 = rnd.uniform(0, 1, (N_PARAMS, POP_SIZE))

    velocity = (w  * velocity
              + c1 * r1 * (pbest_positions - particles)
              + c2 * r2 * (swarm_best_position[:, np.newaxis] - particles))

    particles += velocity

    # Boundary clamping
    out = (particles < x_lo) | (particles > x_hi)
    particles = np.clip(particles, x_lo, x_hi)
    velocity[out] *= -0.5

    new_gains = np.sum(function(particles) * objective_weights[:, np.newaxis], axis=0)

    # Update personal bests
    idx = np.where(new_gains > pbest_gains)[0]
    pbest_positions[:, idx] = particles[:, idx]
    pbest_gains[idx]        = new_gains[idx]

    # Update global best
    if np.max(new_gains) > swarm_best_gain:
        swarm_best_position = particles[:, np.argmax(new_gains)].copy()
        swarm_best_gain     = np.max(new_gains)

    best_gain_history.append(swarm_best_gain)   # at end of each iteration
    print(f'Iteration {i+1}\tGain: {swarm_best_gain:.6f}')




# --- assumes l (N_PARAMS, POP_SIZE, MAX_ITER) is already computed ---
# --- and swarm_best_gain history is stored ---

# Store best gain per iteration (add this inside your loop):

# ── 1. Static grid: particle positions at each iteration ──────────────────
def plot_position_grid(l, MAX_ITER, x_lo=-10, x_hi=3, cols=6):
    rows = int(np.ceil(MAX_ITER / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2.5, rows * 2.5))
    axes = axes.flatten()

    for i in range(MAX_ITER):
        ax = axes[i]
        ax.scatter(l[0, :, i], l[1, :, i], s=18, color='steelblue', alpha=0.7)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(x_lo, x_hi)
        ax.set_title(f'Iter {i+1}', fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])

    for j in range(MAX_ITER, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle('Particle Positions per Iteration', fontsize=13, y=1.01)
    plt.tight_layout()
    plt.show()


# ── 2. Convergence curve ───────────────────────────────────────────────────
def plot_convergence(best_gain_history):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(1, len(best_gain_history) + 1), best_gain_history,
            color='tomato', linewidth=2, marker='o', markersize=4)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Best Gain (weighted)')
    ax.set_title('PSO Convergence')
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()


# ── 3. Animated scatter (particle movement over iterations) ───────────────
def animate_swarm(l, best_gain_history, x_lo=-10, x_hi=3, interval=200):
    fig = plt.figure(figsize=(10, 5))
    gs  = GridSpec(1, 2, figure=fig, width_ratios=[1, 1])

    ax_swarm = fig.add_subplot(gs[0])
    ax_conv  = fig.add_subplot(gs[1])

    # Swarm panel
    scat = ax_swarm.scatter(l[0, :, 0], l[1, :, 0],
                             s=40, color='steelblue', alpha=0.8, zorder=3)
    ax_swarm.set_xlim(x_lo, x_hi)
    ax_swarm.set_ylim(x_lo, x_hi)
    ax_swarm.set_xlabel('x1')
    ax_swarm.set_ylabel('x2')
    ax_swarm.grid(True, linestyle='--', alpha=0.4)
    iter_text = ax_swarm.set_title('Iteration 1')

    # Convergence panel
    conv_line, = ax_conv.plot([], [], color='tomato', linewidth=2)
    ax_conv.set_xlim(1, len(best_gain_history))
    ax_conv.set_ylim(min(best_gain_history) * 1.1,
                     max(best_gain_history) * 0.9 if max(best_gain_history) < 0
                     else max(best_gain_history) * 1.1)
    ax_conv.set_xlabel('Iteration')
    ax_conv.set_ylabel('Best Gain')
    ax_conv.set_title('Convergence')
    ax_conv.grid(True, linestyle='--', alpha=0.4)

    def update(frame):
        scat.set_offsets(l[:, :, frame].T)          # (POP_SIZE, 2)
        iter_text.set_text(f'Iteration {frame + 1}')
        conv_line.set_data(range(1, frame + 2), best_gain_history[:frame + 1])
        return scat, iter_text, conv_line

    ani = animation.FuncAnimation(fig, update, frames=l.shape[2],
                                  interval=interval, blit=True)
    plt.tight_layout()
    plt.show()
    return ani   # keep reference alive


# ── 4. Objective function contour + particle trails ───────────────────────
def plot_contour_with_trails(l, x_lo=-10, x_hi=3, resolution=200):
    x1 = np.linspace(x_lo, x_hi, resolution)
    x2 = np.linspace(x_lo, x_hi, resolution)
    X1, X2 = np.meshgrid(x1, x2)

    # Weighted gain on the grid (uses your objective weights)
    Z1 = -(X1**2 + 0.5*X2**2)
    Z2 = -(2*X1 + X2**2)
    Z  = 0.6 * Z1 + 0.4 * Z2    # matches objective_weights

    fig, ax = plt.subplots(figsize=(7, 6))
    cf = ax.contourf(X1, X2, Z, levels=40, cmap='RdYlGn')
    plt.colorbar(cf, ax=ax, label='Weighted gain')

    # Particle trails: faint lines showing movement over all iterations
    for p in range(l.shape[1]):
        ax.plot(l[0, p, :], l[1, p, :],
                color='white', alpha=0.15, linewidth=0.8)
        ax.scatter(l[0, p, -1], l[1, p, -1],
                   s=30, color='steelblue', zorder=5)

    # Mark start and end of first particle for reference
    ax.scatter(l[0, :, 0],  l[1, :, 0],  s=20, color='black',
               zorder=6, label='Start', marker='x')

    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_title('Particle Trails over Objective Landscape')
    ax.legend()
    plt.tight_layout()
    plt.show()


# Then call:
plot_position_grid(l, MAX_ITER)
plot_convergence(best_gain_history)
ani = animate_swarm(l, best_gain_history)
plot_contour_with_trails(l)