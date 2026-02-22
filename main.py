import math
import matplotlib.pyplot as plt

G = 6.67430e-11
M = 5.972e24
R = 6.371e6

dt = 0.5
t = 0.0

x = R
y = 0.0
vx = 0.0
vy = 0.0

dry_mass = 12000.0
fuel_mass = 22000.0
mass = dry_mass + fuel_mass

thrust = 900_000.0
fuel_burn_rate = 3.0

target_altitude = 200_000.0

phase = "ASCENT"

times = []
alts_km = []
speeds = []
xs_km = []
ys_km = []

prev_r = None
prev_vr = None
apo_r = None
apo_vx = None
apo_vy = None

def accel(x, y, vx, vy, mass, phase, fuel_mass, apo_r):
    r = math.hypot(x, y)
    ax_g = -G * M * x / (r**3)
    ay_g = -G * M * y / (r**3)

    ax_t = 0.0
    ay_t = 0.0
    thrusting = False

    if fuel_mass > 0 and phase in ("ASCENT", "CIRCULARIZE"):
        rx = x / r
        ry = y / r
        tx = -y / r
        ty = x / r

        if phase == "ASCENT":
            alt = r - R
            pitch = max(0.0, min(1.0, alt / 150_000.0))
            ux = (1 - pitch) * rx + pitch * tx
            uy = (1 - pitch) * ry + pitch * ty
            um = math.hypot(ux, uy)
            ux /= um
            uy /= um
            ax_t = thrust * ux / mass
            ay_t = thrust * uy / mass
            thrusting = True

        elif phase == "CIRCULARIZE":
            v_circ = math.sqrt(G * M / apo_r)
            v = math.hypot(vx, vy)
            if v < 0.9997 * v_circ:
                ux = tx
                uy = ty
                ax_t = thrust * ux / mass
                ay_t = thrust * uy / mass
                thrusting = True

    return ax_g + ax_t, ay_g + ay_t, thrusting

def rk4_step(x, y, vx, vy, mass, phase, fuel_mass, apo_r):
    ax1, ay1, th1 = accel(x, y, vx, vy, mass, phase, fuel_mass, apo_r)
    k1x, k1y, k1vx, k1vy = vx, vy, ax1, ay1

    ax2, ay2, th2 = accel(
        x + 0.5 * dt * k1x,
        y + 0.5 * dt * k1y,
        vx + 0.5 * dt * k1vx,
        vy + 0.5 * dt * k1vy,
        mass, phase, fuel_mass, apo_r
    )
    k2x, k2y, k2vx, k2vy = vx + 0.5 * dt * k1vx, vy + 0.5 * dt * k1vy, ax2, ay2

    ax3, ay3, th3 = accel(
        x + 0.5 * dt * k2x,
        y + 0.5 * dt * k2y,
        vx + 0.5 * dt * k2vx,
        vy + 0.5 * dt * k2vy,
        mass, phase, fuel_mass, apo_r
    )
    k3x, k3y, k3vx, k3vy = vx + 0.5 * dt * k2vx, vy + 0.5 * dt * k2vy, ax3, ay3

    ax4, ay4, th4 = accel(
        x + dt * k3x,
        y + dt * k3y,
        vx + dt * k3vx,
        vy + dt * k3vy,
        mass, phase, fuel_mass, apo_r
    )
    k4x, k4y, k4vx, k4vy = vx + dt * k3vx, vy + dt * k3vy, ax4, ay4

    nx = x + (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
    ny = y + (dt / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y)
    nvx = vx + (dt / 6.0) * (k1vx + 2 * k2vx + 2 * k3vx + k4vx)
    nvy = vy + (dt / 6.0) * (k1vy + 2 * k2vy + 2 * k3vy + k4vy)

    thrusting = th1 or th2 or th3 or th4
    return nx, ny, nvx, nvy, thrusting

for step in range(200000):
    r = math.hypot(x, y)
    alt = r - R
    v = math.hypot(vx, vy)

    if r <= R:
        break

    vr = (x * vx + y * vy) / r

    if phase == "ASCENT" and alt >= target_altitude:
        phase = "COAST"
        prev_r = r
        prev_vr = vr

    if phase == "COAST":
        if prev_vr is not None and prev_vr > 0 and vr <= 0:
            apo_r = r
            apo_vx = vx
            apo_vy = vy
            phase = "CIRCULARIZE"

        prev_r = r
        prev_vr = vr

    x, y, vx, vy, thrusting = rk4_step(x, y, vx, vy, mass, phase, fuel_mass, apo_r if apo_r else r)

    if thrusting and fuel_mass > 0:
        fuel_mass -= fuel_burn_rate * dt
        if fuel_mass < 0:
            fuel_mass = 0.0
        mass = dry_mass + fuel_mass

    if phase == "CIRCULARIZE" and apo_r is not None:
        v_circ = math.sqrt(G * M / apo_r)
        if math.hypot(vx, vy) >= 0.9997 * v_circ:
            phase = "ORBIT"

    times.append(t)
    alts_km.append(alt / 1000.0)
    speeds.append(v)
    xs_km.append(x / 1000.0)
    ys_km.append(y / 1000.0)

    t += dt

plt.figure()
plt.plot(times, alts_km)
plt.xlabel("Time (s)")
plt.ylabel("Altitude (km)")
plt.title("Altitude vs Time")
plt.savefig("altitude_vs_time.png")

plt.figure()
plt.plot(times, speeds)
plt.xlabel("Time (s)")
plt.ylabel("Speed (m/s)")
plt.title("Speed vs Time")
plt.savefig("speed_vs_time.png")

theta = [i * 2 * math.pi / 600 for i in range(601)]
earth_x = [R * math.cos(a) / 1000.0 for a in theta]
earth_y = [R * math.sin(a) / 1000.0 for a in theta]

plt.figure()
plt.plot(earth_x, earth_y)
plt.plot(xs_km, ys_km)
plt.gca().set_aspect("equal", adjustable="box")
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Orbit Path (2D)")
plt.savefig("orbit_path.png")