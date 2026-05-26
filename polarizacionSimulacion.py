from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import scienceplots  

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "imagenes"
OUTPUT_DIR.mkdir(exist_ok=True)

plt.style.use(["science", str(BASE_DIR / "aip.mplstyle")])

AIP_WIDTH = 3.37
ANGLE_UNC_DEG = 1.0
POWER_UNC_MW = 0.01

POLARIZER_HORIZONTAL_MARK_DEG = 52.0
POLARIZER_VERTICAL_MARK_DEG = 142.0
QWP_CIRCULAR_MARK_DEG = 54.5
HWP_MINIMUM_MARK_DEG = 190.0
SECOND_QWP_STOKES_MARK_DEG = 300.0

P_MALUS_MAX_MW = 2.105
P_MALUS_MIN_MW = 0.102
RETARDER_RELATIVE_FAST_AXIS_DEG = 45.0


def polarization(elected, campoMatriz):
    polarizado = np.zeros((campoMatriz.shape[0], campoMatriz.shape[1], 2), dtype=complex)
    if elected == "H":
        jones = np.array([1, 0])
    elif elected == "V":
        jones = np.array([0, 1])
    elif elected == "R":
        jones = (1/np.sqrt(2)) * np.array([1, 1j])
    elif elected == "L":
        jones = (1/np.sqrt(2)) * np.array([1, -1j])
    elif elected == "D":
        jones = (1/np.sqrt(2)) * np.array([1, 1])
    elif elected == "LD":
        jones = (1/np.sqrt(2)) * np.array([1, -1])
    else:
        raise ValueError("invalido")
    polarizado[:, :, 0] = campoMatriz * jones[0]
    polarizado[:, :, 1] = campoMatriz * jones[1]
    return polarizado


r = np.linspace(0, 2, 256)
phi = np.linspace(0, 2*np.pi, 256)
R, Phi = np.meshgrid(r, phi)
w0 = 1
m = 0

def campo(r, phi):
    return (r/w0)**m * np.exp(-r**2/w0**2) * np.exp(1j*m*phi)

campoMatriz = campo(R, Phi)
X = R*np.cos(Phi)
Y = R*np.sin(Phi)

E = polarization("V", campoMatriz)
Ex = E[:, :, 0]
Ey = E[:, :, 1]

theta = np.linspace(0, 2*np.pi, 360)
I_malus = np.zeros(theta.shape)
for k in range(len(theta)):
    th = theta[k] + np.pi/2
    M = np.array([[np.cos(th)**2, np.cos(th)*np.sin(th)],
                  [np.sin(th)*np.cos(th), np.sin(th)**2]])
    Ex_out = M[0, 0]*Ex + M[0, 1]*Ey
    Ey_out = M[1, 0]*Ex + M[1, 1]*Ey
    I_malus[k] = np.sum(np.abs(Ex_out)**2 + np.abs(Ey_out)**2)

I_malus_norm = I_malus / np.max(I_malus)


angulos_exp = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
                        150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270,
                        280, 290, 300, 310, 320, 330, 340, 350, 360])
intensidad_exp = np.array([2.003, 2.006, 1.807, 1.604, 1.209, 0.905, 0.608, 0.304, 0.106,
                           0.103, 0.107, 0.304, 0.508, 0.803, 1.106, 1.508, 1.803, 2.009,
                           2.105, 2.002, 1.804, 1.609, 1.305, 0.908, 0.602, 0.306, 0.105,
                           0.102, 0.109, 0.307, 0.506, 0.902, 1.307, 1.605, 1.809, 2.004,
                           2.007])
P0_exp = np.max(intensidad_exp)
I_exp_norm = intensidad_exp / P0_exp
theta_exp_rad = np.deg2rad(angulos_exp)
theta_unc_rad = np.deg2rad(ANGLE_UNC_DEG)
unc_malus = np.sqrt((np.cos(theta_exp_rad)**2 * POWER_UNC_MW)**2
                    + (P0_exp*np.sin(2*theta_exp_rad) * theta_unc_rad)**2)
unc_malus_norm = unc_malus / P0_exp

fig, ax = plt.subplots(figsize=(AIP_WIDTH, 2.25))
ax.plot(np.degrees(theta), I_malus_norm, "-",
        linewidth=1.8, label=r"Sim.")
ax.errorbar(angulos_exp, I_exp_norm,
            yerr=unc_malus_norm,
            fmt="o", markersize=1.7, markerfacecolor="white",
            markeredgewidth=0.55, capsize=2, elinewidth=0.75,
            linewidth=0.8,
            label=r"Exp.")
ax.set_xlabel(r"\'Angulo del polarizador $\theta$ (deg)")
ax.set_ylabel(r"Intensidad normalizada $I/I_0$")
ax.set_xlim([0, 360])
ax.set_ylim([-0.05, 1.1])
ax.set_xticks(np.arange(0, 361, 45))
ax.grid(True, alpha=0.3)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
          borderaxespad=0.0, columnspacing=1.2, handlelength=1.8)
fig.tight_layout(pad=0.25)
fig.savefig(OUTPUT_DIR / "sim_leyMalus.pdf", bbox_inches="tight")
plt.close(fig)

I_malus_teo = np.cos(np.deg2rad(angulos_exp))**2
residuos_malus = I_exp_norm - I_malus_teo
rmse_malus_acum = np.sqrt(np.cumsum(residuos_malus**2) / np.arange(1, residuos_malus.size + 1))
rmse_malus = rmse_malus_acum[-1]

fig, ax = plt.subplots(figsize=(AIP_WIDTH, 2.25))
ax.plot(angulos_exp, rmse_malus_acum, "o-",
        markersize=3, linewidth=1.2, label=r"RMSE acum.")
ax.axhline(rmse_malus, color="black", linestyle="--",
           linewidth=0.9, label=rf"RMSE total = {rmse_malus:.3f}")
ax.set_xlabel(r"\'Angulo del polarizador $\theta$ (deg)")
ax.set_ylabel(r"RMSE de $I/I_0$")
ax.set_xlim([0, 360])
ax.set_ylim([0, max(rmse_malus_acum)*1.15])
ax.set_xticks(np.arange(0, 361, 45))
ax.grid(True, alpha=0.3)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
          borderaxespad=0.0, columnspacing=1.2, handlelength=1.8)
fig.tight_layout(pad=0.25)
fig.savefig(OUTPUT_DIR / "sim_rmse_malus.pdf", bbox_inches="tight")
plt.close(fig)


def stokes_de_jones(J):
    a = J[0]
    b = J[1]
    return np.array([
        np.abs(a)**2 + np.abs(b)**2,
        np.abs(a)**2 - np.abs(b)**2,
        2*np.real(a*np.conj(b)),
        2*np.imag(a*np.conj(b))
    ])


J_CD = (1/np.sqrt(2)) * np.array([1, -1j])
S_CD = stokes_de_jones(J_CD)


def stokes_de_proyecciones(I0, IH, IV, ID, IA, IR, IL):
    return np.array([I0, IH - IV, ID - IA, IR - IL])


def incertidumbre_stokes_normalizada(S):
    return np.full(4, np.sqrt(2)*POWER_UNC_MW / S[0])


def grado_polarizacion(S):
    return np.sqrt(S[1]**2 + S[2]**2 + S[3]**2) / S[0]


def HWP(theta_h):
    c = np.cos(2*theta_h)
    s = np.sin(2*theta_h)
    return np.array([[c, s], [s, -c]])

J_CI_pred = HWP(np.deg2rad(RETARDER_RELATIVE_FAST_AXIS_DEG)) @ J_CD
S_CI_pred = stokes_de_jones(J_CI_pred)

S_CD_exp = stokes_de_proyecciones(
    I0=1.99, IH=0.99, IV=1.00, ID=1.01, IA=0.98, IR=1.95, IL=0.04
)
S_CD_exp_norm = S_CD_exp / S_CD_exp[0]
unc_S_CD = incertidumbre_stokes_normalizada(S_CD_exp)

S_CI_exp = stokes_de_proyecciones(
    I0=1.99, IH=0.99, IV=1.00, ID=1.00, IA=0.99, IR=0.05, IL=1.92
)
S_CI_exp_norm = S_CI_exp / S_CI_exp[0]
unc_S_CI = incertidumbre_stokes_normalizada(S_CI_exp)

u = np.linspace(0, 2*np.pi, 48)
v = np.linspace(0, np.pi, 24)
xs = np.outer(np.cos(u), np.sin(v))
ys = np.outer(np.sin(u), np.sin(v))
zs = np.outer(np.ones_like(u), np.cos(v))

def stokes_poincare(S):
    return S[1:] / S[0]

poincare_states = [
    (r"CD sim.", stokes_poincare(S_CD), "o", "C0", 24),
    (r"CD exp.", stokes_poincare(S_CD_exp_norm), "^", "C1", 54),
    (r"CI sim.", stokes_poincare(S_CI_pred), "o", "C2", 24),
    (r"CI exp.", stokes_poincare(S_CI_exp_norm), "^", "C3", 54),
]

fig = plt.figure(figsize=(AIP_WIDTH, 3.05))
ax = fig.add_subplot(111, projection="3d")
ax.plot_wireframe(xs, ys, zs, color="0.72", linewidth=0.25, alpha=0.55)
ax.plot(np.cos(u), np.sin(u), 0*u, color="0.35", linewidth=0.7)
ax.plot(np.cos(u), 0*u, np.sin(u), color="0.55", linewidth=0.5)
ax.plot(0*u, np.cos(u), np.sin(u), color="0.55", linewidth=0.5)
ax.quiver(0, 0, 0, 1.05, 0, 0, color="black", linewidth=0.6,
          arrow_length_ratio=0.06)
ax.quiver(0, 0, 0, 0, 1.05, 0, color="black", linewidth=0.6,
          arrow_length_ratio=0.06)
ax.quiver(0, 0, 0, 0, 0, 1.05, color="black", linewidth=0.6,
          arrow_length_ratio=0.06)

for label, point, marker, color, size in poincare_states:
    ax.scatter(point[0], point[1], point[2], marker=marker, s=size,
               color=color, edgecolor="black", linewidth=0.45,
               depthshade=False, label=label)
    if "exp." in label:
        ax.plot([0, point[0]], [0, point[1]], [0, point[2]],
                color=color, linestyle=":", linewidth=0.75)
        ax.text(point[0] + 0.08, point[1] + 0.05, point[2],
                label, color=color, fontsize=6)

ax.set_xlabel(r"$S_1/S_0$", labelpad=-2)
ax.set_ylabel(r"$S_2/S_0$", labelpad=-2)
ax.set_zlabel("")
ax.text2D(0.05, 0.43, r"$S_3/S_0$", transform=ax.transAxes,
          rotation=90, va="center", ha="center")
ax.set_xlim([-1.12, 1.12])
ax.set_ylim([-1.12, 1.12])
ax.set_zlim([-1.12, 1.12])
ax.set_xticks([-1, 0, 1])
ax.set_yticks([-1, 0, 1])
ax.set_zticks([-1, 0, 1])
ax.set_box_aspect((1, 1, 1))
ax.tick_params(pad=0)
ax.view_init(elev=24, azim=48)
ax.grid(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
          borderaxespad=0.0, columnspacing=0.9, handlelength=1.2)
fig.tight_layout(pad=0.35)
fig.savefig(OUTPUT_DIR / "sim_poincare.pdf", bbox_inches="tight")
plt.close(fig)

labels = [r"$S_0$", r"$S_1$", r"$S_2$", r"$S_3$"]
x = np.arange(len(labels))
width = 0.28

fig, axs = plt.subplots(1, 2, figsize=(AIP_WIDTH, 2.35), sharey=True)
axs[0].bar(x - width/2, S_CD, width,
           edgecolor="black", linewidth=0.6, label=r"Simulaci\'on")
axs[0].bar(x + width/2, S_CD_exp_norm, width,
           edgecolor="black", linewidth=0.6, label=r"Experimental",
           yerr=unc_S_CD, capsize=3)
axs[0].set_xticks(x)
axs[0].set_xticklabels(labels)
axs[0].set_ylabel(r"Valor normalizado")
axs[0].set_title(r"CD")
axs[0].axhline(0, color="black", linewidth=0.6)
axs[0].grid(True, alpha=0.3, axis="y")
axs[0].set_ylim([-1.2, 1.3])

axs[1].bar(x - width/2, S_CI_pred, width,
           edgecolor="black", linewidth=0.6, label=r"Simulaci\'on")
axs[1].bar(x + width/2, S_CI_exp_norm, width,
           edgecolor="black", linewidth=0.6, label=r"Experimental",
           yerr=unc_S_CI, capsize=3)
axs[1].set_xticks(x)
axs[1].set_xticklabels(labels)
axs[1].set_title(r"Retardador $\lambda/2$ (CI)")
axs[1].axhline(0, color="black", linewidth=0.6)
axs[1].grid(True, alpha=0.3, axis="y")
handles, legend_labels = axs[0].get_legend_handles_labels()
fig.legend(handles, legend_labels, loc="upper center", bbox_to_anchor=(0.5, 1.02),
           ncol=2, borderaxespad=0.0, columnspacing=1.2, handlelength=1.8)

fig.tight_layout(pad=0.25, rect=[0, 0, 1, 0.92])
fig.savefig(OUTPUT_DIR / "sim_stokes.pdf", bbox_inches="tight")
plt.close(fig)


theta_a = np.linspace(0, 2*np.pi, 360)
I_QWP = np.zeros(theta_a.shape)
qwp_angle = np.deg2rad(RETARDER_RELATIVE_FAST_AXIS_DEG)
M_QWP = np.array([
    [np.cos(qwp_angle)**2 - 1j*np.sin(qwp_angle)**2,
     (1 + 1j)*np.cos(qwp_angle)*np.sin(qwp_angle)],
    [(1 + 1j)*np.cos(qwp_angle)*np.sin(qwp_angle),
     np.sin(qwp_angle)**2 - 1j*np.cos(qwp_angle)**2]
])
Ex_qwp = M_QWP[0, 0]*Ex + M_QWP[0, 1]*Ey
Ey_qwp = M_QWP[1, 0]*Ex + M_QWP[1, 1]*Ey

for k in range(len(theta_a)):
    th = theta_a[k]
    M2 = np.array([[np.cos(th)**2, np.cos(-th)*np.sin(-th)],
                   [np.sin(-th)*np.cos(th), np.sin(-th)**2]])
    Ex_o = M2[0, 0]*Ex_qwp + M2[0, 1]*Ey_qwp
    Ey_o = M2[1, 0]*Ex_qwp + M2[1, 1]*Ey_qwp
    I_QWP[k] = np.sum(np.abs(Ex_o)**2 + np.abs(Ey_o)**2)

I_QWP_norm = I_QWP / np.max(I_QWP)

fig, ax = plt.subplots(figsize=(AIP_WIDTH, 2.25))
ax.plot(np.degrees(theta_a), I_malus_norm, "-",
        linewidth=1.8,
        label=r"Sin retardador (V)")
ax.plot(np.degrees(theta_a), I_QWP_norm, "--",
        linewidth=1.8,
        label=r"QWP a $45^{\circ}$")
ax.set_xlabel(r"\'Angulo del polarizador $\theta$ (deg)")
ax.set_ylabel(r"Intensidad normalizada $I/I_0$")
ax.set_xlim([0, 360])
ax.set_ylim([-0.05, 1.1])
ax.set_xticks(np.arange(0, 361, 45))
ax.grid(True, alpha=0.3)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2,
          borderaxespad=0.0, columnspacing=1.2, handlelength=1.8)
fig.tight_layout(pad=0.25)
fig.savefig(OUTPUT_DIR / "sim_QWP_circular.pdf", bbox_inches="tight")
plt.close(fig)


analyzers = np.array([
    [1, 0],
    [0, 1],
    [1, 1],
    [1, -1],
    [1, 1j],
    [1, -1j]
], dtype=complex)

I_an = np.zeros(analyzers.shape[0])
Ex_in = Ex_qwp
Ey_in = Ey_qwp
I0 = np.sum(np.abs(Ex_in)**2 + np.abs(Ey_in)**2)

for k in range(analyzers.shape[0]):
    a = analyzers[k, 0]
    b = analyzers[k, 1]
    norm = np.sqrt(np.abs(a)**2 + np.abs(b)**2)
    a /= norm
    b /= norm
    Eproj = np.conj(a)*Ex_in + np.conj(b)*Ey_in
    I_an[k] = np.sum(np.abs(Eproj)**2) / I0

IH, IV, ID, IA, IR, IL = I_an
S0_sim = IH + IV
S1_sim = IH - IV
S2_sim = ID - IA
S3_sim = IR - IL
DoP_sim = np.sqrt(S1_sim**2 + S2_sim**2 + S3_sim**2) / S0_sim
DoP_CD_exp = grado_polarizacion(S_CD_exp)
DoP_CI_exp = grado_polarizacion(S_CI_exp)

print(f"S0={S0_sim:.4f}  S1={S1_sim:.4f}  S2={S2_sim:.4f}  S3={S3_sim:.4f}")
print(f"DoP={DoP_sim:.4f}")
print(f"DoP CD exp={DoP_CD_exp:.4f}  DoP CI exp={DoP_CI_exp:.4f}")
