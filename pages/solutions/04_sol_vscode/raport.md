#!/usr/bin/env python3
"""Generowanie wykresow do raportu o rzucie ukosnym (indeks 3697).

Model oporu: kwadratowy opor aerodynamiczny F_d = -1/2 rho Cd S |v| v.
Wyniki zapisywane sa do katalogu 'wykresy_rzut_ukosny' obok tego pliku.

Wymagane pakiety:
    python -m pip install numpy scipy matplotlib
Uruchomienie:
    python generuj_wykresy_rzut_ukosny.py
Opcjonalnie:
    python generuj_wykresy_rzut_ukosny.py --show
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# 1. Parametry modelu i stałe
# ---------------------------------------------------------------------------

V0 = 36.0                     # predkosc poczatkowa [m/s], indeks: 36xx
ALPHA_DEG = 97.0              # kat poczatkowy [stopnie], indeks: xx97
ALPHA = np.deg2rad(ALPHA_DEG) # kat w radianach
G = 9.81                      # przyspieszenie ziemskie [m/s^2]
Y0 = 0.0                      # wysokosc poczatkowa [m]

# Przykladowe cialo: kula o parametrach zblizonych do pilki baseballowej.
# Sa to zalozenia modelowe, a nie dane pomiarowe konkretnego obiektu.
M = 0.145                     # masa [kg]
DIAMETER = 0.073              # srednica [m]
RHO = 1.225                   # gestosc powietrza [kg/m^3]
CD = 0.47                     # staly wspolczynnik oporu [-]
AREA = np.pi * DIAMETER**2 / 4.0  # pole przekroju [m^2]
BETA = 0.5 * RHO * CD * AREA / M  # parametr oporu [1/m]

VX0 = V0 * np.cos(ALPHA)
VY0 = V0 * np.sin(ALPHA)

OUT = Path(__file__).resolve().parent / "wykresy_rzut_ukosny"
OUT.mkdir(parents=True, exist_ok=True)

# Jednolity styl wykresow.
BLUE = "#2563A6"
ORANGE = "#D47722"
GREEN = "#168578"
PURPLE = "#8455A4"
RED = "#BE4B49"
DARK = "#233044"
GREY = "#64748B"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelcolor": DARK,
    "text.color": DARK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.22,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def zapisz(fig: plt.Figure, nazwa: str, show: bool = False) -> None:
    """Zapisz wykres PNG i opcjonalnie wyswietl go na ekranie."""
    sciezka = OUT / f"{nazwa}.png"
    fig.savefig(sciezka, dpi=190, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    print(f"Zapisano: {sciezka}")


# ---------------------------------------------------------------------------
# 2. Rozwiazanie analityczne bez oporu
# ---------------------------------------------------------------------------

T_IDEAL = 2.0 * VY0 / G
T_APEX_IDEAL = VY0 / G
H_MAX_IDEAL = VY0**2 / (2.0 * G)
X_APEX_IDEAL = VX0 * T_APEX_IDEAL
X_LAND_IDEAL = VX0 * T_IDEAL


def rozwiazanie_idealne(t: np.ndarray) -> tuple[np.ndarray, ...]:
    """Zwraca x, y, vx, vy, ax, ay dla modelu bez oporu."""
    x = VX0 * t
    y = Y0 + VY0 * t - 0.5 * G * t**2
    vx = np.full_like(t, VX0, dtype=float)
    vy = VY0 - G * t
    ax = np.zeros_like(t, dtype=float)
    ay = np.full_like(t, -G, dtype=float)
    return x, y, vx, vy, ax, ay


# ---------------------------------------------------------------------------
# 3. Numeryczne rozwiazanie z kwadratowym oporem powietrza
# ---------------------------------------------------------------------------


def rownania_z_oporem(_t: float, stan: np.ndarray) -> list[float]:
    """Prawa strona ukladu ODE dla stanu [x, y, vx, vy]."""
    _x, _y, vx, vy = stan
    speed = float(np.hypot(vx, vy))
    ax = -BETA * speed * vx
    ay = -G - BETA * speed * vy
    return [vx, vy, ax, ay]


def zdarzenie_apogeum(_t: float, stan: np.ndarray) -> float:
    """Apogeum: skladowa pionowa predkosci przechodzi przez zero w dol."""
    return float(stan[3])


def zdarzenie_ladowania(_t: float, stan: np.ndarray) -> float:
    """Ladowanie na poziomie startu: y=0 przy opadaniu."""
    return float(stan[1] - Y0)


zdarzenie_apogeum.terminal = False  # type: ignore[attr-defined]
zdarzenie_apogeum.direction = -1    # type: ignore[attr-defined]
zdarzenie_ladowania.terminal = True # type: ignore[attr-defined]
zdarzenie_ladowania.direction = -1  # type: ignore[attr-defined]



def oblicz_ruch_z_oporem(beta: float = BETA):
    """Calka ruchu do chwili ladowania; zwraca solver i zdarzenia."""
    def rhs(t: float, stan: np.ndarray) -> list[float]:
        _x, _y, vx, vy = stan
        speed = float(np.hypot(vx, vy))
        return [vx, vy, -beta * speed * vx, -G - beta * speed * vy]

    sol = solve_ivp(
        rhs,
        (0.0, 30.0),
        [0.0, Y0, VX0, VY0],
        events=(zdarzenie_apogeum, zdarzenie_ladowania),
        dense_output=True,
        rtol=1e-10,
        atol=1e-12,
        max_step=0.01,
    )
    if not sol.success:
        raise RuntimeError(f"Calowanie numeryczne nie powiodlo sie: {sol.message}")
    if len(sol.t_events[0]) == 0 or len(sol.t_events[1]) == 0:
        raise RuntimeError("Nie wykryto apogeum lub chwili ladowania.")
    return sol


SOL = oblicz_ruch_z_oporem()
T_APEX_DRAG = float(SOL.t_events[0][0])
T_DRAG = float(SOL.t_events[1][0])
X_APEX_DRAG, Y_APEX_DRAG, VX_APEX_DRAG, VY_APEX_DRAG = map(
    float, SOL.sol(T_APEX_DRAG)
)
X_LAND_DRAG, Y_LAND_DRAG, VX_LAND_DRAG, VY_LAND_DRAG = map(
    float, SOL.y_events[1][0]
)
V_LAND_DRAG = float(np.hypot(VX_LAND_DRAG, VY_LAND_DRAG))

# Wspolna gestsza siatka czasu dla wykresow.
t_ideal = np.linspace(0.0, T_IDEAL, 1600)
xi, yi, vxi, vyi, axi, ayi = rozwiazanie_idealne(t_ideal)
t_drag = np.linspace(0.0, T_DRAG, 1400)
xd, yd, vxd, vyd = SOL.sol(t_drag)
xd, yd, vxd, vyd = map(np.asarray, (xd, yd, vxd, vyd))
speed_i = np.hypot(vxi, vyi)
speed_d = np.hypot(vxd, vyd)
axd = -BETA * speed_d * vxd
ayd = -G - BETA * speed_d * vyd

# Wielkosci energetyczne i aerodynamiczne.
Ek_i = 0.5 * M * speed_i**2
Ep_i = M * G * yi
E_i = Ek_i + Ep_i
Ek_d = 0.5 * M * speed_d**2
Ep_d = M * G * yd
E_d = Ek_d + Ep_d
Fd_d = 0.5 * RHO * CD * AREA * speed_d**2

# ---------------------------------------------------------------------------
# 4. Rysunek 1: schemat wektorow i ukladu wspolrzednych
# ---------------------------------------------------------------------------


def rysuj_schemat(show: bool) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.0), layout="constrained")
    ax.axhline(0, color=GREY, lw=1.1)
    ax.axvline(0, color=GREY, lw=1.1)
    ax.annotate("", xy=(1.8, 0), xytext=(-1.1, 0),
                arrowprops={"arrowstyle": "->", "color": DARK, "lw": 1.4})
    ax.annotate("", xy=(0, 1.85), xytext=(0, -0.25),
                arrowprops={"arrowstyle": "->", "color": DARK, "lw": 1.4})
    ax.text(1.82, -0.10, "$x$", fontsize=13)
    ax.text(0.10, 1.82, "$y$", fontsize=13)
    ax.scatter([0], [0], color=DARK, zorder=4)

    scale = 1.35 / V0
    vx_draw, vy_draw = VX0 * scale, VY0 * scale
    ax.annotate("", xy=(vx_draw, vy_draw), xytext=(0, 0),
                arrowprops={"arrowstyle": "-|>", "color": BLUE, "lw": 2.8})
    ax.text(vx_draw - 0.38, vy_draw + 0.09, "$\\vec v_0$", color=BLUE, fontsize=13)
    ax.annotate("", xy=(vx_draw, 0), xytext=(0, 0),
                arrowprops={"arrowstyle": "->", "color": ORANGE, "lw": 2.0})
    ax.annotate("", xy=(0, vy_draw), xytext=(0, 0),
                arrowprops={"arrowstyle": "->", "color": GREEN, "lw": 2.0})
    ax.text(vx_draw / 2 - 0.04, -0.20, "$v_{0x}<0$", color=ORANGE, ha="center")
    ax.text(0.10, vy_draw * 0.52, "$v_{0y}>0$", color=GREEN)
    ax.add_patch(Arc((0, 0), 0.92, 0.92, theta1=0, theta2=ALPHA_DEG,
                     color=PURPLE, lw=1.8))
    ax.text(0.44, 0.20, "$\\alpha=97^\\circ$", color=PURPLE, fontsize=11)
    ax.text(-0.98, 1.46, "Początkowy zwrot ruchu: w lewo i do góry", fontsize=10)
    ax.set(xlim=(-1.2, 2.0), ylim=(-0.35, 2.0), aspect="equal")
    ax.set_title("Schemat rzutu ukośnego i przyjęty układ osi", loc="left")
    ax.set_xlabel("Kierunek dodatni osi x: w prawo")
    ax.set_ylabel("Kierunek dodatni osi y: do góry")
    zapisz(fig, "01_schemat_wektorow", show)


# ---------------------------------------------------------------------------
# 5. Rysunek 2: tory lotu bez oporu i z oporem
# ---------------------------------------------------------------------------


def rysuj_tory(show: bool) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.6), layout="constrained")
    ax.plot(xi, yi, color=BLUE, label="Bez oporu powietrza")
    ax.plot(xd, yd, color=ORANGE, label="Z oporem kwadratowym")
    ax.scatter([X_APEX_IDEAL], [H_MAX_IDEAL], color=BLUE, zorder=4)
    ax.scatter([X_APEX_DRAG], [Y_APEX_DRAG], color=ORANGE, zorder=4)
    ax.scatter([X_LAND_IDEAL], [Y0], color=BLUE, marker="x", s=65, zorder=4)
    ax.scatter([X_LAND_DRAG], [Y_LAND_DRAG], color=ORANGE, marker="x", s=65, zorder=4)
    ax.annotate("Apogeum bez oporu", (X_APEX_IDEAL, H_MAX_IDEAL),
                xytext=(-5, 8), textcoords="offset points", color=BLUE)
    ax.annotate("Apogeum z oporem", (X_APEX_DRAG, Y_APEX_DRAG),
                xytext=(10, -18), textcoords="offset points", color=ORANGE)
    ax.axhline(Y0, color=DARK, lw=1.0)
    ax.set(xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]")
    ax.set_title("Porównanie torów lotu dla kąta 97°", loc="left")
    ax.legend()
    ax.text(0.02, 0.96, "x < 0 oznacza ruch w lewo od punktu startu",
            transform=ax.transAxes, va="top", fontsize=9, color=GREY)
    zapisz(fig, "02_tory_lotu_porownanie", show)


# ---------------------------------------------------------------------------
# 6. Rysunek 3: polozenie w funkcji czasu
# ---------------------------------------------------------------------------


def rysuj_polozenie(show: bool) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6.4), sharex=False,
                             layout="constrained")
    axes[0].plot(t_ideal, xi, color=BLUE, label="Bez oporu")
    axes[0].plot(t_drag, xd, color=ORANGE, label="Z oporem")
    axes[0].set(ylabel="Położenie x [m]", title="Ruch poziomy")
    axes[0].legend()
    axes[1].plot(t_ideal, yi, color=BLUE, label="Bez oporu")
    axes[1].plot(t_drag, yd, color=ORANGE, label="Z oporem")
    axes[1].axhline(0, color=DARK, lw=0.8)
    axes[1].set(xlabel="Czas t [s]", ylabel="Wysokość y [m]", title="Ruch pionowy")
    axes[1].legend()
    fig.suptitle("Położenie ciała w funkcji czasu", fontsize=14, fontweight="bold")
    zapisz(fig, "03_polozenie_czas", show)


# ---------------------------------------------------------------------------
# 7. Rysunek 4: skladowe predkosci i szybkosc
# ---------------------------------------------------------------------------


def rysuj_predkosc(show: bool) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.0), sharex=True,
                             layout="constrained")
    for ax, yi_, yd_, label, title in [
        (axes[0], vxi, vxd, "$v_x$ [m/s]", "Składowa pozioma"),
        (axes[1], vyi, vyd, "$v_y$ [m/s]", "Składowa pionowa"),
        (axes[2], speed_i, speed_d, "$|v|$ [m/s]", "Szybkość"),
    ]:
        ax.plot(t_ideal, yi_, color=BLUE, label="Bez oporu")
        ax.plot(t_drag, yd_, color=ORANGE, label="Z oporem")
        ax.axhline(0, color=DARK, lw=0.7, alpha=0.6)
        ax.set_ylabel(label)
        ax.set_title(title, loc="left", fontsize=10)
    axes[0].legend(ncol=2)
    axes[-1].set_xlabel("Czas t [s]")
    fig.suptitle("Prędkość podczas lotu", fontsize=14, fontweight="bold")
    zapisz(fig, "04_predkosc_czas", show)


# ---------------------------------------------------------------------------
# 8. Rysunek 5: przyspieszenie i sila oporu
# ---------------------------------------------------------------------------


def rysuj_przyspieszenie_opor(show: bool) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.0), sharex=True,
                             layout="constrained")
    axes[0].plot(t_drag, axd, color=PURPLE, label="Z oporem")
    axes[0].axhline(0, color=BLUE, ls="--", label="Bez oporu")
    axes[0].set_ylabel("$a_x$ [m/s²]")
    axes[0].set_title("Przyspieszenie poziome", loc="left", fontsize=10)
    axes[1].plot(t_drag, ayd, color=ORANGE, label="Z oporem")
    axes[1].axhline(-G, color=BLUE, ls="--", label="Bez oporu: $-g$")
    axes[1].set_ylabel("$a_y$ [m/s²]")
    axes[1].set_title("Przyspieszenie pionowe", loc="left", fontsize=10)
    axes[2].plot(t_drag, Fd_d, color=GREEN)
    axes[2].set_ylabel("$|F_d|$ [N]")
    axes[2].set_xlabel("Czas t [s]")
    axes[2].set_title("Wartość siły oporu aerodynamicznego", loc="left", fontsize=10)
    axes[0].legend(ncol=2)
    axes[1].legend()
    fig.suptitle("Skutki działania oporu powietrza", fontsize=14, fontweight="bold")
    zapisz(fig, "05_przyspieszenie_i_sila_oporu", show)


# ---------------------------------------------------------------------------
# 9. Rysunek 6: energia kinetyczna, potencjalna i mechaniczna
# ---------------------------------------------------------------------------


def rysuj_energie(show: bool) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 5.6), layout="constrained")
    ax.plot(t_ideal, Ek_i, color=ORANGE, label="$E_k$ bez oporu")
    ax.plot(t_ideal, Ep_i, color=BLUE, label="$E_p$ bez oporu")
    ax.plot(t_ideal, E_i, color=GREEN, ls="--", label="$E_m$ bez oporu")
    ax.plot(t_drag, Ek_d, color=ORANGE, ls=":", lw=2.7, label="$E_k$ z oporem")
    ax.plot(t_drag, Ep_d, color=BLUE, ls=":", lw=2.7, label="$E_p$ z oporem")
    ax.plot(t_drag, E_d, color=RED, lw=2.7, label="$E_m$ z oporem")
    ax.set(xlabel="Czas t [s]", ylabel="Energia [J]")
    ax.set_title("Bilans energii mechanicznej podczas lotu", loc="left")
    ax.legend(ncol=2, fontsize=9)
    ax.text(0.02, 0.03,
            f"Ubytek energii mechanicznej w modelu z oporem: "
            f"{100*(E_d[0]-E_d[-1])/E_d[0]:.1f}%",
            transform=ax.transAxes, fontsize=9, color=RED)
    zapisz(fig, "06_energia", show)


# ---------------------------------------------------------------------------
# 10. Rysunek 7: wrazliwosc toru na wspolczynnik oporu Cd
# ---------------------------------------------------------------------------


def rozwiaz_dla_cd(cd: float):
    beta = 0.5 * RHO * cd * AREA / M

    def rhs(t: float, stan: np.ndarray) -> list[float]:
        _x, _y, vx, vy = stan
        sp = float(np.hypot(vx, vy))
        return [vx, vy, -beta * sp * vx, -G - beta * sp * vy]

    event = zdarzenie_ladowania
    sol = solve_ivp(
        rhs, (0.0, 30.0), [0.0, Y0, VX0, VY0],
        events=event, dense_output=True, rtol=1e-9, atol=1e-11,
        max_step=0.015,
    )
    if not sol.success or len(sol.t_events[0]) == 0:
        raise RuntimeError(f"Nie udalo sie obliczyc toru dla Cd={cd}.")
    return sol


def rysuj_wrazliwosc(show: bool) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.6), layout="constrained")
    warianty = [
        (0.0, BLUE, "Cd = 0, bez oporu"),
        (0.5 * CD, GREEN, "Cd = 0,235"),
        (CD, ORANGE, "Cd = 0,470"),
        (1.5 * CD, PURPLE, "Cd = 0,705"),
    ]
    for cd, color, label in warianty:
        if cd == 0:
            ax.plot(xi, yi, color=color, label=label)
        else:
            sol = rozwiaz_dla_cd(cd)
            tt = np.linspace(0.0, float(sol.t_events[0][0]), 800)
            q = sol.sol(tt)
            ax.plot(q[0], q[1], color=color, label=label)
    ax.axhline(0, color=DARK, lw=0.9)
    ax.set(xlabel="Położenie poziome x [m]", ylabel="Wysokość y [m]")
    ax.set_title("Wpływ przyjętego współczynnika oporu na tor lotu", loc="left")
    ax.legend()
    ax.text(0.02, 0.03,
            "Pozostałe parametry, w tym masa i średnica, są stałe.",
            transform=ax.transAxes, fontsize=9, color=GREY)
    zapisz(fig, "07_wrazliwosc_wspolczynnik_oporu", show)


# ---------------------------------------------------------------------------
# 11. Raport liczb i kontrole numeryczne
# ---------------------------------------------------------------------------


def drukuj_wyniki() -> None:
    e0 = float(E_d[0])
    e_land = float(E_d[-1])
    loss = (e0 - e_land) / e0
    print("\nPARAMETRY WEJSCIOWE")
    print(f"v0 = {V0:.3f} m/s; alpha = {ALPHA_DEG:.3f} deg = {ALPHA:.6f} rad")
    print(f"vx0 = {VX0:.6f} m/s; vy0 = {VY0:.6f} m/s")
    print(f"g = {G:.3f} m/s^2; m = {M:.3f} kg; d = {DIAMETER:.3f} m")
    print(f"rho = {RHO:.3f} kg/m^3; Cd = {CD:.3f}; A = {AREA:.6f} m^2")
    print(f"beta = {BETA:.8f} 1/m")
    print("\nWYNIKI: BEZ OPORU")
    print(f"czas lotu = {T_IDEAL:.6f} s")
    print(f"czas apogeum = {T_APEX_IDEAL:.6f} s")
    print(f"maksymalna wysokosc = {H_MAX_IDEAL:.6f} m")
    print(f"x apogeum = {X_APEX_IDEAL:.6f} m; x ladowania = {X_LAND_IDEAL:.6f} m")
    print("\nWYNIKI: Z OPOREM KWADRATOWYM")
    print(f"czas lotu = {T_DRAG:.6f} s")
    print(f"czas apogeum = {T_APEX_DRAG:.6f} s")
    print(f"maksymalna wysokosc = {Y_APEX_DRAG:.6f} m")
    print(f"x apogeum = {X_APEX_DRAG:.6f} m; x ladowania = {X_LAND_DRAG:.6f} m")
    print(f"predkosc przy ladowaniu = {V_LAND_DRAG:.6f} m/s")
    print(f"energia mechaniczna na starcie = {e0:.6f} J")
    print(f"energia mechaniczna przy ladowaniu = {e_land:.6f} J")
    print(f"ubytek energii mechanicznej = {100*loss:.3f} %")


def kontrole() -> None:
    # Kontrola wzorow analitycznych oraz znaku skladowej poziomej.
    assert VX0 < 0.0 and VY0 > 0.0, "Dla 97 stopni vx0<0 i vy0>0."
    assert abs(yi[-1] - Y0) < 1e-9, "Model idealny powinien konczyc sie na poziomie startu."
    assert abs(E_i[0] - E_i[-1]) < 1e-9, "Energia idealna powinna byc stala."
    assert E_d[-1] < E_d[0], "Opory powinny zmniejszac energie mechaniczna."
    assert np.max(np.diff(E_d)) < 1e-7, "Energia mechaniczna z oporem nie powinna rosnac."
    assert abs(Y_LAND_DRAG - Y0) < 1e-8, "Zdarzenie ladowania powinno spelnic y=y0."
    print("\nKontrole numeryczne: OK")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true", help="pokaz wykresy podczas generowania")
    args = parser.parse_args()

    rysuj_schemat(args.show)
    rysuj_tory(args.show)
    rysuj_polozenie(args.show)
    rysuj_predkosc(args.show)
    rysuj_przyspieszenie_opor(args.show)
    rysuj_energie(args.show)
    rysuj_wrazliwosc(args.show)
    kontrole()
    drukuj_wyniki()
    print(f"\nGotowe. Pliki PNG znajduja sie w: {OUT}")


if __name__ == "__main__":
    main()
