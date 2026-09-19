#!/usr/bin/env python3
"""
Figura: Ebf2 w komórkach Purkinjego — gradient płatowy bez pasów parasagittalnych.

Wejście:  processed/abc_ebf2_komorki.csv, abc_ebf2_per_lobule.csv, lobule_composition.csv
Wyjście:  figures/abc_ebf2.png

Panele:
  A, B  dwie sąsiednie sekcje (~200 µm) — magnituda Ebf2, rampa sekwencyjna
  C     zgodność sąsiadów w obu testach, z kontrolą pozytywną i negatywną
  D     Ebf2 vs udział Aldoc+ po płacikach

Kolor przypisany po zadaniu (dataviz):
  magnituda  -> jedna rampa błękitna, jasny = bliski zeru
  identyczność w panelu C -> emfaza: sygnał #2a78d6, kontrola pozytywna #eb6834,
                             kontrole negatywne neutralny szary
  panel D ma jedną serię -> bez legendy, tytuł ją nazywa
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

P = Path(__file__).resolve().parent.parent / "processed"
FIG = Path(__file__).resolve().parent.parent / "figures"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
INK3 = "#8a8880"
SYGNAL = "#2a78d6"      # slot 1
KONTROLA_POZ = "#eb6834"  # slot 2
NEUTRAL = "#b8b6ae"
GRID = "#e6e5df"
# rampa sekwencyjna: blue 100 -> 700
RAMPA = LinearSegmentedColormap.from_list(
    "blue_seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#123a6b"])

MAPA = {"LING": "I", "CENT2": "II", "CENT3": "III", "CUL4, 5": "CUL", "DEC": "VI",
        "FOTU": "VII", "PYR": "VIII", "UVU": "IX", "NOD": "X", "SIM": "SIM",
        "ANcr1": "AN1", "ANcr2": "AN2", "PRM": "PRM", "COPY": "COP", "PFL": "PF", "FL": "F"}


def styl(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8, length=3, width=0.8)


def panel_mapa(ax, d, sekcja, cmap_norm, tytul):
    g = d[d["brain_section_label"] == sekcja]
    ax.scatter(g["z_ccf"], g["y_ccf"], c=g["Ebf2"], cmap=RAMPA,
               vmin=cmap_norm[0], vmax=cmap_norm[1], s=5, linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(tytul, fontsize=10.5, color=INK, loc="left", pad=8, fontweight="bold")
    ax.set_xlabel("oś przyśrodkowo-boczna (mm)", fontsize=8, color=INK2)
    ax.set_ylabel("oś grzbietowo-brzuszna (mm)", fontsize=8, color=INK2)
    styl(ax)


def main():
    d = pd.read_csv(P / "abc_ebf2_komorki.csv")
    lob = pd.read_csv(P / "abc_ebf2_per_lobule.csv")
    koz = pd.read_csv(P / "lobule_composition.csv")
    koz["lobule"] = koz["atlas"].map(MAPA)

    fig = plt.figure(figsize=(11.6, 8.8), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.24,
                          left=0.085, right=0.955, top=0.875, bottom=0.155)

    fig.suptitle("Ebf2 w komórkach Purkinjego myszy: gradient płatowy, nie pasy zebrinowe",
                 fontsize=14, color=INK, x=0.075, ha="left", y=0.965, fontweight="bold")
    fig.text(0.075, 0.925,
             "ABC Atlas MERFISH (Yao 2023), 16 626 komórek Purkinjego w 14 sekcjach koronalnych, "
             "współrzędne CCFv3 · skrypty 49–50",
             fontsize=9, color=INK2, ha="left")

    # --- A, B: dwie sąsiednie sekcje ---
    lim = (0, float(np.percentile(d["Ebf2"], 98)))
    par = ["C57BL6J-638850.14", "C57BL6J-638850.15"]
    n0 = int((d["brain_section_label"] == par[0]).sum())
    n1 = int((d["brain_section_label"] == par[1]).sum())
    opisy = [f"A · sekcja 14   ({n0} komórek)",
             f"B · sekcja 15 — sąsiednia   ({n1} komórek)"]
    for ax_i, (sek, t) in enumerate(zip(par, opisy)):
        ax = fig.add_subplot(gs[0, ax_i])
        panel_mapa(ax, d, sek, lim, t)
        if ax_i == 1:
            sm = plt.cm.ScalarMappable(cmap=RAMPA,
                                       norm=plt.Normalize(vmin=lim[0], vmax=lim[1]))
            cb = fig.colorbar(sm, ax=ax, fraction=0.036, pad=0.03)
            cb.set_label("Ebf2 (log₂, na jednostkę objętości)", fontsize=8, color=INK2)
            cb.ax.tick_params(labelsize=7.5, colors=INK2, length=2.5, width=0.7)
            cb.outline.set_visible(False)

    # --- C: zgodność sąsiadów ---
    ax = fig.add_subplot(gs[1, 0])
    dane = [
        ("pasy 500 µm", 91.2, KONTROLA_POZ),
        ("pasy 300 µm", 82.9, KONTROLA_POZ),
        ("Ebf2", 52.6, SYGNAL),
        ("Ebf1", 51.1, NEUTRAL),
        ("Ebf3", 49.0, NEUTRAL),
        (None, None, None),                # przerwa między testami
        ("Ebf2", 57.7, SYGNAL),
        ("Ebf1", 52.4, NEUTRAL),
        ("Ebf3", 50.5, NEUTRAL),
    ]
    y = np.arange(len(dane))[::-1]
    etyk_y, etyk_t = [], []
    for yi, (nazwa, v, kolor) in zip(y, dane):
        if nazwa is None:
            continue
        ax.barh(yi, v - 50, left=50, height=0.62, color=kolor,
                edgecolor=SURFACE, linewidth=1.4)
        ax.text(v + 0.8, yi, f"{v:.1f} %", va="center", fontsize=8.5,
                color=INK, fontweight="bold" if kolor != NEUTRAL else "normal")
        etyk_y.append(yi)
        etyk_t.append(nazwa)
    ax.axvline(50, color=INK3, linewidth=1.2, linestyle=(0, (4, 3)), zorder=3)
    ax.text(50.3, -0.62, "permutacja = brak struktury", fontsize=7.5,
            color=INK3, va="center", ha="left")
    # nagłówki grup
    ax.text(46.6, y[0] + 0.62, "TEST B — ciągłość między sekcjami", fontsize=8,
            color=INK, fontweight="bold", va="center", ha="left")
    ax.text(46.6, y[6] + 0.62, "TEST A — skupianie w płaszczyźnie", fontsize=8,
            color=INK, fontweight="bold", va="center", ha="left")
    ax.set_yticks(etyk_y)
    ax.set_yticklabels(etyk_t, fontsize=8.5, color=INK2)
    ax.set_xlim(46, 97)
    ax.set_ylim(-0.7, len(dane) + 0.15)
    ax.set_xlabel("zgodność etykiet między sąsiadującymi komórkami (%)",
                  fontsize=8, color=INK2)
    ax.set_title("C · test wykrywa narzucone pasy 300 µm — Ebf2 ich nie ma",
                 fontsize=10.5, color=INK, loc="left", pad=10, fontweight="bold")
    ax.xaxis.grid(True, color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    styl(ax)
    ax.legend(handles=[Line2D([], [], marker="s", linestyle="", markersize=7.5,
                              color=c, label=l)
                       for c, l in [(SYGNAL, "Ebf2"), (KONTROLA_POZ, "kontrola pozytywna"),
                                    (NEUTRAL, "Ebf1 / Ebf3 (kontrola negatywna)")]],
              fontsize=7.8, frameon=False, loc="lower right", labelcolor=INK2,
              handletextpad=0.5, borderaxespad=0.4)

    # --- D: Ebf2 vs Aldoc+ po płacikach ---
    ax = fig.add_subplot(gs[1, 1])
    t = lob.merge(koz[["lobule", "Aldoc+_frac"]], on="lobule")
    x, yv = t["Aldoc+_frac"] * 100, t["Ebf2_pct"]
    r, p = pearsonr(x, yv)
    b, a = np.polyfit(x, yv, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, a + b * xs, color=INK3, linewidth=1.6, zorder=2)
    ax.scatter(x, yv, s=np.clip(t["n"] / 9, 22, 190), color=SYGNAL, alpha=0.85,
               linewidths=1.6, edgecolors=SURFACE, zorder=3)
    # ręczne offsety tam, gdzie punkty leżą blisko siebie
    OFF = {"PRM": (-16, -2), "AN2": (2, 11), "COP": (-3, -14), "VII": (-13, 4),
           "IX": (10, 3), "X": (-2, -15), "F": (7, 4), "PF": (10, -8),
           "VIII": (0, 11), "SIM": (-1, 11), "CUL": (0, 11), "I": (0, -14)}
    for _, row in t.iterrows():
        dx, dy = OFF.get(row["lobule"], (0, 11))
        ax.annotate(row["lobule"], (row["Aldoc+_frac"] * 100, row["Ebf2_pct"]),
                    textcoords="offset points", xytext=(dx, dy),
                    ha="center" if dx == 0 else ("right" if dx < 0 else "left"),
                    fontsize=7.8, color=INK2, zorder=4)
    ax.set_xlabel("udział komórek Aldoc+ w płaciku (%, Kozareva 2021)", fontsize=8, color=INK2)
    ax.set_ylabel("komórki Ebf2-dodatnie (%)", fontsize=8, color=INK2)
    ax.set_title("D · im więcej zebriny w płaciku, tym mniej Ebf2",
                 fontsize=10.5, color=INK, loc="left", pad=10, fontweight="bold")
    ax.text(0.975, 0.955, f"Pearson r = {r:+.3f}   p = {p:.3f}   n = {len(t)} płacików",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=INK,
            fontweight="bold")
    ax.text(0.975, 0.895, "wielkość punktu = liczba komórek",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK3)
    ax.set_ylim(15, 95)
    ax.grid(True, color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    styl(ax)

    fig.text(0.085, 0.018,
             "Panele A i B: dwie kolejne sekcje serii, odległe o ~200 µm — wzór Ebf2 nie powtarza się między nimi.\n"
             "Etykiety hi/lo wyznaczane osobno w każdej parze (sekcja × płacik), więc gradient płatowy z panelu D\n"
             "jest usunięty z testów w panelu C. Test A: 6 najbliższych sąsiadów w płaszczyźnie sekcji. Test B: "
             "najbliższa komórka w sekcji sąsiedniej. Istotność z 500 permutacji etykiet.",
             fontsize=7.6, color=INK3, ha="left", va="bottom", linespacing=1.6)

    FIG.mkdir(exist_ok=True)
    out = FIG / "abc_ebf2.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"zapisano: {out}")
    print(f"panel D: r = {r:+.4f}, p = {p:.4f}, n = {len(t)}")


if __name__ == "__main__":
    main()
