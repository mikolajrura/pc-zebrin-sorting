#!/usr/bin/env python3
"""
Film: obracająca się warstwa Purkinjego + panel składu fenotypowego płacika.

Wejście:  /mnt/data1t/pc_rebuild/scene_pc_layer.npz   geometria (2 279 886 wokseli)
          /mnt/data1t/pc_rebuild/web_cloud.json       metadane 16 płacików
Wyjście:  figures/film_mozdzek.mp4  (+ klatki w katalogu roboczym)

Rasteryzer własny (numpy): rotacja → projekcja perspektywiczna → sortowanie po
głębi → bufor pikseli. Panel rysowany w PIL. Bez zależności od GPU i przeglądarki.

Paleta: pierwsze 9 kolorów 16-kolorowej listy z artefaktu v1, w kolejności
`subtypes` — zweryfikowane wobec zrzutu panelu (Anti_Aldoc_2 #6fb7c4,
Anti_Aldoc_1 #d96a63, Aldoc_5 #c97ba6).

    python 58_film_mozdzek.py --test          3 klatki kontrolne
    python 58_film_mozdzek.py --frames 1200   pełny render
"""
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont

BASE = Path(__file__).resolve().parent.parent
SCENE = Path("/mnt/data1t/pc_rebuild/scene_pc_layer.npz")
META = Path("/mnt/data1t/pc_rebuild/web_cloud.json")
FIG = BASE / "figures"
WORK = Path("/mnt/data1t/film_mozdzek")

W, H = 1920, 1080
CX = 700                         # środek bryły w poziomie (panel jest po prawej)
VP_W = W                         # widok 3D na pełną szerokość; panel to nakładka
PANEL_X = W - 530
PANEL_W = 490

FONTS = Path("/usr/share/fonts/texlive-otf/ibm/plex")
F_SANS = FONTS / "IBMPlexSans-Regular.otf"
F_SANS_SB = FONTS / "IBMPlexSans-SemiBold.otf"
F_MONO = FONTS / "IBMPlexMono-Regular.otf"
F_MONO_SB = FONTS / "IBMPlexMono-SemiBold.otf"

BG = (13, 15, 19)                # #0d0f13 — viewport artefaktu
PANEL_BG = (22, 26, 32)
PANEL_EDGE = (42, 48, 58)
TRACK = (38, 43, 52)             # tło paska
INK = (238, 241, 245)
INK2 = (168, 176, 187)
INK3 = (121, 129, 141)

PALETA = ["#4d8ede", "#e08442", "#38ab86", "#d9a83c", "#c97ba6",
          "#5fb35f", "#8b7fd4", "#d96a63", "#6fb7c4"]
SUBTYPES = ["Aldoc_1", "Aldoc_2", "Aldoc_3", "Aldoc_4", "Aldoc_5",
            "Aldoc_6", "Aldoc_7", "Anti_Aldoc_1", "Anti_Aldoc_2"]


def hex2rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


COL = np.array([hex2rgb(c) for c in PALETA], dtype=np.float32)


def wczytaj(max_pts):
    z = np.load(SCENE, allow_pickle=True)
    i, j, k = z["i"].astype(np.float32), z["j"].astype(np.float32), z["k"].astype(np.float32)
    # UWAGA: `lobule` w npz jest 1-based (1..16), a `names` 0-based.
    # Rozstrzyga kolumna `lobule_name`; bez tej korekty kolory bryły są
    # przesunięte o jeden płacik względem panelu.
    lob = z["lobule"].astype(np.int16) - 1
    sub = z["subtype_idx"].astype(np.int16)
    names = [str(s) for s in z["names"]]
    lnm = z["lobule_name"]
    for v in np.unique(lob):
        got = {str(x) for x in lnm[lob == v][:200]}
        assert got == {names[v]}, f"mapowanie płacika {v}: {got} != {names[v]}"
    n = len(i)
    if n > max_pts:                                  # podpróbkowanie deterministyczne
        idx = np.random.default_rng(0).choice(n, max_pts, replace=False)
        idx.sort()
        i, j, k, lob, sub = i[idx], j[idx], k[idx], lob[idx], sub[idx]
    P = np.stack([k, i, j], 1)                       # k=ML(x), i=AP(y), j=DV(z)
    P -= P.mean(0)
    P /= np.abs(P).max()
    print(f"  punktów: {len(P)} (z {n}), płacików: {len(names)}")
    return P.astype(np.float32), lob, sub, names


def render3d(P, lob, sub, aktywny, waga, kat_y, elew, przesuw):
    """Jedna klatka widoku 3D. `waga` 0..1 = siła podświetlenia aktywnego płacika."""
    ky, ke = np.deg2rad(kat_y), np.deg2rad(elew)
    cy, sy = np.cos(ky), np.sin(ky)
    ce, se = np.cos(ke), np.sin(ke)
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], np.float32)
    Rx = np.array([[1, 0, 0], [0, ce, -se], [0, se, ce]], np.float32)
    Q = P @ Ry.T @ Rx.T

    d = 2.55                                         # odległość kamery
    zz = Q[:, 2] + d
    f = 1288.0
    px = (Q[:, 0] * f / zz) + CX + przesuw
    py = (-Q[:, 1] * f / zz) + H * 0.5

    ok = (px >= 0) & (px < VP_W - 2) & (py >= 0) & (py < H - 2) & (zz > 0.25)
    px, py, zz = px[ok], py[ok], zz[ok]
    lo, su = lob[ok], sub[ok]
    jest = (lo == aktywny)

    kol = COL[np.clip(su, 0, 8)]                     # kolor podtypu
    tlo_kol = np.array([96, 106, 122], np.float32)   # nieaktywna tkanka
    glebia = np.clip(1.30 - 0.40 * (zz - d + 1.0), 0.50, 1.20)[:, None]

    kol_buf = np.zeros((H, VP_W, 3), np.float32)
    kol_buf[:] = BG

    # 1) cała warstwa jako przygaszone tło anatomiczne
    b_tlo = np.clip(tlo_kol[None, :] * glebia * 0.62, 0, 255)
    kolejnosc = np.argsort(-zz)
    xi = px[kolejnosc].astype(np.int32)
    yi = py[kolejnosc].astype(np.int32)
    kol_buf[yi, xi] = b_tlo[kolejnosc]
    kol_buf[yi + 1, xi] = b_tlo[kolejnosc] * 0.80

    # 2) aktywny płacik na wierzchu, grubiej i w pełnym kolorze
    if jest.any() and waga > 0.01:
        ja = np.where(jest)[0]
        ja = ja[np.argsort(-zz[ja])]
        b = np.clip(kol[ja] * glebia[ja] * (0.55 + 0.45 * waga), 0, 255)
        ax_, ay_ = px[ja].astype(np.int32), py[ja].astype(np.int32)
        for dx, dy, sc in ((0, 0, 1.0), (1, 0, .92), (0, 1, .92), (1, 1, .80), (2, 1, .60), (1, 2, .60)):
            kol_buf[ay_ + dy, ax_ + dx] = b * sc
    return kol_buf.astype(np.uint8)


def czcionki():
    return {
        "hdr": ImageFont.truetype(str(F_SANS_SB), 15),
        "tytul": ImageFont.truetype(str(F_SANS_SB), 40),
        "meta": ImageFont.truetype(str(F_MONO), 17),
        "wiersz": ImageFont.truetype(str(F_MONO), 18),
        "wart": ImageFont.truetype(str(F_MONO_SB), 18),
        "stopka": ImageFont.truetype(str(F_SANS), 15),
    }


def panel(img, m, fnt, waga):
    warstwa = Image.new("RGBA", img.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(warstwa)
    x0, y0 = PANEL_X, 96
    x1, y1 = PANEL_X + PANEL_W, 96 + 560
    d.rounded_rectangle([x0, y0, x1, y1], 8, fill=PANEL_BG + (232,),
                        outline=PANEL_EDGE + (255,), width=1)

    pad = 26
    d.text((x0 + pad, y0 + 24), "L O B U L E", font=fnt["hdr"], fill=INK3)
    d.text((x0 + pad, y0 + 56), m.get("koz", m["name"]), font=fnt["tytul"], fill=INK)
    podpis = f"{m['name']} · n = {m['n']} cells · Aldoc+ {m['aldoc']*100:.1f}%"
    d.text((x0 + pad, y0 + 112), podpis, font=fnt["meta"], fill=INK3)

    pozycje = sorted(SUBTYPES, key=lambda s: -m["subs"].get(s, 0.0))
    ty = y0 + 158
    bx = x0 + pad + 168
    bw = PANEL_W - pad * 2 - 168 - 74
    for s in pozycje:
        v = m["subs"].get(s, 0.0) * 100
        d.text((x0 + pad, ty), s, font=fnt["wiersz"], fill=INK2)
        d.rounded_rectangle([bx, ty + 4, bx + bw, ty + 18], 3, fill=TRACK)
        if v > 0:
            szer = max(4, int(bw * v / 100 * waga))
            k = hex2rgb(PALETA[SUBTYPES.index(s)])
            d.rounded_rectangle([bx, ty + 4, bx + szer, ty + 18], 3, fill=k)
        d.text((x1 - pad, ty), f"{v:.1f}", font=fnt["wart"], fill=INK if v > 0 else INK3,
               anchor="ra")
        ty += 42
    img.paste(Image.alpha_composite(img.convert("RGBA"), warstwa).convert("RGB"), (0, 0))
    return img


def klatka(P, lob, sub, names, meta, t, n_klatek, fnt, per):
    """t = numer klatki. `per` = klatek na płacik."""
    idx = min(int(t // per), 15)
    faza = (t % per) / per
    waga = float(np.clip(faza / 0.18, 0, 1))          # wejście podświetlenia
    if faza > 0.88:                                   # wyjście
        waga = float(np.clip((1 - faza) / 0.12, 0, 1))

    post = t / max(1, n_klatek - 1)
    kat_y = -70 + 180 * post                          # 180° przez całość
    elew = 16 + 17 * np.sin(2 * np.pi * post * 1.25)  # kamera w górę/dół
    przesuw = 70 * np.sin(2 * np.pi * post * 0.7)     # i w bok

    nazwa = names[idx]
    m = next(v for v in meta.values() if v["name"] == nazwa)

    buf = render3d(P, lob, sub, idx, waga, kat_y, elew, przesuw)
    img = Image.new("RGB", (W, H), BG)
    img.paste(Image.fromarray(buf), (0, 0))
    panel(img, m, fnt, waga)

    d = ImageDraw.Draw(img)
    d.text((46, H - 52), f"{nazwa}  ·  {idx+1}/16", font=fnt["meta"], fill=INK3)
    return img


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=1200)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--points", type=int, default=650_000)
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--wznow", action="store_true",
                    help="pomiń klatki, które już są na dysku")
    a = ap.parse_args()

    WORK.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(exist_ok=True)
    meta = json.load(open(META))["meta"]
    P, lob, sub, names = wczytaj(a.points)
    fnt = czcionki()

    if a.test:
        per = a.frames / 16
        for n, t in enumerate([5, int(a.frames * 0.42), a.frames - 8]):
            klatka(P, lob, sub, names, meta, t, a.frames, fnt, per).save(
                WORK / f"test_{n}.png")
            print(f"  test_{n}.png  (klatka {t})")
        return

    per = a.frames / 16
    pominiete = 0
    for t in range(a.frames):
        cel = WORK / f"f_{t:05d}.png"
        if a.wznow and cel.exists() and cel.stat().st_size > 50_000:
            pominiete += 1
            continue
        klatka(P, lob, sub, names, meta, t, a.frames, fnt, per).save(cel)
        if t % 50 == 0:
            print(f"  {t}/{a.frames}", flush=True)
    if pominiete:
        print(f"  pominięto {pominiete} gotowych klatek")

    out = FIG / "film_mozdzek.mp4"
    subprocess.run([
        "ffmpeg", "-y", "-framerate", str(a.fps), "-i", str(WORK / "f_%05d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "17",
        "-preset", "slow", str(out)], check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"zapisano: {out}  ({out.stat().st_size/1e6:.1f} MB, "
          f"{a.frames} klatek, {a.frames/a.fps:.1f} s)")


if __name__ == "__main__":
    main()
