#!/usr/bin/env python3
"""
Czy komórki Ebf2-wysokie tworzą pasy parasagittalne, czy tylko gradient płatowy?

Wejście:  processed/abc_ebf2_komorki.csv  (16 626 PC z ABC MERFISH, współrzędne CCFv3)
Wyjście:  processed/abc_ebf2_pasy.csv     wynik per (sekcja, płacik)
          processed/abc_ebf2_domeny.csv   rozmiary spójnych domen

Konstrukcja testu — trzy zabezpieczenia przed potwierdzeniem samego siebie:

1. Dychotomizacja Ebf2 na hi/lo robiona ODDZIELNIE w każdej parze (sekcja, płacik).
   To usuwa gradient płatowy (26–84 %) i różnice jakości między sekcjami, więc test
   mierzy wyłącznie strukturę WEWNĄTRZ płacika.
2. Istotność z permutacji etykiet w obrębie tej samej pary — rozkład zerowy ma
   dokładnie ten sam skład i tę samą geometrię, różni się tylko przypisaniem.
3. KONTROLA NEGATYWNA: identyczny test dla Ebf1, który jest neutralny wobec osi
   Aldoc (log2FC +0.03) i płaski po płacikach (87–96 %). Jeśli Ebf1 też pokaże
   skupianie, mierzymy artefakt techniczny, nie biologię.

Test A — autokorelacja w płaszczyźnie sekcji (k najbliższych sąsiadów).
Test B — ciągłość między sąsiednimi sekcjami: pas parasagittalny biegnie wzdłuż osi
         przednio-tylnej, więc powinien być widoczny w tym samym miejscu (y,z)
         na sąsiedniej sekcji. Plama — nie.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix

P = Path(__file__).resolve().parent.parent / "processed"
GENY = ["Ebf2", "Ebf1", "Ebf3"]
K = 6            # sąsiadów w teście A
N_PERM = 500
MIN_N = 50       # minimalna liczność pary (sekcja, płacik)
RNG = np.random.default_rng(0)


def zgodnosc(xy, hi, k=K):
    """Odsetek par (komórka, jej k sąsiadów) o tej samej etykiecie."""
    n = len(hi)
    kk = min(k, n - 1)
    if kk < 1:
        return np.nan, None
    _, idx = cKDTree(xy).query(xy, k=kk + 1)
    sasiedzi = idx[:, 1:]
    return float((hi[:, None] == hi[sasiedzi]).mean()), sasiedzi


def test_A(d, gen):
    """Autokorelacja przestrzenna wewnątrz (sekcja, płacik), z permutacją."""
    wyniki = []
    for (sec, lob), g in d.groupby(["brain_section_label", "parcellation_substructure"]):
        if len(g) < MIN_N:
            continue
        v = g[gen].values
        prog = np.median(v)
        hi = v > prog if prog > 0 else v > 0
        if hi.mean() < 0.15 or hi.mean() > 0.85:   # zbyt niezbalansowane
            continue
        xy = g[["y_ccf", "z_ccf"]].values
        obs, sasiedzi = zgodnosc(xy, hi)
        if sasiedzi is None:
            continue
        perm = np.empty(N_PERM)
        for i in range(N_PERM):
            h = RNG.permutation(hi)
            perm[i] = (h[:, None] == h[sasiedzi]).mean()
        sd = perm.std()
        wyniki.append(dict(sekcja=sec, placik=lob, n=len(g), gen=gen,
                           frakcja_hi=float(hi.mean()), zgodnosc=obs,
                           zgodnosc_perm=float(perm.mean()),
                           z=float((obs - perm.mean()) / sd) if sd > 0 else np.nan))
    return pd.DataFrame(wyniki)


def test_B(d, gen, prom=0.15):
    """Ciągłość między sąsiednimi sekcjami w tym samym płaciku.

    Dla każdej komórki szukamy najbliższej komórki w sąsiedniej sekcji (w promieniu
    `prom` mm w płaszczyźnie y,z) i sprawdzamy zgodność etykiet. Permutacja tasuje
    etykiety w sekcji docelowej.
    """
    d = d.copy()
    for (sec, lob), g in d.groupby(["brain_section_label", "parcellation_substructure"]):
        v = g[gen].values
        prog = np.median(v)
        d.loc[g.index, "_hi"] = (v > prog if prog > 0 else v > 0).astype(float)
    sekcje = sorted(d["brain_section_label"].unique(),
                    key=lambda s: d.loc[d["brain_section_label"] == s, "x_ccf"].median())
    wyniki = []
    for a, b in zip(sekcje[:-1], sekcje[1:]):
        for lob in d["parcellation_substructure"].unique():
            ga = d[(d["brain_section_label"] == a) & (d["parcellation_substructure"] == lob)]
            gb = d[(d["brain_section_label"] == b) & (d["parcellation_substructure"] == lob)]
            if len(ga) < MIN_N or len(gb) < MIN_N:
                continue
            ha, hb = ga["_hi"].values, gb["_hi"].values
            if min(ha.mean(), hb.mean()) < 0.15 or max(ha.mean(), hb.mean()) > 0.85:
                continue
            dist, idx = cKDTree(gb[["y_ccf", "z_ccf"]].values).query(
                ga[["y_ccf", "z_ccf"]].values, k=1)
            ok = dist <= prom
            if ok.sum() < MIN_N:
                continue
            obs = float((ha[ok] == hb[idx[ok]]).mean())
            perm = np.array([(ha[ok] == RNG.permutation(hb)[idx[ok]]).mean()
                             for _ in range(N_PERM)])
            sd = perm.std()
            wyniki.append(dict(sekcja_a=a, sekcja_b=b, placik=lob, n_par=int(ok.sum()),
                               gen=gen, zgodnosc=obs, zgodnosc_perm=float(perm.mean()),
                               z=float((obs - perm.mean()) / sd) if sd > 0 else np.nan))
    return pd.DataFrame(wyniki)


def domeny(d, gen, prom=0.12):
    """Spójne skupiska komórek hi w obrębie (sekcja, płacik)."""
    rozmiary = []
    for (sec, lob), g in d.groupby(["brain_section_label", "parcellation_substructure"]):
        if len(g) < MIN_N:
            continue
        v = g[gen].values
        prog = np.median(v)
        hi = v > prog if prog > 0 else v > 0
        xy = g[["y_ccf", "z_ccf"]].values[hi]
        if len(xy) < 3:
            continue
        pary = cKDTree(xy).query_pairs(prom, output_type="ndarray")
        if len(pary) == 0:
            rozmiary.extend([1] * len(xy))
            continue
        m = coo_matrix((np.ones(len(pary)), (pary[:, 0], pary[:, 1])),
                       shape=(len(xy), len(xy)))
        _, etyk = connected_components(m, directed=False)
        rozmiary.extend(np.bincount(etyk).tolist())
    return np.array(rozmiary)


def main():
    d = pd.read_csv(P / "abc_ebf2_komorki.csv")
    print(f"komórek: {len(d)}, sekcji: {d['brain_section_label'].nunique()}")

    print(f"\n{'='*72}\nTEST A — skupianie w płaszczyźnie sekcji (k={K}, {N_PERM} permutacji)\n{'='*72}")
    wszystkie_A = []
    for gen in GENY:
        r = test_A(d, gen)
        wszystkie_A.append(r)
        w = r["n"] / r["n"].sum()
        print(f"\n{gen}:  par (sekcja,płacik) = {len(r)}, komórek = {int(r['n'].sum())}")
        print(f"  zgodność sąsiadów:  {(r['zgodnosc']*w).sum()*100:.1f} %   "
              f"permutacja: {(r['zgodnosc_perm']*w).sum()*100:.1f} %")
        print(f"  mediana z = {r['z'].median():+.2f}   par z z>3: {(r['z']>3).sum()}/{len(r)}"
              f"   par z z>2: {(r['z']>2).sum()}/{len(r)}")
    A = pd.concat(wszystkie_A, ignore_index=True)
    A.to_csv(P / "abc_ebf2_pasy.csv", index=False)

    print(f"\n{'='*72}\nTEST B — ciągłość między sąsiednimi sekcjami (odstęp ~200 µm)\n{'='*72}")
    for gen in GENY:
        r = test_B(d, gen)
        if len(r) == 0:
            print(f"\n{gen}: brak par spełniających kryteria")
            continue
        w = r["n_par"] / r["n_par"].sum()
        print(f"\n{gen}:  par sekcji×płacik = {len(r)}, porównań = {int(r['n_par'].sum())}")
        print(f"  zgodność:  {(r['zgodnosc']*w).sum()*100:.1f} %   "
              f"permutacja: {(r['zgodnosc_perm']*w).sum()*100:.1f} %")
        print(f"  mediana z = {r['z'].median():+.2f}   par z z>3: {(r['z']>3).sum()}/{len(r)}")

    print(f"\n{'='*72}\nROZMIARY DOMEN (spójne skupiska hi, promień 0.12 mm)\n{'='*72}")
    wier = []
    for gen in GENY:
        r = domeny(d, gen)
        duze = r[r >= 3]
        print(f"\n{gen}:  domen razem {len(r)}, ≥3 komórek: {len(duze)}")
        if len(duze):
            print(f"  mediana {np.median(duze):.0f}, p75 {np.percentile(duze,75):.0f}, "
                  f"maks {duze.max()}, komórek w domenach ≥3: {duze.sum()}")
        wier.append(pd.DataFrame({"gen": gen, "rozmiar": r}))
    pd.concat(wier, ignore_index=True).to_csv(P / "abc_ebf2_domeny.csv", index=False)
    print(f"\nzapisano: {P}/abc_ebf2_pasy.csv, abc_ebf2_domeny.csv")


if __name__ == "__main__":
    main()
