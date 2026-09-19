#!/usr/bin/env python3
"""
Rozkład przestrzenny Ebf2 w komórkach móżdżku — ABC Atlas MERFISH-C57BL6J-638850.

Wejście:
  /mnt/data1t/abc_merfish/C57BL6J-638850-log2.h5ad   CSR 4 334 174 x 550, nnz 613 476 175
  <scratch>/abc_cerebellum.csv                        komórki parcellation_division==CB + CCF

Wyjście (processed/):
  abc_ebf2_per_subclass.csv   ekspresja wybranych genów per typ komórki móżdżku
  abc_ebf2_per_lobule.csv     Purkinje: ekspresja per płacik
  abc_ebf2_komorki.csv        pojedyncze komórki Purkinjego + współrzędne CCF

Reguła §5 CLAUDE.md: cała macierz gęsta f32 to 8.87 GiB — NIE wczytujemy.
Czytamy CSR blokami wierszy prosto z h5py, od razu redukując do wybranych kolumn.
Szczyt zużycia = jeden blok, kontrolowany przez BLOK.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import h5py

H5 = Path("/mnt/data1t/abc_merfish/C57BL6J-638850-log2.h5ad")
SCRATCH = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
META = SCRATCH / "abc_cerebellum.csv"
OUT = Path(__file__).resolve().parent.parent / "processed"

GENY = ["Ebf1", "Ebf2", "Ebf3"]
KONTROLE = ["Calb1", "Calb2", "Ppp1r17", "Sorcs3", "Nxph1", "Aqp4", "Grm1", "Eomes"]
BLOK = 200_000  # wierszy na raz

MAPA = {  # CCFv3 -> nomenklatura Kozarevy
    "LING": "I", "CENT2": "II", "CENT3": "III", "CUL4, 5": "CUL", "DEC": "VI",
    "FOTU": "VII", "PYR": "VIII", "UVU": "IX", "NOD": "X", "SIM": "SIM",
    "ANcr1": "AN1", "ANcr2": "AN2", "PRM": "PRM", "COPY": "COP", "PFL": "PF", "FL": "F",
}


def main():
    meta = pd.read_csv(META)
    print(f"komórek móżdżku w metadanych: {len(meta)}")

    f = h5py.File(H5, "r")
    n_wierszy, n_genow = f["X"].attrs["shape"]
    symbole = np.array([s.decode() if isinstance(s, bytes) else s
                        for s in f["var/gene_symbol"][:]])
    etykiety = np.array([s.decode() if isinstance(s, bytes) else s
                         for s in f["obs/cell_label"][:]])
    print(f"macierz: {n_wierszy} x {n_genow}, komórek z etykietą: {len(etykiety)}")

    brak = [g for g in GENY + KONTROLE if g not in symbole]
    if brak:
        print(f"UWAGA: brak w panelu 550: {brak}")
    wanted = [g for g in GENY + KONTROLE if g in symbole]
    kol_idx = np.array([int(np.where(symbole == g)[0][0]) for g in wanted])
    print(f"geny odczytywane ({len(wanted)}): {wanted}")

    # dopasowanie komórek móżdżku do wierszy macierzy
    poz = pd.Series(np.arange(len(etykiety)), index=etykiety)
    wiersz = poz.reindex(meta["cell_label"].astype(str)).values
    ok = ~pd.isna(wiersz)
    print(f"dopasowanych do macierzy: {int(ok.sum())} / {len(meta)}")
    meta = meta.loc[ok].reset_index(drop=True)
    wiersz = wiersz[ok].astype(np.int64)

    chce = np.zeros(n_wierszy, dtype=bool)
    chce[wiersz] = True
    # docelowa pozycja każdego wiersza macierzy w tabeli wynikowej
    docelowa = np.full(n_wierszy, -1, dtype=np.int64)
    docelowa[wiersz] = np.arange(len(wiersz))

    print(f"wynik: {len(wiersz)} x {len(wanted)} = "
          f"{len(wiersz)*len(wanted)*4/2**20:.1f} MiB (f32)")

    indptr = f["X/indptr"][:].astype(np.int64)
    data_ds, ind_ds = f["X/data"], f["X/indices"]
    X = np.zeros((len(wiersz), len(wanted)), dtype=np.float32)
    # mapa: nr kolumny w macierzy -> nr kolumny w wyniku (-1 = pomijamy)
    kol_map = np.full(n_genow, -1, dtype=np.int32)
    kol_map[kol_idx] = np.arange(len(kol_idx))

    n_blokow = -(-n_wierszy // BLOK)
    for b in range(n_blokow):
        r0, r1 = b * BLOK, min((b + 1) * BLOK, n_wierszy)
        if not chce[r0:r1].any():
            continue
        p0, p1 = indptr[r0], indptr[r1]
        dane = data_ds[p0:p1]
        kol = ind_ds[p0:p1]
        # rozwiń numery wierszy dla każdego niezerowego elementu
        licznosci = np.diff(indptr[r0:r1 + 1])
        wiersze_nnz = np.repeat(np.arange(r0, r1), licznosci)
        maska = chce[wiersze_nnz] & (kol_map[kol] >= 0)
        if maska.any():
            X[docelowa[wiersze_nnz[maska]], kol_map[kol[maska]]] = dane[maska]
        print(f"  blok {b+1}/{n_blokow}: wiersze {r0}–{r1}, nnz {p1-p0}", flush=True)
    f.close()

    for i, g in enumerate(wanted):
        meta[g] = X[:, i]

    # --- per subclass ---
    rows = []
    for sc, d in meta.groupby("subclass"):
        r = {"subclass": sc, "n": len(d)}
        for g in wanted:
            r[f"{g}_mean"] = d[g].mean()
            r[f"{g}_pct"] = 100.0 * (d[g] > 0).mean()
        rows.append(r)
    per_sc = pd.DataFrame(rows).sort_values("n", ascending=False)
    per_sc.to_csv(OUT / "abc_ebf2_per_subclass.csv", index=False)

    # --- Purkinje per płacik ---
    pc = meta[meta["subclass"].str.contains("Purkinje", na=False)].copy()
    pc["lobule"] = pc["parcellation_substructure"].map(MAPA)
    rows = []
    for lob, d in pc.dropna(subset=["lobule"]).groupby("lobule"):
        r = {"lobule": lob, "n": len(d)}
        for g in wanted:
            r[f"{g}_mean"] = d[g].mean()
            r[f"{g}_pct"] = 100.0 * (d[g] > 0).mean()
        rows.append(r)
    per_lob = pd.DataFrame(rows).sort_values("n", ascending=False)
    per_lob.to_csv(OUT / "abc_ebf2_per_lobule.csv", index=False)

    kol_out = ["cell_label", "brain_section_label", "subclass", "supertype", "cluster",
               "x_ccf", "y_ccf", "z_ccf", "parcellation_substructure"] + wanted
    pc[[c for c in kol_out if c in pc.columns]].to_csv(OUT / "abc_ebf2_komorki.csv", index=False)

    fmt = lambda v: f"{v:.3f}"
    print("\n=== EBF1/2/3 PER TYP KOMÓRKI MÓŻDŻKU (top 18) ===")
    kols = ["subclass", "n"] + [f"{g}_pct" for g in GENY] + [f"{g}_mean" for g in GENY]
    print(per_sc[[c for c in kols if c in per_sc]].head(18).to_string(index=False, float_format=fmt))

    print("\n=== PURKINJE: EBF2 PER PŁACIK ===")
    kols = ["lobule", "n"] + [f"{g}_pct" for g in GENY] + ["Ebf2_mean", "Calb1_pct"]
    print(per_lob[[c for c in kols if c in per_lob]].to_string(index=False, float_format=fmt))

    print(f"\nzapisano do {OUT}/: abc_ebf2_per_subclass.csv, abc_ebf2_per_lobule.csv, "
          f"abc_ebf2_komorki.csv")


if __name__ == "__main__":
    main()
