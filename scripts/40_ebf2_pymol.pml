# Sesja PyMOL: EBF2 vs EBF1 (modele AlphaFold) z podswietlonym epitopem HPA003954.
# NIE JEST to dokowanie - przeciwciala poliklonalne nie maja struktury, a epitop
# lezy w regionie nieuporzadkowanym (0/138 reszt o pLDDT>=70).
# Uruchomienie: pymol -cq scripts/40_ebf2_pymol.pml

set ray_opaque_background, 0
set orthoscopic, 1
set cartoon_transparency, 0
bg_color white

load references/struktury/AF-O08792-EBF2_mysz.pdb, EBF2
load references/struktury/AF-Q07802-EBF1_mysz.pdb, EBF1

hide everything
show cartoon, EBF2 or EBF1
set cartoon_loop_radius, 0.25

# superpozycja po strukturze (nie po sekwencji)
cealign EBF2, EBF1

# --- selekcje ---
select dbd,      EBF2 and resi 34-244
select ipt,      EBF2 and resi 253-336
select epitope,  EBF2 and resi 413-550
select epi_core, EBF2 and resi 496-528
select znf,      EBF2 and resi 150-169
deselect

# ================= WIDOK 1: EBF2 pokolorowany pLDDT =================
hide everything
show cartoon, EBF2
spectrum b, red_yellow_green, EBF2, 30, 95
orient EBF2
turn x, -15
png figures/pymol_1_ebf2_plddt.png, width=1600, height=1100, dpi=160, ray=1

# ================= WIDOK 2: domeny + epitop =================
hide everything
show cartoon, EBF2
color grey70, EBF2
color 0x2a78d6, dbd
color 0x1baf7a, ipt
color 0xeb6834, epitope
color 0xc23c0a, epi_core
show spheres, znf and name CA
color 0x4a3aa7, znf
set sphere_scale, 0.45
orient EBF2
turn x, -15
png figures/pymol_2_ebf2_domeny.png, width=1600, height=1100, dpi=160, ray=1

# ================= WIDOK 3: naloznie EBF1 (szary) na EBF2 =================
hide everything
show cartoon, EBF2 or EBF1
color grey80, EBF1
set cartoon_transparency, 0.55, EBF1
color 0x2a78d6, dbd
color 0x1baf7a, ipt
color 0xeb6834, epitope
color grey60, EBF2 and not (dbd or ipt or epitope)
orient EBF2 or EBF1
turn x, -15
png figures/pymol_3_naloznie.png, width=1600, height=1100, dpi=160, ray=1

# ================= WIDOK 4: zblizenie na epitop =================
hide everything
show cartoon, EBF2
color grey80, EBF2
color 0xeb6834, epitope
color 0xc23c0a, epi_core
orient epitope
zoom epitope, 4
png figures/pymol_4_epitop_zblizenie.png, width=1600, height=1100, dpi=160, ray=1

# --- sesja do otwarcia interaktywnie ---
hide everything
show cartoon, EBF2 or EBF1
color grey80, EBF1
set cartoon_transparency, 0.55, EBF1
color grey70, EBF2
color 0x2a78d6, dbd
color 0x1baf7a, ipt
color 0xeb6834, epitope
color 0xc23c0a, epi_core
orient EBF2
save processed/EBF2_vs_EBF1_epitop.pse

print "--- RMSD i statystyki ---"
python
from pymol import cmd
r = cmd.align("EBF1","EBF2")
print("align EBF1->EBF2: RMSD %.2f A na %d atomach, %d reszt w dopasowaniu" % (r[0], r[1], r[4]))
for nm, sel in [("domena DNA 34-244","dbd"),("IPT/TIG 253-336","ipt"),
                ("epitop 413-550","epitope"),("rdzen 496-528","epi_core")]:
    n = cmd.count_atoms(sel+" and name CA")
    b = cmd.get_model(sel+" and name CA").get_residues()
    bs = [a.b for a in cmd.get_model(sel+" and name CA").atom]
    print("  %-22s %3d reszt, sredni pLDDT %.1f" % (nm, n, sum(bs)/len(bs)))
python end
