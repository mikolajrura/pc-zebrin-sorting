# Figura do raportu: EBF2 na CZARNYM tle, epitop przeciwciala wyrozniony.
# Render headless (-cq) - nie dotyka sterownika GL, dziala zawsze.
#   ~/omics-data-highway/.venv/bin/python -m pymol -cq scripts/45_fig_epitop_czarne.pml

load references/struktury/AF-O08792-EBF2_mysz.pdb, EBF2
bg_color black
set ray_opaque_background, 1
set orthoscopic, 1
set depth_cue, 0
set specular, 0.15
set cartoon_loop_radius, 0.3
set ray_shadows, 0
set antialias, 2
dss

select dbd,      resi 34-244
select ipt,      resi 253-336
select dimer,    resi 337-412
select epitope,  resi 413-550
select epi_core, resi 496-528
deselect

hide everything
show cartoon, EBF2
color grey40,   EBF2
color 0x3987e5, dbd
color 0x199e70, ipt
color 0x6f6f6f, dimer
color 0xff7a45, epitope
color 0xffc300, epi_core

orient EBF2
turn x, -15
png figures/raport_ebf2_epitop.png, width=2000, height=1250, dpi=200, ray=1

# drugi kadr: zblizenie na epitop
orient epitope
zoom epitope, 5
png figures/raport_ebf2_epitop_zoom.png, width=2000, height=1250, dpi=200, ray=1
