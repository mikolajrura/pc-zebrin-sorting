# Sesja EBF2 vs EBF1 w CIEMNYM motywie (tlo sceny + paleta Qt6).
#
# URUCHOMIENIE - QT_QPA_PLATFORM=xcb jest OBOWIAZKOWE:
#   QT_QPA_PLATFORM=xcb ~/omics-data-highway/.venv/bin/python -m pymol scripts/41_ebf2_pymol_dark.pml
#
# Bez tego Qt6 startuje natywnie na Wayland (QT_QPA_PLATFORM=wayland;xcb w srodowisku),
# idzie przez EGL, dostaje kontekst "OpenGL ES 3.2", shader default sie nie kompiluje
# (gl_FrontColor zarezerwowane w ES) i SCENA 3D JEST PUSTA. Przez xcb/XWayland jest GLX
# i pelny OpenGL 4.6 - zweryfikowane 2026-09-05.

load references/struktury/AF-O08792-EBF2_mysz.pdb, EBF2
load references/struktury/AF-Q07802-EBF1_mysz.pdb, EBF1

# --- ciemne tlo sceny 3D (obejmuje tez wewnetrzne GUI PyMOL-a) ---
bg_color grey10
set ray_opaque_background, 1
set orthoscopic, 1
set depth_cue, 0
set ray_trace_fog, 0
set specular, 0.2
set cartoon_loop_radius, 0.25
set internal_gui_width, 260

hide everything
show cartoon, EBF2 or EBF1
cealign EBF2, EBF1

select dbd,      EBF2 and resi 34-244
select ipt,      EBF2 and resi 253-336
select epitope,  EBF2 and resi 413-550
select epi_core, EBF2 and resi 496-528
select znf,      EBF2 and resi 150-169
deselect

# --- kolory dobrane pod CIEMNE tlo (nie te z renderow na bialym) ---
color grey45, EBF1
set cartoon_transparency, 0.6, EBF1
color grey65, EBF2
color 0x3987e5, dbd
color 0x199e70, ipt
color 0xff7a45, epitope
color 0xffc300, epi_core
show spheres, znf and name CA
color 0x9085e9, znf
set sphere_scale, 0.45

# --- struktura drugorzedowa liczona, nie z nagłowka pliku ---
dss

# =====================================================================
# SCENY - przelaczanie klawiszami F1-F4 albo "scene <nazwa>, recall"
# =====================================================================

# F1: oba bialka nalozone (domyslny widok)
hide everything
show cartoon, EBF2 or EBF1
show spheres, znf and name CA
orient EBF2
turn x, -15
scene naloznie, store, message="EBF2 w kolorach + EBF1 szary przezroczysty (RMSD 0.89 A)"

# F2: sam EBF2 z domenami
hide everything
show cartoon, EBF2
show spheres, znf and name CA
orient EBF2
turn x, -15
scene domeny, store, message="niebieski=DNA-binding 34-244 | zielony=IPT/TIG 253-336 | szary 337-412=dimeryzacja (3 helisy) | pomaranczowy/zolty=epitop 413-550, 100% petla"

# F3: zblizenie na epitop, domena dimeryzacji przygaszona
hide everything
show cartoon, EBF2
color grey30, EBF2 and resi 337-412
orient epitope
zoom epitope, 4
scene epitop, store, message="epitop HPA003954 413-550: pLDDT 36.6, 0/138 reszt o pLDDT>=70 - tu przeciwcialo rozroznia EBF2 od EBF1 (62% tozsamosci)"

# F4: kolorowanie wg pewnosci predykcji AlphaFold
hide everything
show cartoon, EBF2
spectrum b, red_yellow_green, EBF2, 30, 95
orient EBF2
turn x, -15
scene plddt, store, message="pLDDT: czerwony=30 (chaos) -> zielony=95 (pewna struktura). Srednia 71.9, mediana 87.7"

scene naloznie, recall

# --- ciemna paleta Qt dla menu, panelu i pol tekstowych ---
python
try:
    from pymol.Qt import QtWidgets, QtGui
    app = QtWidgets.QApplication.instance()
    if app is None:
        print(" dark mode Qt: brak QApplication (tryb -cq?) - pomijam")
    else:
        app.setStyle("Fusion")
        C = QtGui.QColor
        R = QtGui.QPalette.ColorRole
        text = C(232, 231, 225)
        p = QtGui.QPalette()
        p.setColor(R.Window,          C(32, 32, 30))
        p.setColor(R.WindowText,      text)
        p.setColor(R.Base,            C(26, 26, 25))
        p.setColor(R.AlternateBase,   C(38, 38, 36))
        p.setColor(R.ToolTipBase,     C(38, 38, 36))
        p.setColor(R.ToolTipText,     text)
        p.setColor(R.Text,            text)
        p.setColor(R.Button,          C(45, 45, 43))
        p.setColor(R.ButtonText,      text)
        p.setColor(R.BrightText,      C(255, 122, 69))
        p.setColor(R.Link,            C(57, 135, 229))
        p.setColor(R.Highlight,       C(57, 135, 229))
        p.setColor(R.HighlightedText, C(16, 16, 15))
        g = QtGui.QPalette.ColorGroup.Disabled
        for role in (R.Text, R.ButtonText, R.WindowText):
            p.setColor(g, role, C(130, 129, 122))
        app.setPalette(p)
        app.setStyleSheet(
            "QMenu{background:#262624;color:#e8e7e1;border:1px solid #3a3a37;}"
            "QMenu::item:selected{background:#3987e5;color:#101010;}"
            "QMenuBar{background:#202020;color:#e8e7e1;}"
            "QMenuBar::item:selected{background:#3987e5;color:#101010;}"
            "QLineEdit,QPlainTextEdit,QTextEdit{background:#1a1a19;color:#e8e7e1;"
            "border:1px solid #3a3a37;}"
            "QTabBar::tab{background:#262624;color:#e8e7e1;padding:4px;}"
            "QTabBar::tab:selected{background:#3987e5;color:#101010;}"
            "QToolTip{background:#262624;color:#e8e7e1;border:1px solid #3a3a37;}")
        print(" dark mode Qt: zastosowany")
except Exception as e:
    print(" dark mode Qt NIE zadzialal:", e)
python end
