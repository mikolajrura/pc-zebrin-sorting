#!/bin/bash
# Odpala PyMOL z sesja EBF2 vs EBF1 (ciemny motyw, sceny F1-F4).
#
#   ./scripts/42_pymol_ebf2.sh            # GUI do obracania
#   ./scripts/42_pymol_ebf2.sh --render   # tylko PNG, bez okna (omija sterownik GL)
#
# QT_QPA_PLATFORM=xcb jest OBOWIAZKOWE na tej maszynie: w srodowisku jest
# QT_QPA_PLATFORM="wayland;xcb", wiec Qt6 startuje natywnie na Wayland, idzie przez
# EGL i dostaje kontekst "OpenGL ES 3.2". Shadery PyMOL-a sa pisane pod desktop GL,
# nie kompiluja sie (gl_FrontColor zarezerwowane w ES) i SCENA 3D JEST PUSTA.
# Przez xcb/XWayland jest GLX i pelny OpenGL 4.6 - zweryfikowane 2026-09-05.
set -euo pipefail

P="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="$HOME/omics-data-highway/.venv/bin/python"

[ -x "$PY" ] || { echo "BRAK interpretera: $PY" >&2; exit 1; }
"$PY" -c "import pymol" 2>/dev/null || {
  echo "BRAK pymol w venv. Zainstaluj:" >&2
  echo "  uv pip install --python $PY pymol-open-source" >&2; exit 1; }

cd "$P"
for f in references/struktury/AF-O08792-EBF2_mysz.pdb references/struktury/AF-Q07802-EBF1_mysz.pdb; do
  [ -f "$f" ] || { echo "BRAK struktury: $f" >&2; exit 1; }
done

if [ "${1:-}" = "--render" ]; then
  # tryb bez GUI - wlasny ray tracer, nie dotyka OpenGL, dziala zawsze
  exec "$PY" -m pymol -cq scripts/40_ebf2_pymol.pml
fi

echo "PyMOL: sceny F1=naloznie  F2=domeny  F3=epitop  F4=pLDDT"
exec env QT_QPA_PLATFORM=xcb "$PY" -m pymol scripts/41_ebf2_pymol_dark.pml
