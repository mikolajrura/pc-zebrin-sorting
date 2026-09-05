"""Pobranie poziomu 10 um. NIE dotykamy .annotation (to by wczytalo 4.81 GiB)."""
import time, os, subprocess
t0=time.time()
from brainglobe_atlasapi import BrainGlobeAtlas
print("pobieram ccfv3augmented_mouse_10um (bez wczytywania do RAM) ...", flush=True)
a=BrainGlobeAtlas("ccfv3augmented_mouse_10um", check_latest=False)
print(f"metadane gotowe w {time.time()-t0:.0f} s")
print("rozdzielczosc:", a.resolution, "ksztalt wg metadanych:", a.metadata.get('shape'))
Z=os.path.expanduser("~/.brainglobe/brainglobe-atlasapi/annotation-sets/ccfv3augmented_mouse-annotation/3_0")
print("\nrozmiar katalogu adnotacji po pobraniu:")
print(subprocess.run(["du","-sh",Z],capture_output=True,text=True).stdout.strip())
