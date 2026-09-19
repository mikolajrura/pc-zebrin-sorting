#!/usr/bin/env bash
# Odtworzenie WERSJI 1 modelu 3D `figures/purkinje3d.html`.
#
# Kontekst: v1 opublikowano 2026-09-04 18:22:28 UTC, v2 — 2026-09-05 00:24:28 UTC
# z TEJ SAMEJ ścieżki pliku, więc v2 nadpisała artefakt pod tym samym URL-em,
# a skrypt generujący nadpisał plik lokalny. Lokalna kopia v1 przestała istnieć.
#
# Odzyskano 2026-09-18: generator wyciągnięty z transkryptu sesji
# ~/.claude/projects/-home-mikolajrurad/635bdc59-…jsonl (linia 1398), dane wejściowe
# `/mnt/data1t/pc_rebuild/web_cloud.json` (1 219 120 B, mtime 2026-09-04 20:20) —
# nietknięte, bo v2 czytała dodatkowo `web_registered.json` i nie ruszała tego pliku.
#
# Wyjście: figures/purkinje3d_v1.html (NIE purkinje3d.html — żeby nie nadpisać v2).
# Opublikowane jako osobny artefakt: https://claude.ai/artifact/FiX9nPVWm6yDy1vYQbsQaw
#
# v1 zawiera: 2 279 886 wokseli warstwy, 2,28 mm³, 16 634 komórek Purkinjego,
# 16 płacików; tryby „Frakcja Aldoc+" / „Dominujący podtyp" / „Płacik".
# NIE zawiera 244 438 komórek Stereo-seq — te doszły dopiero w v2.

set -euo pipefail
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cat > "$TMP"/mk_html.py <<'PYEOF'
import json, os
OUT="/mnt/data1t/pc_rebuild"
D=json.load(open(f"{OUT}/web_cloud.json"))
DEST="/home/mikolajrurad/omics-data-highway/ichb-purkinje-analysis/pc-zebrin-sorting/figures/purkinje3d_v1.html"

HTML = r"""<title>Warstwa Purkinjego myszy</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<style>
:root{
  --ground:#f7f7f5; --panel:#ffffff; --edge:#e0e0dc; --edge-soft:#eceae6;
  --ink:#14161a; --ink-2:#4b5058; --ink-3:#787f88;
  --accent:#2a6fd6; --accent-soft:#e8f0fc;
  --viewport:#0d0f13; --viewport-edge:#232833;
  --warn:#9a5b12;
  --r:3px;
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){
  --ground:#0f1115; --panel:#161a20; --edge:#2a303a; --edge-soft:#222831;
  --ink:#eef1f5; --ink-2:#a8b0bb; --ink-3:#79818d;
  --accent:#6aa5f0; --accent-soft:#1b2735;
  --viewport:#080a0d; --viewport-edge:#242a34;
  --warn:#d9a05b;
}}
:root[data-theme="dark"]{
  --ground:#0f1115; --panel:#161a20; --edge:#2a303a; --edge-soft:#222831;
  --ink:#eef1f5; --ink-2:#a8b0bb; --ink-3:#79818d;
  --accent:#6aa5f0; --accent-soft:#1b2735;
  --viewport:#080a0d; --viewport-edge:#242a34;
  --warn:#d9a05b;
}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);
  font-family:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
  font-size:14px;line-height:1.5}
.wrap{max-width:1420px;margin:0 auto;padding:26px 22px 40px}
header{display:flex;flex-wrap:wrap;gap:18px 26px;align-items:flex-end;justify-content:space-between;
  padding-bottom:16px;border-bottom:1px solid var(--edge)}
h1{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:26px;margin:0 0 5px;
  letter-spacing:-0.01em;text-wrap:balance}
.sub{color:var(--ink-2);margin:0;max-width:62ch}
.stats{display:flex;gap:22px;flex-wrap:wrap}
.stat .v{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:19px;font-weight:500;
  font-variant-numeric:tabular-nums;display:block;line-height:1.2}
.stat .l{font-size:10.5px;text-transform:uppercase;letter-spacing:.09em;color:var(--ink-3);
  display:block;margin-top:3px}
main{display:grid;grid-template-columns:minmax(0,1fr) 330px;gap:18px;margin-top:18px}
@media(max-width:960px){main{grid-template-columns:minmax(0,1fr)}}
.viewport{position:relative;background:var(--viewport);border:1px solid var(--viewport-edge);
  border-radius:var(--r);overflow:hidden;min-height:520px;aspect-ratio:4/3;cursor:grab}
.viewport.drag{cursor:grabbing}
canvas{display:block;width:100%;height:100%}
.hud{position:absolute;left:12px;bottom:11px;font-family:"IBM Plex Mono",monospace;font-size:11px;
  color:#7d8794;pointer-events:none;letter-spacing:.02em}
.tip{position:absolute;right:12px;top:11px;font-family:"IBM Plex Mono",monospace;font-size:11px;
  color:#7d8794;pointer-events:none;text-align:right}
.rail{display:flex;flex-direction:column;gap:14px;min-width:0}
.card{background:var(--panel);border:1px solid var(--edge);border-radius:var(--r);padding:13px 14px}
.card h2{font-size:10.5px;text-transform:uppercase;letter-spacing:.1em;color:var(--ink-3);
  margin:0 0 10px;font-weight:600}
.modes{display:flex;flex-direction:column;gap:5px}
.modes button{all:unset;cursor:pointer;padding:7px 10px;border-radius:2px;font-size:13px;
  border:1px solid transparent;display:flex;justify-content:space-between;align-items:center;gap:8px}
.modes button:hover{background:var(--edge-soft)}
.modes button[aria-pressed="true"]{background:var(--accent-soft);border-color:var(--accent);color:var(--ink)}
.modes button:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
.modes .k{font-family:"IBM Plex Mono",monospace;font-size:10.5px;color:var(--ink-3)}
.ramp{height:9px;border-radius:1px;margin:8px 0 5px;
  background:linear-gradient(90deg,#12203a,#1d3c6e,#2a63b4,#4d8ede,#8fbcf0,#cfe2fa)}
.rampL{display:flex;justify-content:space-between;font-family:"IBM Plex Mono",monospace;
  font-size:10.5px;color:var(--ink-3);font-variant-numeric:tabular-nums}
.lobs{max-height:290px;overflow-y:auto;display:flex;flex-direction:column;gap:1px;margin:-3px -4px}
.lob{all:unset;cursor:pointer;display:grid;grid-template-columns:52px 1fr 40px;gap:8px;
  align-items:center;padding:5px 6px;border-radius:2px}
.lob:hover{background:var(--edge-soft)}
.lob[aria-pressed="true"]{background:var(--accent-soft)}
.lob:focus-visible{outline:2px solid var(--accent);outline-offset:-1px}
.lob .n{font-family:"IBM Plex Mono",monospace;font-size:11.5px;color:var(--ink-2)}
.lob .bar{height:7px;background:var(--edge-soft);border-radius:1px;overflow:hidden}
.lob .bar i{display:block;height:100%;background:var(--accent)}
.lob .p{font-family:"IBM Plex Mono",monospace;font-size:11px;text-align:right;
  font-variant-numeric:tabular-nums;color:var(--ink-2)}
.detail .big{font-family:"IBM Plex Serif",Georgia,serif;font-size:17px;font-weight:600;margin:0 0 2px}
.detail .meta{font-family:"IBM Plex Mono",monospace;font-size:11px;color:var(--ink-3);margin:0 0 10px;
  font-variant-numeric:tabular-nums}
.comp{display:flex;flex-direction:column;gap:4px}
.crow{display:grid;grid-template-columns:78px 1fr 38px;gap:7px;align-items:center;font-size:11.5px}
.crow .cn{font-family:"IBM Plex Mono",monospace;font-size:10.5px;color:var(--ink-2);
  white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
.crow .cb{height:6px;background:var(--edge-soft);border-radius:1px;overflow:hidden}
.crow .cb i{display:block;height:100%}
.crow .cv{font-family:"IBM Plex Mono",monospace;font-size:10.5px;text-align:right;
  color:var(--ink-3);font-variant-numeric:tabular-nums}
.notes{margin-top:20px;border-top:1px solid var(--edge);padding-top:16px;
  display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:16px 26px}
.note h3{font-size:10.5px;text-transform:uppercase;letter-spacing:.1em;color:var(--ink-3);margin:0 0 6px}
.note p{margin:0 0 7px;color:var(--ink-2);font-size:12.5px;line-height:1.55;max-width:60ch}
.note b{color:var(--ink);font-weight:600}
code{font-family:"IBM Plex Mono",monospace;font-size:11.5px;background:var(--edge-soft);
  padding:1px 4px;border-radius:2px}
.flag{border-left:2px solid var(--warn);padding-left:10px}
.flag h3{color:var(--warn)}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
</style>

<div class="wrap">
<header>
  <div>
    <h1>Warstwa Purkinjego myszy</h1>
    <p class="sub">Geometria z atlasu CCFv3a przy 10&nbsp;µm, pokolorowana składem fenotypowym
    zmierzonym w snRNA-seq. Przeciągnij, żeby obrócić; kółkiem przybliżasz.</p>
  </div>
  <div class="stats">
    <div class="stat"><span class="v">2 279 886</span><span class="l">wokseli warstwy</span></div>
    <div class="stat"><span class="v">2,28 mm³</span><span class="l">objętość</span></div>
    <div class="stat"><span class="v">16 634</span><span class="l">komórek Purkinjego</span></div>
    <div class="stat"><span class="v">16</span><span class="l">płacików</span></div>
  </div>
</header>

<main>
  <div class="viewport" id="vp">
    <canvas id="cv"></canvas>
    <div class="hud" id="hud">—</div>
    <div class="tip">przeciągnij&nbsp;— obrót · kółko&nbsp;— zoom<br>klik płacika po prawej&nbsp;— podświetlenie</div>
  </div>

  <div class="rail">
    <div class="card">
      <h2>Kolor punktów</h2>
      <div class="modes" id="modes">
        <button data-m="aldoc" aria-pressed="true">Frakcja Aldoc+ <span class="k">ciągła</span></button>
        <button data-m="sub" aria-pressed="false">Dominujący podtyp <span class="k">9 kat.</span></button>
        <button data-m="lob" aria-pressed="false">Płacik <span class="k">16 kat.</span></button>
      </div>
      <div id="legend">
        <div class="ramp"></div>
        <div class="rampL"><span>0,07 — sam Aldoc−</span><span>1,00 — sam Aldoc+</span></div>
      </div>
    </div>

    <div class="card">
      <h2>Płaciki — udział Aldoc+</h2>
      <div class="lobs" id="lobs"></div>
    </div>

    <div class="card detail" id="detail"></div>
  </div>
</main>

<section class="notes">
  <div class="note">
    <h3>Skąd bierze się geometria</h3>
    <p>Atlas <b>CCFv3a</b> (Piluso i wsp., <i>Imaging Neuroscience</i> 2025) definiuje w hierarchii
    warstwę Purkinjego, ale <b>wydany wolumen adnotacji ma dla niej zero wokseli</b> — sprawdzone
    przelotem przez wszystkie 1 290 480 000 wokseli. Powód jest fizyczny: warstwa ma ok. 20&nbsp;µm,
    a przy 25&nbsp;µm jest podwokselowa.</p>
    <p>Policzyłem ją więc geometrycznie, jako <b>granicę warstwy ziarnistej i drobinowej</b> —
    tam z definicji leżą ciała komórek Purkinjego. Kontrola spójności: granica stanowi
    <code>7,0–11,8 %</code> objętości warstwy ziarnistej w każdym z 16 płacików osobno, mimo że
    różnią się one objętością czterdziestokrotnie.</p>
  </div>
  <div class="note">
    <h3>Skąd bierze się kolor</h3>
    <p>Dane snRNA-seq: <b>Kozareva i wsp.</b>, <i>Nature</i> 2021, GEO <code>GSE165371</code>.
    Macierz przebudowana od zera ze strumieniowego wycięcia z pliku <code>.mtx.gz</code> —
    <b>16 634 komórki</b>, walidacja per komórka (16 538 zgodnych co do liczby niezerowych genów).</p>
    <p>Mapowanie 16 płacików Kozarevy na struktury atlasu jest <b>pełne, 16/16</b>.
    Klasyfikator 9 podtypów na 30 składowych głównych: dokładność <code>0,917</code>,
    zbalansowana <code>0,929</code>; po redukcji do osi Aldoc+/− <code>0,974</code>.</p>
  </div>
  <div class="note flag">
    <h3>Czego ten obraz nie pokazuje</h3>
    <p>Kolor jest przypisany <b>per płacik, nie per komórka</b>. Kozareva daje etykietę płacika,
    nie współrzędną, więc wszystkie punkty w obrębie jednego płacika mają ten sam kolor.</p>
    <p><b>To nie są pasy zebriny.</b> Pasy to wąskie pasma parasagittalne wewnątrz płacików —
    drobniejsze niż rozdzielczość tych danych. Widać tu gradient płatowy, który jest realny
    i zmierzony, ale to inny poziom organizacji.</p>
  </div>
</section>
</div>

<script id="payload" type="application/json">__DATA__</script>
<script>
(function(){
const D=JSON.parse(document.getElementById('payload').textContent);
const b64=s=>{const b=atob(s),u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);return u;};
const posBuf=new Int16Array(b64(D.pos).buffer);
const lob=b64(D.lob);
const N=D.n;
const meta=D.meta, subs=D.subtypes;

// --- palety ---
const rampStops=[[0x12,0x20,0x3a],[0x1d,0x3c,0x6e],[0x2a,0x63,0xb4],[0x4d,0x8e,0xde],[0x8f,0xbc,0xf0],[0xcf,0xe2,0xfa]];
function ramp(t){t=Math.max(0,Math.min(1,t));const x=t*(rampStops.length-1),i=Math.floor(x),f=x-i;
  const a=rampStops[i],b=rampStops[Math.min(i+1,rampStops.length-1)];
  return [(a[0]+(b[0]-a[0])*f)/255,(a[1]+(b[1]-a[1])*f)/255,(a[2]+(b[2]-a[2])*f)/255];}
const CAT=["#4d8ede","#e08442","#38ab86","#d9a83c","#c97ba6","#5fb35f","#8b7fd4","#d96a63",
           "#6fb7c4","#b0894a","#7f96b5","#c59acb","#63a878","#cf9a6a","#8fa2d8","#a8b36a"];
const hex2rgb=h=>[parseInt(h.slice(1,3),16)/255,parseInt(h.slice(3,5),16)/255,parseInt(h.slice(5,7),16)/255];
const lobIds=Object.keys(meta).map(Number).sort((a,b)=>a-b);
const subIdx={}; subs.forEach((s,i)=>subIdx[s]=i);

// --- three ---
const cv=document.getElementById('cv'), vp=document.getElementById('vp');
const renderer=new THREE.WebGLRenderer({canvas:cv,antialias:true,alpha:false});
renderer.setPixelRatio(Math.min(devicePixelRatio,2));
const scene=new THREE.Scene(); scene.background=new THREE.Color(0x0d0f13);
const cam=new THREE.PerspectiveCamera(35,1,0.1,400000);
const geo=new THREE.BufferGeometry();
const pos=new Float32Array(N*3);
for(let i=0;i<N*3;i++) pos[i]=posBuf[i];
geo.setAttribute('position',new THREE.BufferAttribute(pos,3));
const col=new Float32Array(N*3);
geo.setAttribute('color',new THREE.BufferAttribute(col,3));
const mat=new THREE.PointsMaterial({size:260,vertexColors:true,sizeAttenuation:true});
const cloud=new THREE.Points(geo,mat); scene.add(cloud);

let mode='aldoc', sel=null;
function paint(){
  const ca=geo.attributes.color.array;
  for(let i=0;i<N;i++){
    const L=lob[i], m=meta[L]; let c;
    if(!m){c=[0.3,0.3,0.3];}
    else if(mode==='aldoc') c=ramp(m.aldoc);
    else if(mode==='lob')   c=hex2rgb(CAT[(L-1)%CAT.length]);
    else                    c=hex2rgb(CAT[subIdx[m.dominant]%CAT.length]);
    let d=1.0;
    if(sel!==null) d = (L===sel)?1.0:0.16;
    ca[i*3]=c[0]*d; ca[i*3+1]=c[1]*d; ca[i*3+2]=c[2]*d;
  }
  geo.attributes.color.needsUpdate=true;
}

// orbit
let rx=0.30, ry=-0.55, dist=95000, dragging=false, px=0, py=0;
function place(){
  cam.position.set(dist*Math.cos(rx)*Math.sin(ry), dist*Math.sin(rx), dist*Math.cos(rx)*Math.cos(ry));
  cam.lookAt(0,0,0);
}
vp.addEventListener('pointerdown',e=>{dragging=true;px=e.clientX;py=e.clientY;vp.classList.add('drag');vp.setPointerCapture(e.pointerId);});
vp.addEventListener('pointermove',e=>{if(!dragging)return;
  ry-=(e.clientX-px)*0.0065; rx+=(e.clientY-py)*0.0065;
  rx=Math.max(-1.5,Math.min(1.5,rx)); px=e.clientX; py=e.clientY; place(); draw();});
vp.addEventListener('pointerup',e=>{dragging=false;vp.classList.remove('drag');});
vp.addEventListener('wheel',e=>{e.preventDefault();dist*=(e.deltaY>0?1.09:0.917);
  dist=Math.max(22000,Math.min(260000,dist));place();draw();},{passive:false});

function resize(){const r=vp.getBoundingClientRect();
  renderer.setSize(r.width,r.height,false);cam.aspect=r.width/r.height;cam.updateProjectionMatrix();}
function draw(){renderer.render(scene,cam);
  document.getElementById('hud').textContent=
    `${N.toLocaleString('pl-PL')} punktów · podpróbka z 2 279 886 · tryb: ${
      mode==='aldoc'?'frakcja Aldoc+':mode==='sub'?'dominujący podtyp':'płacik'}`;}
addEventListener('resize',()=>{resize();draw();});

// UI
document.getElementById('modes').addEventListener('click',e=>{
  const b=e.target.closest('button'); if(!b)return;
  mode=b.dataset.m;
  [...e.currentTarget.children].forEach(x=>x.setAttribute('aria-pressed',x===b));
  document.getElementById('legend').style.display = mode==='aldoc'?'':'none';
  paint(); draw();});

const lobsEl=document.getElementById('lobs');
lobIds.forEach(id=>{const m=meta[id];
  const b=document.createElement('button'); b.className='lob'; b.dataset.id=id;
  b.setAttribute('aria-pressed','false');
  b.innerHTML=`<span class="n">${m.koz}</span><span class="bar"><i style="width:${(m.aldoc*100).toFixed(0)}%"></i></span><span class="p">${(m.aldoc*100).toFixed(0)}%</span>`;
  lobsEl.appendChild(b);});
lobsEl.addEventListener('click',e=>{const b=e.target.closest('.lob'); if(!b)return;
  const id=Number(b.dataset.id);
  sel = (sel===id)?null:id;
  [...lobsEl.children].forEach(x=>x.setAttribute('aria-pressed', Number(x.dataset.id)===sel));
  detail(sel); paint(); draw();});

function detail(id){
  const el=document.getElementById('detail');
  if(id===null||!meta[id]){
    el.innerHTML='<h2>Szczegóły płacika</h2><p style="margin:0;color:var(--ink-3);font-size:12.5px">Wybierz płacik z listy, żeby zobaczyć jego skład podtypowy i podświetlić go w modelu.</p>';
    return;}
  const m=meta[id];
  const rows=Object.entries(m.subs).sort((a,b)=>b[1]-a[1])
    .map(([k,v])=>`<div class="crow"><span class="cn">${k.replace('Purkinje_','')}</span>
      <span class="cb"><i style="width:${(v*100).toFixed(1)}%;background:${CAT[subIdx[k]%CAT.length]}"></i></span>
      <span class="cv">${(v*100).toFixed(1)}</span></div>`).join('');
  el.innerHTML=`<h2>Szczegóły płacika</h2>
    <p class="big">${m.koz}</p>
    <p class="meta">${m.name} · n = ${m.n.toLocaleString('pl-PL')} komórek · Aldoc+ ${(m.aldoc*100).toFixed(1)}%</p>
    <div class="comp">${rows}</div>`;
}
detail(null);
resize(); place(); paint(); draw();
})();
</script>
"""
HTML = HTML.replace("__DATA__", json.dumps(D))
os.makedirs(os.path.dirname(DEST), exist_ok=True)
open(DEST,"w",encoding="utf-8").write(HTML)
print("zapisano", DEST, os.path.getsize(DEST), "B")
PYEOF
/home/mikolajrurad/omics-data-highway/.venv/bin/python "$TMP"/mk_html.py