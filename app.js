// Desktop viewer: wraps OccultationScene with a camera, timeline and telemetry HUD.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { OccultationScene, loadOccultation, COLORS } from './sim.js';

const $ = (id) => document.getElementById(id);
const canvas = $('view');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true, logarithmicDepthBuffer: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.setClearColor(0x05070c, 1);
renderer.toneMapping = THREE.ACESFilmicToneMapping;

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(35, 1, 1, 2e6);
const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.minDistance = 50;
controls.maxDistance = 2e5;

const sim = new OccultationScene({ assetBase: 'assets/' });
scene.add(sim.root);

const stars = (() => {
  const n = 4000, pos = new Float32Array(n * 3), col = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) {
    const v = new THREE.Vector3().randomDirection().multiplyScalar(8e5);
    pos.set([v.x, v.y, v.z], 3 * i);
    const b = 0.35 + 0.65 * Math.random() ** 3;
    col.set([b, b, b * (0.9 + 0.2 * Math.random())], 3 * i);
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  return new THREE.Points(g, new THREE.PointsMaterial({ size: 1.6, sizeAttenuation: false, vertexColors: true }));
})();
scene.add(stars);

let index = [];
let data = null;
let t = 0;
let playing = true;
let speed = 20;
let view = 'overview';
let layers = {};        // current on/off state of the layer buttons
let follow = null;      // world point the camera is translated with each frame
let lastFollow = new THREE.Vector3();

// ---------------------------------------------------------------- loading
async function init() {
  await sim.init();
  index = await (await fetch('data/index.json')).json();
  const sel = $('occSelect');
  for (const e of index) {
    const opt = document.createElement('option');
    const lat = `${Math.abs(e.lat).toFixed(1)}°${e.lat >= 0 ? 'N' : 'S'}`;
    opt.value = e.id;
    opt.textContent = `#${String(e.id).padStart(3, '0')}  ${e.occUTC.replace('T', ' ').slice(0, 19)}  ${e.type}  ${lat} ${e.lon.toFixed(1)}°E`;
    sel.append(opt);
  }
  sel.addEventListener('change', () => selectOcc(+sel.value));
  const fromHash = +new URLSearchParams(location.hash.slice(1)).get('occ');
  await selectOcc(index.some((e) => e.id === fromHash) ? fromHash : index[0].id);
  $('loading').remove();
  requestAnimationFrame(frame);
}

async function selectOcc(id) {
  const entry = index.find((e) => e.id === id);
  $('occSelect').value = id;
  history.replaceState(null, '', `#occ=${id}`);
  data = await loadOccultation(`data/${entry.file}`);
  sim.setData(data);
  sim.setVisible(layers);
  const sza = data.atOcc.sza;
  $('occNote').textContent = [
    `${data.type === 'ingress' ? 'Ingress' : 'Egress'} · SZA ${sza.toFixed(0)}° (${sza < 90 ? 'dayside' : 'nightside'})`,
    data.forecast.Notes ? `note: ${data.forecast.Notes}` : null,
  ].filter(Boolean).join(' · ');
  t = Math.max(0, data.etStart - data.et0 - 30);
  sim.update(t);
  setView(view, true);
  drawChart();
}

// ---------------------------------------------------------------- camera views
function linkPlaneNormal() {
  // normal of the Mars-centre / MEX / TGO plane at the occultation moment: an edge-on view of the ray grazing the limb
  const i = Math.round((data.etOcc - data.et0) / data.step);
  const m = new THREE.Vector3().fromArray(data.mex, 3 * i).applyMatrix4(sim.frame.matrix);
  const g = new THREE.Vector3().fromArray(data.tgo, 3 * i).applyMatrix4(sim.frame.matrix);
  const n = new THREE.Vector3().crossVectors(m, g).normalize();
  // tilt slightly towards the sun so the planet reads as a lit sphere rather than a silhouette
  return n.addScaledVector(sim.sun.position.clone().normalize(), 0.25).normalize();
}

function setView(v, reset = false) {
  view = v;
  for (const b of document.querySelectorAll('#views button')) b.classList.toggle('on', b.dataset.view === v);
  if (!data) return;
  const n = linkPlaneNormal();
  const fit = (r) => r / Math.sin(THREE.MathUtils.degToRad(camera.fov / 2)) * 1.05;
  follow = null;
  if (v === 'overview') {
    let r = data.radii[0] * 1.2;
    const i0 = Math.max(0, Math.floor((data.etStart - data.et0) / data.step));
    const i1 = Math.min(data.n - 1, Math.ceil((data.etEnd - data.et0) / data.step));
    for (let i = i0; i <= i1; i++) {
      r = Math.max(r, Math.hypot(data.tgo[3 * i], data.tgo[3 * i + 1], data.tgo[3 * i + 2]));
      r = Math.max(r, Math.min(Math.hypot(data.mex[3 * i], data.mex[3 * i + 1], data.mex[3 * i + 2]), r * 2.2));
    }
    controls.target.set(0, 0, 0);
    camera.position.copy(n).multiplyScalar(fit(r));
  } else if (v === 'limb') {
    const tan = sim.tanMarker.visible ? sim.worldPos('tan') : sim.worldPos('tgo');
    controls.target.copy(tan);
    camera.position.copy(tan).addScaledVector(n, 2600).addScaledVector(tan.clone().normalize(), 600);
    follow = 'tan';
  } else if (v === 'mex') {
    follow = 'mexview';
  }
  if (follow && follow !== 'mexview') lastFollow.copy(sim.worldPos(follow === 'tan' ? 'tan' : follow));
  controls.update();
}

function updateCamera() {
  if (follow === 'mexview') {
    const mex = sim.worldPos('mex'), tgo = sim.worldPos('tgo');
    const dir = tgo.clone().sub(mex).normalize();
    const up = mex.clone().normalize();
    camera.position.copy(mex).addScaledVector(dir, -2000).addScaledVector(up, 450);
    controls.target.copy(tgo);
  } else if (follow === 'tan') {
    // the marker keeps its last valid position while the ray is away from the limb
    const p = sim.worldPos('tan');
    const d = p.clone().sub(lastFollow);
    camera.position.add(d);
    controls.target.add(d);
    lastFollow.copy(p);
  }
  sim.craft.mex.label.visible = layers.labels && follow !== 'mexview';
}

controls.addEventListener('start', () => {
  if (follow === 'mexview') { follow = null; setViewButton('free'); }
});
function setViewButton(v) {
  view = v;
  for (const b of document.querySelectorAll('#views button')) b.classList.toggle('on', b.dataset.view === v);
}

// ---------------------------------------------------------------- HUD
const fmt = (v, d, unit = '') => (v === null || v === undefined || Number.isNaN(v) ? '—' : `${v.toFixed(d)}${unit}`);
function hms(sec) {
  const s = Math.abs(sec);
  const m = Math.floor(s / 60);
  return `T${sec < 0 ? '−' : '+'}${String(m).padStart(2, '0')}:${(s % 60).toFixed(1).padStart(4, '0')}`;
}
function lst(h) {
  if (h === null) return '—';
  const hh = Math.floor(h), mm = Math.floor((h - hh) * 60);
  return `${String(hh).padStart(2, '0')}:${String(mm).padStart(2, '0')}`;
}

function updateHud(s) {
  $('utc').textContent = `${s.utc.toISOString().replace('T', ' ').slice(0, 21)} UTC`;
  $('trel').textContent = `${hms(s.tRel)} from ${data.type}`;
  const near = s.alt !== null && s.alt < 1500;
  $('telemetry').classList.toggle('blocked', s.blocked);
  $('linkState').textContent = s.blocked ? 'TGO behind Mars' : s.inTest ? 'Sounding atmosphere' : 'Link clear';
  $('vAlt').textContent = near ? fmt(s.alt, 1, ' km') : '—';
  $('vLat').textContent = near ? `${fmt(Math.abs(s.lat), 2, '°')}${s.lat >= 0 ? 'N' : 'S'}` : '—';
  $('vLon').textContent = near ? fmt(s.lon, 2, '°E') : '—';
  $('vSza').textContent = near ? fmt(s.sza, 1, '°') : '—';
  $('vLst').textContent = near ? lst(s.lst) : '—';
  $('vRange').textContent = fmt(s.range, 0, ' km');
  $('vDop').textContent = fmt(s.doppler / 1000, 3, ' kHz');
  $('vMexOff').textContent = fmt(s.mexOff, 1, '°');
  $('vTgoOff').textContent = fmt(s.tgoOff, 1, '°');
}

// ---------------------------------------------------------------- altitude chart / scrubber
const chart = $('chart');
const ALT_MIN = -200, ALT_MAX = 600;
function drawChart() {
  if (!data) return;
  const dpr = window.devicePixelRatio || 1;
  const w = chart.clientWidth, h = chart.clientHeight;
  chart.width = w * dpr; chart.height = h * dpr;
  const ctx = chart.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  const dur = sim.duration;
  const x = (tt) => (tt / dur) * w;
  const y = (a) => h - ((Math.min(Math.max(a, ALT_MIN), ALT_MAX) - ALT_MIN) / (ALT_MAX - ALT_MIN)) * h;

  ctx.fillStyle = 'rgba(255,255,255,0.07)';
  ctx.fillRect(x(data.etStart - data.et0), 0, x(data.etEnd - data.et0) - x(data.etStart - data.et0), h);
  ctx.fillStyle = 'rgba(122,168,255,0.16)';
  ctx.fillRect(0, y(200), w, y(100) - y(200));
  ctx.strokeStyle = 'rgba(255,255,255,0.25)';
  ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(0, y(0)); ctx.lineTo(w, y(0)); ctx.stroke();
  ctx.setLineDash([]);

  ctx.font = '11px Inter, sans-serif';
  ctx.fillStyle = 'rgba(200,210,230,0.6)';
  ctx.fillText('surface', 6, y(0) - 4);
  ctx.fillText('ray altitude (km)', 6, 12);
  for (const a of [200, 400]) {
    ctx.fillStyle = 'rgba(200,210,230,0.35)';
    ctx.fillText(`${a}`, w - 30, y(a) + 4);
    ctx.fillStyle = 'rgba(255,255,255,0.06)';
    ctx.fillRect(0, y(a), w, 1);
  }

  ctx.lineWidth = 2;
  let prev = null;
  for (let i = 0; i < data.n; i++) {
    const a = data.alt[i];
    if (a === null) { prev = null; continue; }
    const p = [x(i * data.step), y(a)];
    if (prev) {
      ctx.strokeStyle = a < 0 ? '#ff4d5e' : '#7cf2d6';
      ctx.beginPath(); ctx.moveTo(...prev); ctx.lineTo(...p); ctx.stroke();
    }
    prev = p;
  }
  const xo = x(data.etOcc - data.et0);
  ctx.strokeStyle = '#ff9a3c';
  ctx.beginPath(); ctx.moveTo(xo, 0); ctx.lineTo(xo, h); ctx.stroke();
  ctx.fillStyle = '#ff9a3c';
  ctx.fillText(data.type, xo + 5, 12);
  chartCursor();
}

const cursor = document.createElement('div');
Object.assign(cursor.style, { position: 'absolute', top: '0', width: '2px', background: '#fff', pointerEvents: 'none',
  boxShadow: '0 0 6px rgba(255,255,255,0.8)' });
$('timeline').style.position = 'fixed';
$('timeline').append(cursor);
function chartCursor() {
  if (!data) return;
  const r = chart.getBoundingClientRect(), pr = $('timeline').getBoundingClientRect();
  cursor.style.left = `${r.left - pr.left + (t / sim.duration) * r.width - 1}px`;
  cursor.style.top = `${r.top - pr.top}px`;
  cursor.style.height = `${r.height}px`;
}

let scrubbing = false;
function scrubTo(ev) {
  const r = chart.getBoundingClientRect();
  t = THREE.MathUtils.clamp((ev.clientX - r.left) / r.width, 0, 1) * sim.duration;
}
chart.addEventListener('pointerdown', (ev) => { scrubbing = true; chart.setPointerCapture(ev.pointerId); scrubTo(ev); });
chart.addEventListener('pointermove', (ev) => { if (scrubbing) scrubTo(ev); });
chart.addEventListener('pointerup', () => { scrubbing = false; });

// ---------------------------------------------------------------- controls
function setPlaying(p) { playing = p; $('playBtn').textContent = p ? 'Pause' : 'Play'; }
$('playBtn').addEventListener('click', () => {
  if (!playing && t >= sim.duration) t = 0;
  setPlaying(!playing);
});
$('speed').addEventListener('change', (e) => { speed = +e.target.value; });
for (const b of document.querySelectorAll('#views button')) b.addEventListener('click', () => setView(b.dataset.view));

function layerState() {
  const s = {};
  for (const b of document.querySelectorAll('#layers button')) s[b.dataset.layer] = b.classList.contains('on');
  return s;
}
for (const b of document.querySelectorAll('#layers button')) {
  b.addEventListener('click', () => {
    b.classList.toggle('on');
    const s = layers = layerState();
    stars.visible = s.stars;
    sim.setVisible(s);
  });
}
layers = layerState();
stars.visible = layers.stars;

function setScale() {
  const s = Math.round(10 ** +$('scScale').value);
  sim.setScScale(s);
  $('scScaleVal').textContent = s <= 1 ? 'true scale' : `×${s.toLocaleString('en-GB').replace(/,/g, ' ')}`;
}
$('scScale').addEventListener('input', setScale);

async function savePng(transparent) {
  const k = 3;
  const pr = renderer.getPixelRatio();
  const showStars = stars.visible;
  renderer.setPixelRatio(pr * k);
  sim.setResolution(canvas.clientWidth * pr * k, canvas.clientHeight * pr * k);
  if (transparent) { renderer.setClearColor(0x000000, 0); stars.visible = false; }
  renderer.render(scene, camera);
  const url = canvas.toDataURL('image/png');
  renderer.setPixelRatio(pr);
  renderer.setClearColor(0x05070c, 1);
  stars.visible = showStars;
  resize();
  const a = document.createElement('a');
  a.href = url;
  a.download = `occultation_${String(data.id).padStart(3, '0')}_${hms(sim.state.tRel).replace(/[:.]/g, '')}${transparent ? '_alpha' : ''}.png`;
  a.click();
}
$('savePng').addEventListener('click', () => savePng(false));
$('saveAlpha').addEventListener('click', () => savePng(true));
$('hideUi').addEventListener('click', () => document.body.classList.add('hide-ui'));

window.addEventListener('keydown', (e) => {
  if (e.target.tagName === 'SELECT') return;
  if (e.code === 'Space') { e.preventDefault(); setPlaying(!playing); }
  if (e.code === 'ArrowRight') t = Math.min(sim.duration, t + (e.shiftKey ? 60 : 10));
  if (e.code === 'ArrowLeft') t = Math.max(0, t - (e.shiftKey ? 60 : 10));
  if (e.code === 'KeyO') t = data.etOcc - data.et0;
  if (e.code === 'KeyH') document.body.classList.toggle('hide-ui');
});

// ---------------------------------------------------------------- loop
function resize() {
  const w = canvas.clientWidth, h = canvas.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  sim.setResolution(w * renderer.getPixelRatio(), h * renderer.getPixelRatio());
  drawChart();
}
window.addEventListener('resize', resize);

let last = performance.now();
function frame(now) {
  const dt = Math.min((now - last) / 1000, 0.1);
  last = now;
  if (playing && !scrubbing) {
    t += dt * speed;
    if (t >= sim.duration) { t = sim.duration; setPlaying(false); }
  }
  const s = sim.update(t);
  updateCamera();
  controls.update();
  updateHud(s);
  chartCursor();
  renderer.render(scene, camera);
  requestAnimationFrame(frame);
}

resize();
init().catch((err) => {
  console.error(err);
  $('loading').textContent = `Failed to load: ${err.message}`;
});
window.__occ = { sim, camera, controls, get t() { return t; }, set t(v) { t = v; }, setPlaying, setView, COLORS };
