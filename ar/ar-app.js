// AR viewer: the OccultationScene from ../sim.js anchored to an ARToolKit pattern marker via AR.js (three.js build).
//
// URL options (all optional):
//   occ=6           forecast row to show
//   patt=markers/patt.hiro   pattern file of the printed marker
//   mars=1.0        Mars radius in marker widths
//   lift=0.3        Mars centre height in front of the marker, in marker widths
//   focus=0.7       what sits over the marker: 0 = Mars' centre, 1 = the limb point where the ray grazes
//   speed=20        playback speed (× real time)
//   sc=30000        spacecraft size exaggeration
//   smooth=1        pose smoothing strength (2 = steadier but laggier, 0.5 = snappier, 0 = off)
//   res=960         long side (px) of the image ARToolKit analyses; higher = steadier corners, more CPU
//   debug           show ARToolKit's detection canvas
import * as THREE from 'three';
import { ArToolkitSource, ArToolkitContext, ArMarkerControls } from './lib/ar-threex.mjs';
import { OccultationScene, loadOccultation } from '../sim.js';

const params = new URLSearchParams(location.search);
const cfg = {
  occ: +(params.get('occ') ?? 6),
  patt: params.get('patt') ?? 'markers/patt.hiro',
  mars: +(params.get('mars') ?? 1.0),
  lift: +(params.get('lift') ?? 0.3),
  focus: +(params.get('focus') ?? 0.7),
  speed: +(params.get('speed') ?? 20),
  sc: +(params.get('sc') ?? 30000),
  smooth: +(params.get('smooth') ?? 1),
  res: +(params.get('res') ?? 960),
  debug: params.has('debug'),
};
const PAD = 30;          // seconds shown either side of the test window
const END_HOLD = 2.5;    // seconds to hold the final frame before looping

const $ = (id) => document.getElementById(id);

// ---------------------------------------------------------------- three.js
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setClearColor(0x000000, 0);
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.domElement.id = 'gl';
document.body.prepend(renderer.domElement);

const scene = new THREE.Scene();
const camera = new THREE.Camera();      // projection comes from the ARToolKit calibration
scene.add(camera);

const markerRoot = new THREE.Group();   // raw pose written by ArMarkerControls
const anchor = new THREE.Group();       // smoothed pose that the content hangs from
anchor.visible = false;
scene.add(markerRoot, anchor);

// ---------------------------------------------------------------- pose smoothing
// One Euro filter (Casiez et al., CHI 2012): a low-pass whose cutoff rises with speed, so the model is held
// steady while the phone is nearly still (where jitter shows) but still follows deliberate movement promptly.
class OneEuro {
  constructor(minCutoff, beta, dCutoff = 1) {
    Object.assign(this, { minCutoff, beta, dCutoff, x: null, dx: null });
  }
  static alpha(cutoff, dt) { return 1 / (1 + 1 / (2 * Math.PI * cutoff * dt)); }
  reset() { this.x = null; }
  filter(x, dt) {
    // phones can deliver two frames with the same timestamp; dt = 0 would turn the state into NaN for good
    dt = Math.max(dt, 1 / 240);
    if (!this.x?.every(Number.isFinite)) { this.x = x.slice(); this.dx = x.map(() => 0); return this.x; }
    const ad = OneEuro.alpha(this.dCutoff, dt);
    this.dx = this.dx.map((d, i) => d + ad * ((x[i] - this.x[i]) / dt - d));
    const a = OneEuro.alpha(this.minCutoff + this.beta * Math.hypot(...this.dx), dt);
    this.x = this.x.map((v, i) => v + a * (x[i] - v));
    return this.x;
  }
}

const pose = {
  // positions are in marker widths (camera space), rotations as quaternion components
  pos: new OneEuro(1.2 / cfg.smooth, 0.8),
  rot: new OneEuro(0.5 / cfg.smooth, 1.2),
  lostFor: 0,
  jumpFrames: 0,
  p: new THREE.Vector3(), q: new THREE.Quaternion(), s: new THREE.Vector3(),
};
const HOLD_LOST = 0.6;     // s to keep showing the last pose after the marker drops out
const FLIP_DEG = 20;       // a sudden rotation larger than this ...
const FLIP_FRAMES = 6;     // ... must persist this many frames before it is believed (square-marker pose flips)

function updatePose(dt) {
  if (!markerRoot.visible) {
    pose.lostFor += dt;
    if (pose.lostFor > HOLD_LOST) { anchor.visible = false; pose.pos.reset(); pose.rot.reset(); }
    return;
  }
  pose.lostFor = 0;
  markerRoot.matrix.decompose(pose.p, pose.q, pose.s);
  if (cfg.smooth <= 0 || !anchor.visible) {
    pose.pos.reset(); pose.rot.reset();
  } else if (pose.rot.x) {
    const prev = new THREE.Quaternion().fromArray(pose.rot.x);
    if (prev.dot(pose.q) < 0) pose.q.set(-pose.q.x, -pose.q.y, -pose.q.z, -pose.q.w);  // same hemisphere
    if (THREE.MathUtils.radToDeg(prev.angleTo(pose.q)) > FLIP_DEG) {
      if (++pose.jumpFrames < FLIP_FRAMES) { anchor.visible = true; return; }   // ignore a one-off flip
      pose.pos.reset(); pose.rot.reset();                                       // genuine: snap to it
    }
    pose.jumpFrames = 0;
  }
  anchor.position.fromArray(pose.pos.filter(pose.p.toArray(), dt));
  anchor.quaternion.fromArray(pose.rot.filter(pose.q.toArray(), dt)).normalize();
  anchor.visible = true;
}

const stage = new THREE.Group();        // turns the occultation geometry to face out of the poster
anchor.add(stage);
const sim = new OccultationScene({ assetBase: '../assets/', modelTexturePath: '../assets/models/mobile/' });
stage.add(sim.root);

// ---------------------------------------------------------------- content
let data = null;
let t = 0, tStart = 0, tEnd = 0, hold = 0;
let playing = true;

async function loadContent() {
  await sim.init();
  const index = await (await fetch('../data/index.json')).json();
  const entry = index.find((e) => e.id === cfg.occ) ?? index[0];
  data = await loadOccultation(`../data/${entry.file}`);
  sim.setData(data);
  sim.setVisible({ orbits: true, iono: true, boresights: false, labels: true, trace: true });
  sim.setLabelSize(0.05);
  sim.setScScale(cfg.sc);
  for (const m of sim.lineMaterials) m.linewidth *= 1.6;   // thicker lines read better on a phone
  tStart = Math.max(0, data.etStart - data.et0 - PAD);
  tEnd = Math.min(sim.duration, data.etEnd - data.et0 + PAD);
  t = tStart;
  orientStage();
  $('date').textContent = `${data.occUTC.slice(0, 10)} · ${data.type}`;
  $('occMark').style.left = `${((data.etOcc - data.et0 - tStart) / (tEnd - tStart)) * 100}%`;
}

// Marker axes (AR.js): +Y out of the marker towards the viewer, -Z towards the top edge of the printed image.
// Point the link-plane normal (the edge-on "overview" direction, on the sunlit side) out of the poster and
// keep Mars' north pole as close to the top of the marker as possible.
function orientStage() {
  sim.update(data.etOcc - data.et0);
  const i = Math.round((data.etOcc - data.et0) / data.step);
  const mex = new THREE.Vector3().fromArray(data.mex, 3 * i).applyMatrix4(sim.frame.matrix);
  const tgo = new THREE.Vector3().fromArray(data.tgo, 3 * i).applyMatrix4(sim.frame.matrix);
  const sun = sim.sun.position.clone().normalize();
  const out = new THREE.Vector3().crossVectors(mex, tgo).normalize();
  if (out.dot(sun) < 0) out.negate();
  out.addScaledVector(sun, 0.25).normalize();
  const north = new THREE.Vector3(0, 1, 0).projectOnPlane(out).normalize();
  const z = north.clone().negate();
  const x = new THREE.Vector3().crossVectors(out, z);
  stage.quaternion.setFromRotationMatrix(new THREE.Matrix4().makeBasis(x, out, z)).invert();
  const scale = cfg.mars / data.radii[0];
  sim.setRootScale(1);
  // centre the composition on the limb where the ray grazes, not on Mars' centre (cfg.focus of the way there)
  sim.root.updateMatrixWorld(true);
  const tan = sim.root.worldToLocal(sim.tanMarker.getWorldPosition(new THREE.Vector3()));
  const shift = tan.multiplyScalar(cfg.focus * scale).applyQuaternion(stage.quaternion);
  stage.position.set(-shift.x, cfg.lift - shift.y, -shift.z);
  sim.setRootScale(scale);
}

// ---------------------------------------------------------------- AR.js
const source = new ArToolkitSource({ sourceType: 'webcam', sourceWidth: 1280, sourceHeight: 960 });
let context = null;
let trackedAspect = 0;

function startTracking() {
  // size ARToolKit's analysis canvas to the real video aspect (AR.js assumes 4:3 landscape, which misplaces the
  // overlay on phones held upright) and pick the matching calibration
  const video = source.domElement;
  const w = video.videoWidth, h = video.videoHeight;
  const k = Math.min(1, cfg.res / Math.max(w, h));
  context?.dispose();
  context = new ArToolkitContext({
    cameraParametersUrl: h > w ? 'markers/camera_para_portrait.dat' : 'markers/camera_para.dat',
    detectionMode: 'mono',
    maxDetectionRate: 30,
    canvasWidth: Math.round(w * k),
    canvasHeight: Math.round(h * k),
    debug: cfg.debug,
  });
  context.init(() => {
    camera.projectionMatrix.copy(context.getProjectionMatrix());
    camera.projectionMatrixInverse.copy(camera.projectionMatrix).invert();
  });
  new ArMarkerControls(context, markerRoot, { type: 'pattern', patternUrl: cfg.patt, changeMatrixMode: 'modelViewMatrix' });
  trackedAspect = w / h;
}

function onResize() {
  const video = source.domElement;
  if (!source.ready || !video?.videoWidth) return;
  source.onResizeElement();                     // video covers the window (cropped, aspect kept)
  const s = video.style;
  const el = renderer.domElement;
  // overlay exactly the same box as the video so the projection lines up
  Object.assign(el.style, { width: s.width, height: s.height, marginLeft: s.marginLeft, marginTop: s.marginTop });
  const w = parseFloat(s.width), h = parseFloat(s.height);
  renderer.setSize(w, h, false);
  sim.setResolution(w * renderer.getPixelRatio(), h * renderer.getPixelRatio());
  if (Math.abs(video.videoWidth / video.videoHeight - trackedAspect) > 0.01) startTracking();
}

function startCamera() {
  return new Promise((resolve, reject) => {
    source.init(() => {
      const video = source.domElement;
      const ready = () => { onResize(); resolve(); };
      if (video.videoWidth) ready(); else video.addEventListener('loadedmetadata', ready, { once: true });
      video.addEventListener('resize', onResize);  // camera stream rotates with the phone
    }, (err) => reject(new Error(err?.message || 'camera unavailable')));
  });
}
window.addEventListener('resize', () => setTimeout(onResize, 300));

// ---------------------------------------------------------------- HUD
function updateHud(s, visible) {
  $('hint').classList.toggle('hidden', visible);
  $('hud').classList.toggle('blocked', s.blocked);
  $('state').textContent = s.blocked ? 'TGO behind Mars' : s.inTest ? 'Sounding the atmosphere' : 'Link clear';
  $('alt').textContent = s.alt !== null && s.alt < 1500 ? `${s.alt.toFixed(0)} km` : '—';
  const sec = Math.abs(s.tRel);
  $('trel').textContent = `T${s.tRel < 0 ? '−' : '+'}${String(Math.floor(sec / 60)).padStart(2, '0')}:${String(Math.floor(sec % 60)).padStart(2, '0')} from ${data.type}${playing ? '' : ' · paused'}`;
  $('prog').style.width = `${((t - tStart) / (tEnd - tStart)) * 100}%`;
}

renderer.domElement.addEventListener('click', () => { playing = !playing; });
$('hud').addEventListener('click', () => { playing = !playing; });

// ---------------------------------------------------------------- loop
let last = performance.now();
function frame(now) {
  const dt = THREE.MathUtils.clamp((now - last) / 1000, 0, 0.1);
  last = now;
  if (context?.arController && source.ready) context.update(source.domElement);
  updatePose(dt);

  if (playing) {
    if (t < tEnd) t = Math.min(tEnd, t + dt * cfg.speed);
    else if ((hold += dt) > END_HOLD) { t = tStart; hold = 0; }
  }
  updateHud(sim.update(t), anchor.visible);
  renderer.render(scene, camera);
  requestAnimationFrame(frame);
}

function fail(msg) {
  $('errorMsg').innerHTML = msg;
  $('error').classList.add('show');
}

(async () => {
  if (!window.isSecureContext) {
    fail('The camera needs a secure page.<br>Open this over <b>https://</b> (or localhost).');
    return;
  }
  try {
    await Promise.all([loadContent(), startCamera()]);
    $('hintSub').textContent = 'Use the marker on the poster';
    requestAnimationFrame(frame);
  } catch (err) {
    console.error(err);
    fail(`Could not start: ${err.message}`);
  }
})();

window.__ar = { sim, anchor, markerRoot, stage, camera, cfg, get context() { return context; }, get t() { return t; }, set t(v) { t = v; }, setPlaying(p) { playing = p; } };
