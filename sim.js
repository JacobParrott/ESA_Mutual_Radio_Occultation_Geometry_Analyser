// OccultationScene - self-contained three.js scene for one MEX -> TGO radio occultation.
// Everything lives under `scene.root`, so the same object can be dropped into a desktop scene or an AR marker.
// Data comes from export_occultations.py (km, Mars-centred, +Z = Mars north pole).
import * as THREE from 'three';
import { TDSLoader } from 'three/addons/loaders/TDSLoader.js';
import { Line2 } from 'three/addons/lines/Line2.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';

export const COLORS = {
  mex: 0x6fb7ff,
  tgo: 0xffc857,
  link: 0x7cf2d6,
  blocked: 0xff4d5e,
  trace: 0xff9a3c,
  iono: 0x7aa8ff,
};

// Meshes as Cosmographia configures them: model units -> km scale, and the part each sits on
const MODELS = {
  mex: { scale: 1, bus: 'mex_bus.3ds', parts: { sa_py: 'mex_sa.3ds', sa_my: 'mex_sa.3ds' } },
  tgo: { scale: 0.001, bus: 'tgo_bus.3ds', parts: { sa_pz: 'tgo_samz.3ds', sa_mz: 'tgo_sapz.3ds', hga: 'tgo_hga.3ds' } },
};

export async function loadOccultation(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`${url}: ${res.status}`);
  return res.json();
}

function fatLine(points, color, width, opts = {}) {
  const geom = new LineGeometry();
  geom.setPositions(points);
  const mat = new LineMaterial({ color, linewidth: width, transparent: true, opacity: opts.opacity ?? 1,
    dashed: !!opts.dashed, dashSize: opts.dashSize ?? 60, gapSize: opts.gapSize ?? 40, depthWrite: false });
  const line = new Line2(geom, mat);
  line.computeLineDistances();
  return line;
}

function textSprite(text, color) {
  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');
  const font = '600 64px "Inter", "Helvetica Neue", Arial, sans-serif';
  ctx.font = font;
  canvas.width = Math.ceil(ctx.measureText(text).width) + 32;
  canvas.height = 96;
  ctx.font = font;
  ctx.fillStyle = `#${new THREE.Color(color).getHexString()}`;
  ctx.shadowColor = 'rgba(0,0,0,0.9)';
  ctx.shadowBlur = 12;
  ctx.textBaseline = 'middle';
  ctx.fillText(text, 16, 50);
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, sizeAttenuation: false, depthWrite: false,
    depthTest: false }));
  sprite.userData.aspect = canvas.width / canvas.height;
  sprite.center.set(-0.08, -0.25);
  sprite.renderOrder = 10;
  return sprite;
}

function glowTexture(color) {
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = 128;
  const ctx = canvas.getContext('2d');
  const c = new THREE.Color(color);
  const rgb = `${Math.round(c.r * 255)},${Math.round(c.g * 255)},${Math.round(c.b * 255)}`;
  const g = ctx.createRadialGradient(64, 64, 0, 64, 64, 64);
  g.addColorStop(0, 'rgba(255,255,255,1)');
  g.addColorStop(0.18, `rgba(${rgb},1)`);
  g.addColorStop(0.45, `rgba(${rgb},0.35)`);
  g.addColorStop(1, `rgba(${rgb},0)`);
  ctx.fillStyle = g;
  ctx.fillRect(0, 0, 128, 128);
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

function ionosphereMaterial(color) {
  return new THREE.ShaderMaterial({
    uniforms: { glowColor: { value: new THREE.Color(color) }, power: { value: 2.5 }, strength: { value: 0.6 } },
    vertexShader: /* glsl */`
      #include <common>
      #include <logdepthbuf_pars_vertex>
      varying vec3 vNormal; varying vec3 vView;
      void main() {
        vec4 mv = modelViewMatrix * vec4(position, 1.0);
        vNormal = normalize(normalMatrix * normal);
        vView = normalize(-mv.xyz);
        gl_Position = projectionMatrix * mv;
        #include <logdepthbuf_vertex>
      }`,
    fragmentShader: /* glsl */`
      #include <logdepthbuf_pars_fragment>
      uniform vec3 glowColor; uniform float power; uniform float strength;
      varying vec3 vNormal; varying vec3 vView;
      void main() {
        #include <logdepthbuf_fragment>
        float rim = pow(1.0 - abs(dot(vNormal, vView)), power);
        gl_FragColor = vec4(glowColor, rim * strength);
      }`,
    transparent: true,
    blending: THREE.AdditiveBlending,
    depthWrite: false,
  });
}

export class OccultationScene {
  // modelTexturePath: where the 3DS models' textures are read from (AR uses downscaled copies for phones)
  constructor({ assetBase = 'assets/', modelTexturePath = null } = {}) {
    this.assetBase = assetBase;
    this.modelTexturePath = modelTexturePath ?? `${assetBase}models/`;
    this.root = new THREE.Group();
    this.root.name = 'OccultationScene';
    // OCC frame (+Z north) -> three.js (+Y up)
    this.frame = new THREE.Group();
    this.frame.rotation.x = -Math.PI / 2;
    this.root.add(this.frame);

    this.lineMaterials = [];
    this.scScale = 30000;
    this.data = null;
    this.state = {};
    this._q = [new THREE.Quaternion(), new THREE.Quaternion()];
  }

  async init() {
    const texLoader = new THREE.TextureLoader();
    const marsTex = await texLoader.loadAsync(`${this.assetBase}textures/2k_mars.jpg`);
    marsTex.colorSpace = THREE.SRGBColorSpace;
    marsTex.anisotropy = 8;

    // Mars: body-fixed node rotated per sample; mesh rotated so the sphere's +Y pole lies on IAU_MARS +Z
    this.marsNode = new THREE.Group();
    this.frame.add(this.marsNode);
    this.marsMesh = new THREE.Mesh(new THREE.SphereGeometry(1, 192, 96),
      new THREE.MeshStandardMaterial({ map: marsTex, roughness: 1, metalness: 0 }));
    this.marsMesh.rotation.x = Math.PI / 2;
    this.marsNode.add(this.marsMesh);

    this.iono = new THREE.Mesh(new THREE.SphereGeometry(1, 128, 64), ionosphereMaterial(COLORS.iono));
    this.iono.rotation.x = Math.PI / 2;
    this.marsNode.add(this.iono);

    // tangent point: constant on-screen size so it reads at every zoom level
    this.tanMarker = new THREE.Sprite(new THREE.SpriteMaterial({ map: glowTexture(COLORS.trace), sizeAttenuation: false,
      depthTest: false, depthWrite: false, blending: THREE.AdditiveBlending }));
    this.tanMarker.scale.setScalar(0.045);
    this.tanMarker.renderOrder = 5;
    this.marsNode.add(this.tanMarker);

    this.sun = new THREE.DirectionalLight(0xffffff, 3.2);
    this.ambient = new THREE.AmbientLight(0xffffff, 0.12);
    this.root.add(this.sun, this.sun.target, this.ambient);

    const loader = new TDSLoader();
    loader.setResourcePath(this.modelTexturePath);
    const cache = {};
    const load = async (file) => {
      cache[file] ??= loader.loadAsync(`${this.assetBase}models/${file}`);
      return (await cache[file]).clone();
    };

    this.craft = {};
    for (const [key, cfg] of Object.entries(MODELS)) {
      const node = new THREE.Group();          // position + attitude in OCC frame
      const body = new THREE.Group();          // km units, scaled by the exaggeration factor
      node.add(body);
      const bus = await load(cfg.bus);
      bus.scale.setScalar(cfg.scale);
      body.add(bus);
      const parts = {};
      for (const [part, file] of Object.entries(cfg.parts)) {
        const pnode = new THREE.Group();
        const mesh = await load(file);
        mesh.scale.setScalar(cfg.scale);
        pnode.add(mesh);
        body.add(pnode);
        parts[part] = pnode;
      }
      node.traverse((o) => {
        if (o.isMesh) {
          o.material = [o.material].flat().map((m) => {
            if (m.map) m.map.colorSpace = THREE.SRGBColorSpace;
            m.side = THREE.DoubleSide;
            if (m.emissive) m.emissive.setHex(0x1c1f26);  // keep dark solar cells readable against space
            return m;
          });
          if (o.material.length === 1) o.material = o.material[0];
        }
      });
      // boresight (MEX +Z, TGO -Y), the antenna axes the forecast's pointing angles refer to
      const dir = key === 'mex' ? new THREE.Vector3(0, 0, 1) : new THREE.Vector3(0, -1, 0);
      const bore = new THREE.ArrowHelper(dir, new THREE.Vector3(), 700, COLORS[key], 120, 60);
      node.add(bore);
      // the true-scale-exaggerated models are still only a few pixels in wide views, so mark them with a glow
      const beacon = new THREE.Sprite(new THREE.SpriteMaterial({ map: glowTexture(COLORS[key]), sizeAttenuation: false,
        depthWrite: false, blending: THREE.AdditiveBlending, opacity: 0.8 }));
      beacon.scale.setScalar(0.03);
      node.add(beacon);
      const label = textSprite(key === 'mex' ? 'Mars Express' : 'TGO', COLORS[key]);
      this.root.add(label);
      this.frame.add(node);
      this.craft[key] = { node, body, parts, bore, label, beacon };
    }

    this.labels = [this.craft.mex.label, this.craft.tgo.label];
    this.setScScale(this.scScale);
    this.setLabelSize(0.02);
  }

  setData(data) {
    this.data = data;
    const [a, b, c] = data.radii;
    this.marsMesh.scale.set(a, c, b);
    const ionoAlt = 200;
    this.iono.scale.set(a + ionoAlt, c + ionoAlt, b + ionoAlt);

    const sun = new THREE.Vector3(...data.sunDir);
    this.sunDirOcc = sun;

    for (const k of ['mex', 'tgo']) {
      const { parts } = this.craft[k];
      for (const [part, node] of Object.entries(parts)) node.position.fromArray(data.offsets[`${k}.${part}`]);
    }

    // dispose old line objects
    for (const o of [this.orbitMex, this.orbitTgo, this.link, this.traceAll, this.traceNow]) {
      if (o) { o.parent.remove(o); o.geometry.dispose(); o.material.dispose(); }
    }
    this.lineMaterials = [];
    const add = (line, parent) => { parent.add(line); this.lineMaterials.push(line.material); return line; };
    this.orbitMex = add(fatLine(data.orbits.mex, COLORS.mex, 2, { opacity: 0.55 }), this.frame);
    this.orbitTgo = add(fatLine(data.orbits.tgo, COLORS.tgo, 2, { opacity: 0.55 }), this.frame);
    this.link = add(fatLine([0, 0, 0, 1, 1, 1], COLORS.link, 3), this.frame);
    this.link.frustumCulled = false;

    // tangent-point profile across the occultation test window, in body-fixed coordinates
    const pts = [];
    this.traceIdx = [];
    const iStart = Math.max(0, Math.floor((data.etStart - data.et0) / data.step));
    const iEnd = Math.min(data.n - 1, Math.ceil((data.etEnd - data.et0) / data.step));
    for (let i = iStart; i <= iEnd; i++) {
      if (data.alt[i] === null || data.alt[i] > 1500) continue;
      pts.push(data.tan[3 * i], data.tan[3 * i + 1], data.tan[3 * i + 2]);
      this.traceIdx.push(i);
    }
    if (pts.length < 6) pts.push(...pts.slice(0, 3), ...pts.slice(0, 3));
    this.traceAll = add(fatLine(pts, COLORS.trace, 2, { opacity: 0.3 }), this.marsNode);
    this.traceNow = add(fatLine(pts, COLORS.trace, 4), this.marsNode);
    if (this._res) this.setResolution(...this._res);
    this.update(data.etOcc - data.et0);
  }

  setResolution(w, h) {
    this._res = [w, h];
    for (const m of this.lineMaterials) m.resolution.set(w, h);
  }

  setScScale(s) {
    this.scScale = s;
    for (const { body } of Object.values(this.craft ?? {})) body.scale.setScalar(s);
  }

  setLabelSize(s) {
    this.labelSize = s;
    const k = 1 / this.root.scale.x;  // screen-sized sprites inherit the root scale, so undo it
    for (const l of this.labels ?? []) l.scale.set(s * l.userData.aspect * k, s * k, 1);
  }

  // Scale the whole km-unit scene (e.g. into AR marker units) while keeping sprites the same size on screen
  setRootScale(s) {
    this.root.scale.setScalar(s);
    this.tanMarker.scale.setScalar(0.045 / s);
    for (const c of Object.values(this.craft)) c.beacon.scale.setScalar(0.03 / s);
    this.setLabelSize(this.labelSize);
  }

  setVisible({ orbits, iono, boresights, labels, trace } = {}) {
    if (orbits !== undefined) { this.orbitMex.visible = this.orbitTgo.visible = orbits; }
    if (iono !== undefined) this.iono.visible = iono;
    if (boresights !== undefined) for (const c of Object.values(this.craft)) c.bore.visible = boresights;
    if (labels !== undefined) for (const c of Object.values(this.craft)) c.label.visible = c.beacon.visible = labels;
    if (trace !== undefined) { this.traceAll.visible = this.traceNow.visible = trace; }
  }

  get duration() { return (this.data.n - 1) * this.data.step; }

  // interpolation helpers over the column arrays; t is seconds since data.et0
  _idx(t) {
    const x = Math.min(Math.max(t / this.data.step, 0), this.data.n - 1);
    const i = Math.min(Math.floor(x), this.data.n - 2);
    return [i, x - i];
  }

  _scalar(key, i, f) {
    const a = this.data[key][i], b = this.data[key][i + 1];
    if (a === null || b === null) return f < 0.5 ? a : b;
    return a + (b - a) * f;
  }

  _vec(key, i, f, out) {
    const d = this.data[key];
    if (d[3 * i] === null || d[3 * i + 3] === null) return null;
    return out.set(d[3 * i] + (d[3 * i + 3] - d[3 * i]) * f,
      d[3 * i + 1] + (d[3 * i + 4] - d[3 * i + 1]) * f,
      d[3 * i + 2] + (d[3 * i + 5] - d[3 * i + 2]) * f);
  }

  _quat(arr, i, f, out) {
    const [qa, qb] = this._q;
    qa.fromArray(arr, 4 * i);
    qb.fromArray(arr, 4 * i + 4);
    return out.slerpQuaternions(qa, qb, f);
  }

  update(t) {
    const d = this.data;
    if (!d) return;
    const [i, f] = this._idx(t);
    const mex = this._vec('mex', i, f, this.craft.mex.node.position);
    const tgo = this._vec('tgo', i, f, this.craft.tgo.node.position);
    this._quat(d.mexQ, i, f, this.craft.mex.node.quaternion);
    this._quat(d.tgoQ, i, f, this.craft.tgo.node.quaternion);
    this._quat(d.marsQ, i, f, this.marsNode.quaternion);
    for (const [k, c] of Object.entries(this.craft)) {
      for (const [part, node] of Object.entries(c.parts)) this._quat(d.parts[k][part], i, f, node.quaternion);
    }

    this.link.geometry.setPositions([mex.x, mex.y, mex.z, tgo.x, tgo.y, tgo.z]);
    const alt = this._scalar('alt', i, f);
    const blocked = alt !== null && alt < 0;
    this.link.material.color.set(blocked ? COLORS.blocked : COLORS.link);
    this.link.material.dashed = blocked;
    this.link.material.opacity = blocked ? 0.7 : 1;
    this.link.computeLineDistances();

    const tan = this._vec('tan', i, f, this.tanMarker.position);
    this.tanMarker.visible = !!tan && alt !== null && alt > -20 && alt < 1500;

    // profile drawn up to "now"
    let k = 0;
    while (k < this.traceIdx.length && this.traceIdx[k] <= i) k++;
    this.traceNow.geometry.instanceCount = Math.max(0, k - 1);

    // sun light in world space (frame is rotated within root)
    this.frame.updateMatrix();
    const sunWorld = this.sunDirOcc.clone().applyMatrix4(new THREE.Matrix4().extractRotation(this.frame.matrix));
    this.sun.position.copy(sunWorld.multiplyScalar(50000));
    this.sun.target.position.set(0, 0, 0);

    for (const c of Object.values(this.craft)) c.label.position.copy(c.node.position).applyMatrix4(this.frame.matrix);

    this.state = {
      t, utc: new Date(d.utc0ms + t * 1000),
      tRel: t - (d.etOcc - d.et0),
      alt, blocked,
      lat: this._scalar('lat', i, f), lon: this._scalar('lon', i, f),
      sza: this._scalar('sza', i, f), lst: this._scalar('lst', i, f),
      range: this._scalar('range', i, f), doppler: this._scalar('doppler', i, f),
      mexOff: this._scalar('mexOff', i, f), tgoOff: this._scalar('tgoOff', i, f),
      inTest: t >= d.etStart - d.et0 && t <= d.etEnd - d.et0,
    };
    return this.state;
  }

  // world-space helpers for camera presets
  worldPos(key) {
    const v = key === 'tan' ? this.tanMarker.getWorldPosition(new THREE.Vector3())
      : this.craft[key].node.getWorldPosition(new THREE.Vector3());
    return v;
  }
}
