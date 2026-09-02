"""
Interactive viewer and statistics for IASW runs, for use from a python console.

Reads both output layouts automatically:
  SPIC_1D (C++)    raw doubles,          domain/ phi/ charged_particles/
  PIC1D  (Fortran) unformatted records,  domainGrid.dat Phi/ Density/

Typical session
    from iasw_analysis import *
    ls()                    # numbered list of every run found under BASE_DIRS
    stats(0)                # index, name fragment, or full path all work
    compare(0, 1)           # profile difference table
    energy_budget(0, 1)     # kinetic and field energy separately
    profiles(0, 1)          # interactive: slider, or left/right arrow keys
    traces(0, 1)            # conservation traces

Edit BASE_DIRS below to point at your mount, or call set_base(...) at runtime.
Needs only numpy and matplotlib.
"""
import os, glob
import numpy as np
import matplotlib
# Prefer QtAgg so figures open in their own window.  Probe the Qt binding first, because
# matplotlib.use() defers the import and so cannot report an unusable backend.  Falls back
# silently where there is no GUI, and force=True overrides an inline backend that a Jupyter
# or VS Code session has already installed.
try:
    from matplotlib.backends import backend_qtagg as _qt_probe  # noqa: F401
    matplotlib.use("QtAgg", force=True)
except Exception:
    pass

# --------------------------------------------------------------- configuration
# Directories that hold run folders.  Windows paths are fine, use raw strings:
#   BASE_DIRS = [r'Z:\ImplicitPic1D\comment_IASW_2', r'Z:\ImplicitPic1D\comment_IASW']
BASE_DIRS = [
    "//wsl.localhost/Ubuntu/mnt/Supercomputer/home/nsavard/ImplicitPic1D/comment_IASW",
    "//wsl.localhost/Ubuntu/mnt/Supercomputer/home/nsavard/ImplicitPic1D/comment_IASW_2",
    "//wsl.localhost/Ubuntu/mnt/Supercomputer/home/nsavard/ImplicitPic1D/comment_IASW_2/seed_study",
    "//wsl.localhost/Ubuntu/mnt/Supercomputer/home/nsavard/ImplicitPic1D/comment_IASW_2/resolution",
]

_EPS0 = 8.8541878128e-12
_QE = 1.602176634e-19
_cache = {}
_index = []


def set_base(*dirs):
    """Point the viewer at different directories and rescan."""
    global BASE_DIRS
    BASE_DIRS = [str(d) for d in dirs]
    _cache.clear()
    return ls()

# --------------------------------------------------------------------- readers

def _read_cpp(path):
    return np.fromfile(path, dtype=np.float64)


def _read_fortran(path):
    # unformatted sequential: a 4 byte length marker at each end of the record
    raw = np.fromfile(path, dtype=np.uint8)
    return np.frombuffer(raw[4:-4].tobytes(), dtype=np.float64)


def _table(path, skip=1):
    """Whitespace separated numeric table, skipping `skip` header lines."""
    if not os.path.isfile(path):
        return None
    rows = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i < skip or not line.strip():
                continue
            try:
                rows.append([float(v) for v in line.split()])
            except ValueError:
                continue
    if not rows:
        return None
    width = min(len(r) for r in rows)
    return np.array([r[:width] for r in rows])


def _named_columns(path):
    """Read a one-row table whose header names the columns, keyed by name."""
    if not os.path.isfile(path):
        return {}
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    if len(lines) < 2:
        return {}
    names = [n.strip() for n in lines[0].split(',')]
    vals = [float(v) for v in lines[1].split()]
    return {n: vals[i] for i, n in enumerate(names) if i < len(vals)}


def _column(cols, *keys):
    """First column whose name contains any of the given fragments."""
    for k in keys:
        for name, v in cols.items():
            if k.lower() in name.lower():
                return v
    return None


def _properties(path):
    """Property table whose first column is a name rather than a number."""
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path).read().splitlines()[1:]:
        parts = line.split()
        if len(parts) >= 4:
            try:
                out[parts[0]] = [float(v) for v in parts[1:]]
            except ValueError:
                pass
    return out


class Run:
    """One simulation directory, either code's layout."""

    def __init__(self, path):
        self.path = os.path.abspath(str(path).rstrip('/\\'))
        self.name = os.path.basename(self.path)
        if os.path.isdir(os.path.join(self.path, 'domain')):
            self.code = 'SPIC_1D'
            self._load_cpp()
        elif os.path.isfile(os.path.join(self.path, 'domainGrid.dat')):
            self.code = 'PIC1D'
            self._load_fortran()
        else:
            raise IOError('%s does not look like either output layout' % self.path)

    # -- SPIC_1D ---------------------------------------------------------
    def _load_cpp(self):
        self.read = _read_cpp
        par = _table(os.path.join(self.path, 'domain', 'parameters.dat'))[0]
        self.number_nodes, self.number_cells = int(par[1]), int(par[2])
        self.left_bc, self.right_bc, self.length = int(par[3]), int(par[4]), par[5]
        self.grid = self.read(os.path.join(self.path, 'domain', 'grid.dat'))
        ic = _table(os.path.join(self.path, 'initial_condition.dat'))[0]
        self.num_mpi, self.num_threads, self.scheme_id = int(ic[0]), int(ic[1]), int(ic[2])
        self.del_t, self.num_diag = ic[4], int(ic[5])
        self.scheme = {0: 'MC-PIC', 1: 'EC-PIC', 2: 'I-NGP', 3: 'I-CIC'}.get(self.scheme_id, '?')
        self.layout = 'centred' if self.scheme_id == 3 else 'nodal'
        self.dx = np.diff(self.grid)
        self.species = sorted(os.listdir(os.path.join(self.path, 'charged_particles')))
        self._phi = lambda i: os.path.join(self.path, 'phi', 'potential_%d.dat' % i)
        self._den = lambda s, i: os.path.join(self.path, 'charged_particles', s, 'density', 'density_%d.dat' % i)
        self._tmp = lambda s, i: os.path.join(self.path, 'charged_particles', s, 'temperature', 'cell_temp_%d.dat' % i)
        g = _table(os.path.join(self.path, 'global_diagnostic_data.dat'))
        self.trace_diag = np.arange(len(g))
        self.time, self.steps = g[:, 0], g[:, 1]
        self.momentum, self.energy = g[:, 2], g[:, 5]
        fd = _table(os.path.join(self.path, 'field_diagnostics.dat'))
        self.field_energy, self.gauss = (fd[:, 0], fd[:, 1]) if fd is not None else (None, None)
        t = _table(os.path.join(self.path, 'simulation_timing_data.dat'))
        self.wall, self.field_time, self.particle_time = (t[-1, 0], t[-1, 1], t[-1, 2]) if t is not None else (None,) * 3
        self.n_particles, self._ke = {}, {}
        for s in self.species:
            n = _table(os.path.join(self.path, 'charged_particles', s, 'number_diagnostics.dat'))
            self.n_particles[s] = n[:, 0] if n is not None else None
            d = _table(os.path.join(self.path, 'charged_particles', s, 'energy_diagnostics.dat'))
            p = _properties(os.path.join(self.path, 'charged_particles', s, 'particle_properties.dat')).get(s)
            self._ke[s] = 0.5 * p[0] * p[2] * d[:, 3] if (d is not None and p) else None
        # the t=0 row double counts: the velocity sums are set at load and then
        # accumulated again by the first diagnostic, so skip it in any energy work
        self.first_valid = 1

    # -- PIC1D -----------------------------------------------------------
    def _load_fortran(self):
        self.read = _read_fortran
        self.grid = self.read(os.path.join(self.path, 'domainGrid.dat'))
        bc = np.fromfile(os.path.join(self.path, 'domainBoundaryConditions.dat'), dtype=np.int32)
        n_bc = len(bc) - 2 if len(bc) > 3 else len(bc)     # strip the record markers
        self.left_bc, self.right_bc = (int(bc[1]), int(bc[-2])) if len(bc) > 3 else (3, 3)
        # The CIC scheme stores the potential nodes, which are the cell centres, so its grid
        # is one shorter than its boundary condition array.  Every other scheme stores the
        # cell edges, and the two arrays are the same length.
        self.layout = 'centred' if n_bc == len(self.grid) + 1 else 'nodal'
        dxfile = os.path.join(self.path, 'domainDxDl.dat')
        self.dx = self.read(dxfile) if os.path.isfile(dxfile) else np.diff(self.grid)
        if self.layout == 'centred':
            self.number_cells = len(self.grid)
            self.number_nodes = self.number_cells + 1
            self.length = self.grid[-1] - self.grid[0] + 0.5 * (self.dx[0] + self.dx[-1])
            self.scheme = 'Fortran CIC'
        else:
            self.number_nodes = len(self.grid)
            self.number_cells = self.number_nodes - 1
            self.length = self.grid[-1] - self.grid[0]
            self.scheme = 'Fortran nodal'
        cols = _named_columns(os.path.join(self.path, 'InitialConditions.dat'))
        self.del_t = _column(cols, 'Delta t')
        self.num_diag = int(_column(cols, 'numDiag')) + 1
        self.num_threads, self.num_mpi = int(_column(cols, 'numThread')), 1
        self.scheme_id = -1
        self.species = sorted({os.path.basename(p).split('_')[1]
                               for p in glob.glob(os.path.join(self.path, 'Density', 'density_*_*.dat'))
                               if 'Average' not in p})
        self._phi = lambda i: os.path.join(self.path, 'Phi', 'phi_%d.dat' % i)
        self._den = lambda s, i: os.path.join(self.path, 'Density', 'density_%s_%d.dat' % (s, i))
        self._tmp = lambda s, i: os.path.join(self.path, 'Temperature', 'Temp_%s_%d.dat' % (s, i))
        gpath = os.path.join(self.path, 'GlobalDiagnosticData.dat')
        g = _table(gpath)
        # this table has no t=0 row, so its first entry is dump 1
        self.trace_diag = np.arange(1, len(g) + 1)
        # The explicit and CIC variants write different columns in different orders, so the
        # header is read rather than assuming positions.
        names = [n.strip().lower() for n in open(gpath).read().splitlines()[0].split(',')]
        def col(*keys, **kw):
            for k in keys:
                for i, n in enumerate(names):
                    if k in n and i < g.shape[1]:
                        return g[:, i]
            return kw.get('default')
        self.time = col('time')
        self.momentum = col('momentum')
        self.energy = col('energy total', 'totalenergy', 'energytotal')
        if self.energy is None:
            self.energy = g[:, 5]
        self.gauss = col('gausserror', 'gauss error')
        self.charge_error = col('chargeerror', 'charge error')
        self.energy_error = col('energyerror', 'energy error')
        self.field_energy = None
        t = _table(os.path.join(self.path, 'SimulationFinalData.dat'))
        if t is not None:
            self.wall, self.field_time, self.particle_time = t[0, 0], t[0, 1], t[0, 2]
            self.steps = np.linspace(0, t[0, 4], len(self.time))
        else:
            self.wall = self.field_time = self.particle_time = None
            self.steps = self.trace_diag.astype(float)
        props = _properties(os.path.join(self.path, 'ParticleProperties.dat'))
        self.n_particles, self._ke = {}, {}
        for s in self.species:
            d = _table(os.path.join(self.path, 'ParticleDiagnostic_%s.dat' % s))
            self.n_particles[s] = d[:, 5] if d is not None else None
            if d is not None and s in props:
                self._ke[s] = 1.5 * d[:, 6] * _QE * d[:, 5] * props[s][2]
            else:
                self._ke[s] = None
        self.first_valid = 0

    # -- grids -----------------------------------------------------------
    @property
    def field_x(self):
        """Abscissa for potential and density."""
        if self.layout == 'centred':
            # SPIC_1D stores the cell edges, PIC1D stores the centres directly
            return 0.5 * (self.grid[:-1] + self.grid[1:]) if len(self.grid) == self.number_cells + 1 else self.grid
        return self.grid

    @property
    def cell_x(self):
        if self.layout == 'centred' and len(self.grid) == self.number_cells:
            return self.grid
        return 0.5 * (self.grid[:-1] + self.grid[1:])

    # -- data ------------------------------------------------------------
    def n_dumps(self):
        i = 0
        while os.path.isfile(self._phi(i)):
            i += 1
        return i

    def phi(self, i):
        return self.read(self._phi(i))

    def density(self, s, i):
        return self.read(self._den(s, i))

    def temperature(self, s, i):
        p = self._tmp(s, i)
        return self.read(p) if os.path.isfile(p) else None

    def potential_energy(self, i):
        """Field energy from the phi dump, computed the same way for either code.

        For a centred layout the potential sits at the cell centres, so the field lives in
        the gaps between them, with one extra gap across the periodic seam.
        """
        phi, x = self.phi(i), self.field_x
        gaps = np.diff(phi) ** 2 / np.diff(x)
        total = np.sum(gaps)
        if self.layout == 'centred':
            wrap = self.length - (x[-1] - x[0])
            if self.left_bc == 3:
                total += (phi[-1] - phi[0]) ** 2 / wrap
            else:
                half_l, half_r = x[0], self.length - x[-1]
                total += phi[0] ** 2 / half_l + phi[-1] ** 2 / half_r
        return 0.5 * _EPS0 * total

    def kinetic_energy(self):
        """(dump indices, total kinetic energy) taken from the particle diagnostics."""
        total = None
        for s in self.species:
            k = self._ke.get(s)
            if k is None:
                return None, None
            total = k if total is None else total + k
        return self.trace_diag, total

    def __repr__(self):
        return '<Run %s [%s %s]>' % (self.name, self.code, self.scheme)

# ------------------------------------------------------------------ discovery

def ls(quiet=False):
    """List every run found under BASE_DIRS.  Returns the list of paths."""
    global _index
    _index = []
    for base in BASE_DIRS:
        if not os.path.isdir(base):
            continue
        for entry in sorted(os.listdir(base)):
            p = os.path.join(base, entry)
            if os.path.isdir(p) and (os.path.isdir(os.path.join(p, 'domain'))
                                     or os.path.isfile(os.path.join(p, 'domainGrid.dat'))):
                _index.append(p)
    if not quiet:
        if not _index:
            print('no runs found under:')
            for b in BASE_DIRS:
                print('   ', b, '' if os.path.isdir(b) else '  (missing)')
            print('edit BASE_DIRS at the top of the file, or call set_base(...)')
        else:
            print('%-4s %-9s %-15s %-7s %-9s %s' % ('idx', 'code', 'scheme', 'cells', 'dumps', 'name'))
            for i, p in enumerate(_index):
                try:
                    r = get(i)
                    print('%-4d %-9s %-15s %-7d %-9d %s' % (i, r.code, r.scheme, r.number_cells, r.n_dumps(), r.name))
                except Exception as exc:
                    print('%-4d %-9s %s' % (i, 'unreadable', os.path.basename(p)), exc)
    return _index


def get(which):
    """Resolve an index, a fragment of a name, or a path into a Run (cached)."""
    if isinstance(which, Run):
        return which
    if not _index:
        ls(quiet=True)
    if isinstance(which, (int, np.integer)):
        path = _index[which]
    elif os.path.isdir(str(which)):
        path = str(which)
    else:
        hits = [p for p in _index if str(which).lower() in os.path.basename(p).lower()]
        if not hits:
            raise KeyError('no run matching %r; try ls()' % which)
        if len(hits) > 1:
            print('several runs match %r, using the first:' % which)
            for h in hits:
                print('   ', os.path.basename(h))
        path = hits[0]
    if path not in _cache:
        _cache[path] = Run(path)
    return _cache[path]

# --------------------------------------------------------------------- reports

def stats(*which):
    """Print a summary for each run."""
    for w in (which or [0]):
        r = get(w)
        nd = r.n_dumps()
        bar = '-' * 74
        print(bar)
        print(r.name)
        print('  ' + r.path)
        print(bar)
        print('  code                 %s   scheme %s' % (r.code, r.scheme))
        print('  grid                 %d cells, %d nodes, L = %.6f m'
              % (r.number_cells, r.number_nodes, r.length))
        print('  boundaries           left %d, right %d' % (r.left_bc, r.right_bc))
        print('  time step            %.6e s' % r.del_t)
        print('  parallel             %d MPI x %d threads' % (r.num_mpi, r.num_threads))
        print('  diagnostics written  %d of %d' % (nd, r.num_diag))
        if r.wall is not None:
            extra = ('   (field %.1f s, particles %.1f s)' % (r.field_time, r.particle_time)
                     if r.field_time is not None else '')
            print('  wall time            %.1f s%s' % (r.wall, extra))
        if len(r.time) > 1:
            print('  simulated to         %.6e s over %d steps' % (r.time[-1], int(r.steps[-1])))
        v = r.first_valid
        E, P = r.energy[v:], r.momentum[v:]
        if len(E) > 1:
            print('  total energy column  %.9e -> %.9e J/m^2' % (E[0], E[-1]))
            print('    max rel drift      %.3e' % (np.abs(E - E[0]).max() / abs(E[0])))
        if len(P) > 1 and P[0] != 0:
            print('  total momentum       %.6e -> %.6e kg/m/s' % (P[0], P[-1]))
            print('    max rel drift      %.3e' % (np.abs(P - P[0]).max() / abs(P[0])))
        if r.gauss is not None and len(r.gauss) > 1:
            print('  Gauss residual       mean %.3e, max %.3e'
                  % (np.mean(r.gauss[1:]), np.max(r.gauss[1:])))
        for s in r.species:
            n = r.n_particles.get(s)
            head = '  species %-12s' % s
            if n is not None and len(n):
                head += ' N = %.4g (%.0f per cell)' % (n[0], n[0] / r.number_cells)
                if n[-1] != n[0]:
                    head += ' -> %.4g' % n[-1]
            print(head)
            if nd:
                n0, nl = r.density(s, 0), r.density(s, nd - 1)
                print('    density t=0        min %.4e  max %.4e  mean %.4e' % (n0.min(), n0.max(), n0.mean()))
                print('    density final      min %.4e  max %.4e  mean %.4e' % (nl.min(), nl.max(), nl.mean()))
        if nd:
            p0, pl = r.phi(0), r.phi(nd - 1)
            print('  potential t=0        min %+.4e  max %+.4e V' % (p0.min(), p0.max()))
            print('  potential final      min %+.4e  max %+.4e V' % (pl.min(), pl.max()))
        print(bar)


def energy_budget(*which):
    """Kinetic and field energy separately.

    The stored total energy column is not comparable between the two codes: in a
    leapfrog the velocity used for kinetic energy sits half a step away from the
    field, so what that column means depends on where each code samples it.  The
    kinetic and field parts, taken from the particle diagnostics and from the phi
    dumps, are directly comparable.
    """
    for w in (which or [0]):
        r = get(w)
        idx, ke = r.kinetic_energy()
        if ke is None:
            print('%s: no particle energy diagnostic' % r.name)
            continue
        nd = r.n_dumps()
        keep = [(j, i) for j, i in enumerate(idx) if r.first_valid <= i < nd]
        if len(keep) < 2:
            continue
        pe = np.array([r.potential_energy(i) for _, i in keep])
        k = np.array([ke[j] for j, _ in keep])
        tot = k + pe
        print('-' * 74)
        print('%s  [%s %s]' % (r.name, r.code, r.scheme))
        print('%5s %18s %14s %18s' % ('diag', 'kinetic (J/m^2)', 'field (J/m^2)', 'sum'))
        for n, (_, i) in enumerate(keep):
            if n % max(1, len(keep) // 6) == 0 or n == len(keep) - 1:
                print('%5d %18.9e %14.4e %18.9e' % (i, k[n], pe[n], tot[n]))
        print('  d(kinetic) = %+.4e   d(field) = %+.4e   d(sum) = %+.4e'
              % (k[-1] - k[0], pe[-1] - pe[0], tot[-1] - tot[0]))
        print('  sum conserved to %.3e relative' % (np.abs(tot - tot[0]).max() / abs(tot[0])))
        print('-' * 74)


def periodic_interp(x_src, y_src, x_target, length):
    xx = np.concatenate([x_src - length, x_src, x_src + length])
    yy = np.concatenate([y_src, y_src, y_src])
    return np.interp(x_target, xx, yy)


def compare(*which, **kw):
    """Profile difference of every run against the first."""
    species = kw.get('species', 'ion')
    runs = [get(w) for w in which]
    if len(runs) < 2:
        print('compare needs at least two runs, e.g. compare(0, 1)')
        return
    a = runs[0]
    sp = species if species in a.species else a.species[0]
    nd = min(r.n_dumps() for r in runs)
    print('\n%s difference against %s  (L2 as %% of profile amplitude)' % (sp, a.name))
    print('%5s  ' % 'diag' + '  '.join('%26s' % r.name[:26] for r in runs[1:]))
    print('%5s  ' % '' + '  '.join('%26s' % 'density / potential' for _ in runs[1:]))
    for i in range(nd):
        na, pa = a.density(sp, i), a.phi(i)
        amp_n, amp_p = na.max() - na.min(), pa.max() - pa.min()
        cells = []
        for r in runs[1:]:
            nb = periodic_interp(r.field_x, r.density(sp, i), a.field_x, a.length)
            pb = periodic_interp(r.field_x, r.phi(i), a.field_x, a.length)
            dn = np.sqrt(np.mean((na - nb) ** 2)) / amp_n * 100
            if amp_p < 1e-9:      # neutral start, phi is at rounding level
                cells.append('%11.3f%% %12s' % (dn, 'n/a'))
            else:
                dp = np.sqrt(np.mean((pa - pb) ** 2)) / amp_p * 100
                cells.append('%11.3f%% %11.3f%%' % (dn, dp))
        print('%5d  ' % i + '  '.join('%26s' % c for c in cells))
    print()

# ----------------------------------------------------------------------- plots

def _plt():
    import matplotlib.pyplot as plt
    try:
        plt.ion()
    except Exception:
        pass
    return plt


def profiles(*which, **kw):
    """Interactive profile window.  Slider, or left/right arrow keys."""
    plt = _plt()
    from matplotlib.widgets import Slider
    species = kw.get('species', 'ion')
    runs = [get(w) for w in (which or [0])]
    sp = species if species in runs[0].species else runs[0].species[0]
    nd = min(r.n_dumps() for r in runs)
    if nd == 0:
        print('no diagnostics written yet')
        return
    colors = ['#111111', '#c0392b', '#2471a3', '#1e8449']
    fig, ax = plt.subplots(2, 2, figsize=(12, 7.5))
    fig.canvas.manager.set_window_title('IASW profiles')
    fig.subplots_adjust(bottom=0.14, hspace=0.3, wspace=0.22)
    ax = ax.ravel()
    lines = [[] for _ in ax]
    for k, r in enumerate(runs):
        st = dict(color=colors[k % len(colors)], lw=1.0, ls='-' if k == 0 else '--',
                  label='%s [%s]' % (r.name[:32], r.scheme))
        e = 'e' if 'e' in r.species else r.species[0]
        lines[0].append(ax[0].plot(r.field_x, r.density(sp, 0), **st)[0])
        lines[1].append(ax[1].plot(r.field_x, r.density(e, 0), **st)[0])
        lines[2].append(ax[2].plot(r.field_x, r.phi(0), **st)[0])
        T = r.temperature(sp, 0)
        lines[3].append(ax[3].plot(r.cell_x, T, **st)[0] if T is not None else None)
    for a, t in zip(ax, ['%s density (m^-3)' % sp, 'electron density (m^-3)',
                         'potential (V)', '%s temperature (eV)' % sp]):
        a.set_title(t, fontsize=10)
        a.set_xlabel('x (m)')
        a.grid(alpha=0.3, lw=0.5)
    ax[0].legend(fontsize=8)
    slider = Slider(fig.add_axes([0.13, 0.04, 0.75, 0.03]), 'diagnostic',
                    0, max(nd - 1, 1), valinit=0, valstep=1)

    def update(_):
        i = int(slider.val)
        for k, r in enumerate(runs):
            e = 'e' if 'e' in r.species else r.species[0]
            lines[0][k].set_ydata(r.density(sp, i))
            lines[1][k].set_ydata(r.density(e, i))
            lines[2][k].set_ydata(r.phi(i))
            if lines[3][k] is not None:
                lines[3][k].set_ydata(r.temperature(sp, i))
        for a in ax:
            a.relim(); a.autoscale_view(scalex=False)
        r0 = runs[0]
        j = np.searchsorted(r0.trace_diag, i)
        t = r0.time[min(j, len(r0.time) - 1)]
        fig.suptitle('diagnostic %d of %d     t = %.4e s' % (i, nd - 1, t), fontsize=11)
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key in ('right', 'left'):
            slider.set_val(min(nd - 1, max(0, int(slider.val) + (1 if event.key == 'right' else -1))))

    slider.on_changed(update)
    fig.canvas.mpl_connect('key_press_event', on_key)
    fig._slider = slider          # keep the widget alive
    update(0)
    return fig


def traces(*which):
    """Conservation traces.  Kinetic and field energy, plus the stored columns."""
    plt = _plt()
    runs = [get(w) for w in (which or [0])]
    colors = ['#111111', '#c0392b', '#2471a3', '#1e8449']
    fig, ax = plt.subplots(2, 2, figsize=(12, 7))
    fig.canvas.manager.set_window_title('IASW conservation')
    ax = ax.ravel()
    for k, r in enumerate(runs):
        st = dict(color=colors[k % len(colors)], lw=1.2, ls='-' if k == 0 else '--',
                  marker='o', ms=3, label='%s [%s]' % (r.name[:32], r.scheme))
        v = r.first_valid
        sel = r.trace_diag >= v
        t, E, P = r.time[sel], r.energy[sel], r.momentum[sel]
        idx, ke = r.kinetic_energy()
        nd = r.n_dumps()
        if ke is not None:
            keep = [(j, i) for j, i in enumerate(idx) if v <= i < nd]
            kk = np.array([ke[j] for j, _ in keep])
            pe = np.array([r.potential_energy(i) for _, i in keep])
            tt = np.array([r.time[np.searchsorted(r.trace_diag, i)] for _, i in keep])
            ax[0].plot(tt, (kk + pe - (kk[0] + pe[0])) / abs(kk[0] + pe[0]), **st)
            ax[2].plot(tt, pe, **st)
        ax[1].plot(t, (E - E[0]) / abs(E[0]), **st)
        if P[0] != 0:
            ax[3].plot(t, (P - P[0]) / abs(P[0]), **st)
    for a, t in zip(ax, ['relative drift of kinetic + field energy',
                         'relative drift of the stored total energy column',
                         'field energy (J/m^2)', 'relative momentum drift']):
        a.set_title(t, fontsize=10)
        a.set_xlabel('time (s)')
        a.grid(alpha=0.3, lw=0.5)
    ax[0].legend(fontsize=8)
    fig.tight_layout()
    return fig


def help_():
    print(__doc__)


if __name__ == '__main__':
    print(__doc__)
    ls()
