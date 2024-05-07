import os
import sys
import numpy as np


xmin = 1.
rng = np.random.default_rng()


def init_time(dest):
    # iterator and time
    step = np.array(0, dtype=np.uint64)
    time = np.array(0, dtype=np.float64)
    np.save(f"{dest}/step.npy", step)
    np.save(f"{dest}/time.npy", time)
    return


def init_domain(lengths, glsizes, dest):
    # NOTE: cell face has +1 elements
    xf = xmin + np.linspace(0., lengths[0], glsizes[0] + 1, endpoint=True)
    # cell centers are located at the center
    #   of the two neighbouring cell faces,
    #   which are appended by the boundaries
    xc = xmin
    xc = np.append(xc, 0.5 * xf[:-1] + 0.5 * xf[1:])
    xc = np.append(xc, xmin + lengths[0])
    np.save(f"{dest}/xf.npy", np.array(xf, dtype=np.float64))
    np.save(f"{dest}/xc.npy", np.array(xc, dtype=np.float64))
    np.save(f"{dest}/glsizes.npy", np.array(glsizes, dtype=np.uint64))
    np.save(f"{dest}/lengths.npy", np.array(lengths, dtype=np.float64))
    return xf, xc


def init_fluid(lengths, glsizes, xf, xc, dest):
    ux = np.zeros((glsizes[2], glsizes[1], glsizes[0] + 1), dtype=np.float64)
    uy = np.zeros((glsizes[2], glsizes[1], glsizes[0] + 2), dtype=np.float64)
    uz = np.zeros((glsizes[2], glsizes[1], glsizes[0] + 2), dtype=np.float64)
    p  = - 0.5 + np.random.random_sample((glsizes[2], glsizes[1], glsizes[0] + 2))
    ri = xf[ 0]
    ro = xf[-1]
    ui = 1.
    uo = 0.
    a = 1. / (ro**2. - ri**2.) * (+ ro * uo - ri * ui)
    b = 1. / (ro**2. - ri**2.) * (- ri**2. * ro * uo + ro**2. * ri * ui)
    uy1d = a * xc + b / xc
    uy[:, :] = uy1d.reshape(1, glsizes[0] + 2)
    np.save(f"{dest}/ux.npy", ux)
    np.save(f"{dest}/uy.npy", uy)
    np.save(f"{dest}/uz.npy", uz)
    np.save(f"{dest}/p.npy", p)


def init_interface(lengths, glsizes, rf, rc, vfrac, dest):
    shape = (glsizes[2], glsizes[1], glsizes[0] + 2)
    lt = lengths[1]
    lz = lengths[2]
    nt = glsizes[1]
    nz = glsizes[2]
    dt = lt / nt
    dz = lz / nz
    tc = np.linspace(0.5 * dt, lt - 0.5 * dt, nt)
    zc = np.linspace(0.5 * dz, lz - 0.5 * dz, nz)
    zc, tc, rc = np.meshgrid(zc, tc, rc, indexing="ij")
    r0 = 0.5 * rf[0] + 0.5 * rf[-1]
    t0 = 0.5 * lt
    z0 = 0.5 * lz
    x0 = r0 * np.cos(t0)
    y0 = r0 * np.sin(t0)
    xc = rc * np.cos(tc)
    yc = rc * np.sin(tc)
    ri = rf[ 0]
    ro = rf[-1]
    vtot = 0.5 * (ro**2. - ri**2.) * lt * lz
    rad = np.power(vtot * vfrac * 3. / 4. / np.pi, 1. / 3.)
    d = rad - np.sqrt( \
            + np.power(xc - x0, 2.) \
            + np.power(yc - y0, 2.) \
            + np.power(zc - z0, 2.) \
    )
    vof = 1. / (1. + np.exp(- 2. * glsizes[0] * d))
    np.save(f"{dest}/vof.npy", vof)


def main():
    lengths = list()
    lengths.append(float(os.environ["lx"]))
    lengths.append(float(os.environ["ly"]))
    lengths.append(float(os.environ["lz"]))
    glsizes = list()
    glsizes.append(int(os.environ["glisize"]))
    glsizes.append(int(os.environ["gljsize"]))
    glsizes.append(int(os.environ["glksize"]))
    vfrac = float(os.environ["vfrac"])
    dest = sys.argv[1]
    # sanitise
    ndims = len(lengths)
    assert 3 == ndims
    # init and save
    init_time(dest)
    xf, xc = init_domain(lengths, glsizes, dest)
    init_fluid(lengths, glsizes, xf, xc, dest)
    init_interface(lengths, glsizes, xf, xc, vfrac, dest)


main()
