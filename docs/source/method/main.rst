######
Method
######

.. include:: /reference/reference.txt

Using the second-order accurate finite-difference scheme or the first-order upwind scheme is unfavourable due to its unstable and too diffusive nature, respectively.
To integrate the equation stably while preserving interfacial structures which are inherently sharp, we need a procedure called surface reconstruction, which is accomplished by means of the THINC/QQ scheme (|XIE2017|) in cylindrical coordinates.

**********************
Surface reconstruction
**********************

We consider :math:`\hat{H}`: the diffused representation of the indicator function :math:`H` for each cell.
For each cell, we consider a piecewise planar function

.. math::

    P
    \equiv
    \vec{n}
    \cdot
    \vec{x}
    +
    d.

Note that surface reconstructions are performed on the computational coordinate system rather than the original cylindrical coordinate system.
Although this does not give planar representations in the physical coordinate systems, they give similar results to each other since both coordinates are orthogonal, and this greatly simplifies the whole implementation.

To start we focus on finding :math:`\vec{n}`: surface normal.
By definition the direction of :math:`\vec{n}` is given by the gradient of the indicator function:

.. math::

    \pder{H}{x_i}
    =
    \egr
    \frac{1}{\hvr}
    \pder{H}{\gvr}
    +
    \egt
    \frac{1}{\hvt}
    \pder{H}{\gvt}
    +
    \egz
    \frac{1}{\hvz}
    \pder{H}{\gvz}.

Numerically this is approximated as the gradient of :math:`\phi` and, following Youngs' approach, this is first computed at each cell vertex:

.. math::

    \vat{\dder{\phi}{x_i}}{\cpidx{i},\cpidx{j},\cpidx{k}}
    =
    \egr
    \frac{1}{\hvr}
    \dif{\phi}{\gvr}
    +
    \egt
    \frac{1}{\hvt}
    \dif{\phi}{\gvt}
    +
    \egz
    \frac{1}{\hvz}
    \dif{\phi}{\gvz},

which is implemented as follows:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
    :language: c
    :tag: compute surface gradient at each cell vertex

The normal vector :math:`\vec{n}` is obtained by normalising this, i.e.,

.. math::

    \vat{\vec{n}}{\cpidx{i},\cpidx{j},\cpidx{k}}
    =
    \egr
    n_{\gvr}
    +
    \egt
    n_{\gvt}
    +
    \egz
    n_{\gvz},

which is implemented as follows:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
    :language: c
    :tag: compute surface normal at each cell vertex

In addition to the normal vector at cell vertices which are adopted to compute local curvature, we need to compute surface normal at cell centers to find :math:`P`, which are obtained by averaging the normal vector of the surrounding cell vertices, which are implemented as follows:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
    :language: c
    :tag: compute surface normal at cell center

The position of the surface is found by approximating the definition of the volume-of-fluid :math:`\phi`:

.. math::

    \int_{V^\xi} \hat{H} \left( \xi^r, \xi^\theta, \xi^z \right) dV^\xi
    =
    \vat{J}{\ccidx{i},\ccidx{j},\ccidx{k}}
    \vat{\phi}{\ccidx{i},\ccidx{j},\ccidx{k}}

as

.. math::

    \sum_k
    \sum_j
    \sum_i
    w_i
    w_j
    w_k
    \hat{H} \left( \xi_i^r, \xi_j^\theta, \xi_k^z \right)
    =
    \vat{J}{\ccidx{i},\ccidx{j},\ccidx{k}}
    \vat{\phi}{\ccidx{i},\ccidx{j},\ccidx{k}}.

Utilising the so-called mid-point rule (first-order Gauss-Legendre quadrature) for the given computational coordinate system yields

.. math::

    \hat{H} \left( 0, 0, 0 \right)
    =
    \frac{
        1
    }{
        1
        +
        \exp
        \left(
            - 2 \beta d
        \right)
    }
    =
    \vat{\phi}{\ccidx{i},\ccidx{j},\ccidx{k}},

or

.. math::

    d
    =
    -
    \frac{1}{2 \beta}
    \log \left(
        \frac{1}{\vat{\phi}{\ccidx{i},\ccidx{j},\ccidx{k}}}
        -
        1
    \right),

which is implemented here:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
    :language: c
    :tag: compute surface position

*********
Advection
*********

Through the surface reconstruction, we obtain the piecewise linear function :math:`P` and its diffused representation :math:`\hat{H}` for each cell and are ready to evaluate the fluxes integrated on cell faces:

.. math::

  \int_{\partial V^\xi} u_j \hat{H} dS_{\xi^j}
  =
  u_j \int_{\partial V^\xi} \hat{H} dS_{\xi^j}

to update :math:`\phi`.

Note that the velocities on the corresponding cell surfaces are taken out of the integrals as they are assumed to be constant on each cell face in the framework of finite-difference methods.
To evaluate the surface integrals, we adopt the two-point Gaussian quadratures.
Note that, depending on the sign of the local velocity, :math:`\psi` computed in the upwind cell should be used due to the stability requirement.

Finally, by normalising the integrated values by the corresponding cell surface areas,  we obtain the fluxes as follows.

Radial fluxes:

.. math::

    \vat{\mathcal{F}_\vr}{
        i \pm \frac{1}{2},
        j,
        k
    }
    \equiv
    \vat{\ur}{i \pm \frac{1}{2}, j, k}
    \sum_{n = 0}^1
    \sum_{m = 0}^1
    w_{\gvt_m}
    w_{\gvz_n}
    \hat{H}
    \left(
        \gvr,
        \gvt_m,
        \gvz_n
    \right)

.. myliteralinclude:: /../../src/interface/update/fluxx.c
   :language: c
   :tag: evaluate flux

Azimuthal fluxes:

.. math::

    \vat{\mathcal{F}_\vt}{
        i,
        j \pm \frac{1}{2},
        k
    }
    \equiv
    \vat{\ut}{i, j \pm \frac{1}{2}, k}
    \sum_{n = 0}^1
    \sum_{l = 0}^1
    w_{\gvr_l}
    w_{\gvz_n}
    \hat{H}
    \left(
        \gvr_l,
        \gvt,
        \gvz_n
    \right).

.. myliteralinclude:: /../../src/interface/update/fluxy.c
   :language: c
   :tag: evaluate flux

Axial fluxes:

.. math::

    \vat{\mathcal{F}_\vz}{
        i,
        j,
        k \pm \frac{1}{2}
    }
    \equiv
    \vat{\uz}{i, j, k \pm \frac{1}{2}}
    \sum_{m = 0}^1
    \sum_{l = 0}^1
    w_{\gvr_l}
    w_{\gvt_m}
    \hat{H}
    \left(
        \gvr_l,
        \gvt_m,
        \gvz
    \right).

.. myliteralinclude:: /../../src/interface/update/fluxz.c
   :language: c
   :tag: evaluate flux

These values are used to update :math:`\phi` defined at each cell center :math:`\left( \ccidx{i},\ccidx{j},\ccidx{k} \right)` following

.. math::

    \pder{\phi}{t}
    =
    -
    \frac{1}{J}
    \dif{}{\gvr}
    \left(
        \jhvr
        \mathcal{F}_\vr
    \right)
    -
    \frac{1}{J}
    \dif{}{\gvt}
    \left(
        \jhvt
        \mathcal{F}_\vt
    \right)
    -
    \frac{1}{J}
    \dif{}{\gvz}
    \left(
        \jhvz
        \mathcal{F}_\vz
    \right),

which is implemented as follows:

.. myliteralinclude:: /../../src/interface/update/main.c
   :language: c
   :tag: compute right-hand-side of advection equation

