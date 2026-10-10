===========
Model guide
===========

.. contents:: Table of Contents
    :depth: 3

This page explains the magnet *models* available for thick elements
(:class:`xtrack.Magnet`, :class:`xtrack.Bend`, :class:`xtrack.RBend`,
:class:`xtrack.Quadrupole`, ...): what each one integrates exactly, what it
leaves to thin kicks, when it is worth using, and where it can go wrong. The
Hamiltonians and the derivation of each map are in the section "Symplectic
integrators and magnet models" of the :doc:`Physics Guide <physicsguide>`.
How to set models, integrators and the number of kicks is shown in
:doc:`lattice_model`.

How a model is built
====================

Every model splits the Hamiltonian of the magnet into a *propagation part*
(solved exactly, or with a closed-form approximation) and a *kick part* (applied
as thin kicks). The kicks are distributed along the element by the
*integrator* (``uniform``, ``teapot``, ``yoshida-4``, ``yoshida-6``,
``yoshida-8``), and the propagation part is applied between them.

.. list-table::
   :header-rows: 1

   * - Model
     - Propagation part
     - Kick part
   * - ``rot-kick-rot-low-order``
     - R
     - K0 K0h K1 K1h Kn
   * - ``rot-kick-rot``
     - R, K0 (Yoshida-4)
     - K0h K1 K1h Kn
   * - ``rot-kick-rot-high-order``
     - R, K0 (Yoshida-6)
     - K0h K1 K1h Kn
   * - ``bend-kick-bend``
     - B
     - K1 K1h Kn
   * - ``mat-kick-mat``
     - M
     - K1h Kn
   * - ``mat-kick-mat-exact``
     - M + exact path correction
     - K1h Kn
   * - ``drift-kick-drift-exact``
     - D
     - h K0 K0h K1 K1h Kn
   * - ``drift-kick-drift-expanded``
     - De
     - h K0 K0h K1 K1h Kn

``R`` is the exact (polar, if ``h`` is non-zero) drift, ``B`` the exact bend,
``M`` the expanded combined-function (bend-quadrupole) map, ``D`` and ``De`` the
exact and expanded drifts. ``K0``, ``K1``, ``K0h``, ``K1h`` and ``Kn`` are the
dipole, quadrupole, weak-focusing, quadrupole-curvature and general multipole
kicks.

Three rules explain most of the behaviour described below:

* Whatever is in the propagation part is integrated exactly (up to the
  approximation of that map), independently of the number of kicks.
* Whatever is in the kick part carries the splitting error of the integrator:
  it decreases with the integrator order and with the number of slices.
* If a model has *no active kick part* for the given element (e.g. a pure
  dipole with ``bend-kick-bend``), the result is exact for any split and extra
  kicks only add cost and round-off. For ``h != 0`` the integrator still runs,
  so use ``uniform`` with ``num_multipole_kicks=1``.

Default models
--------------

With ``adaptive``, the model and integrator depend on the element class:

.. list-table::
   :header-rows: 1

   * - Element
     - Model
     - Integrator
   * - ``Bend``, ``RBend``
     - ``rot-kick-rot``
     - ``yoshida-6``
   * - ``Quadrupole``
     - ``mat-kick-mat``
     - ``uniform``
   * - ``Sextupole``, ``Octupole``, ``Multipole``
     - ``drift-kick-drift-expanded``
     - ``uniform``

Curvature terms in the kick
---------------------------

In a curved frame (``h != 0``) every multipole picks up terms proportional to
``h`` times a higher power of ``x``. Only the first-order terms of the normal
``K0`` and ``K1`` components are implemented (as thin kicks, following eq. 5.15
of the MAD-8 physics manual):

.. math::

    H_{K0h} = \chi h K_0 \frac{x^2}{2}, \qquad
    H_{K1h} = \chi h K_1 \left(\frac{x^3}{3} - \frac{x y^2}{2}\right)

Skew ``K0s``/``K1s``, all components from ``K2`` up (including ``knl``/``ksl``
entries of order 2 and above) and terms of order :math:`h^2` are **not**
corrected, in any model. A sextupole or octupole inside a bend therefore lacks
its :math:`O(h)` feed-down. The relative error is of order :math:`h x`:
negligible for large rings, but up to percent level in compact, strongly bending
lattices. It **does not decrease with the number of kicks or the integrator order**.

Choosing a model
================

``rot-kick-rot`` (default)
--------------------------

*What it is.* The default for ``Bend`` and ``RBend`` (``adaptive`` selects it).
The drift is the exact polar drift, and ``K0`` is nested inside it: each drift
step of the integrator becomes four polar drifts interleaved with three ``K0``
kicks (Yoshida-4). ``K1``, higher multipoles, K0h and K1h are thin kicks.

*Why use it.* It is the general-purpose model for bends: curvature is exact in
the drift, ``K0`` is integrated to higher order than the other strengths, and
all multipole components are supported.

*Issues.*

* Not exact for a pure dipole: K0h is a kick, not part of the nested step (use
  ``bend-kick-bend``). The design orbit is preserved only up to the nested
  truncation error, which grows quickly with the angle per step.
* The nested order is 4 regardless of the integrator, so with the default
  ``yoshida-6`` the ``K0`` part limits the convergence. See
  ``rot-kick-rot-high-order``.
* Each drift step costs four polar drifts (eight for ``rot-kick-rot-high-order``).
* The nesting needs ``|h| >= 1e-8``. For ``h = 0`` (e.g. a steering dipole) the
  model reduces to ``drift-kick-drift-exact`` and ``K0`` is a plain kick. The
  same applies to the other ``rot-kick-rot`` variants.
* Quadrupole and higher components are kicks, with only the curvature terms
  listed above.

``rot-kick-rot-high-order``
---------------------------

Same as ``rot-kick-rot`` with Yoshida-6 (eight polar drifts, seven ``K0``
kicks) for the nested slices. Use it when the dipole and the curvature dominate
the error, as the order of ``yoshida-6`` is otherwise capped by the nesting.

``rot-kick-rot-low-order``
--------------------------

The exact polar drift only; ``K0`` is a kick like all other strengths, so the
dipole carries the splitting error. Cheaper than ``rot-kick-rot``.

``bend-kick-bend``
------------------

*What it is.* The exact bend map ``B`` (``K0`` and ``h``, including the
weak-focusing term) is the propagation part; ``K1``, ``K1h`` and multipoles are
kicks.

*Why use it.* For a pure dipole (no ``K1``, no multipoles) the map is exact for
any number of slices, so ``nslice=1`` is enough and slicing only adds cost and
round-off. It is also the natural choice for pure bending magnets and is a must
for magnets with strong curvature (``h``).

*Issues.* With a ``K1`` component present, the quadrupole is a kick and
converges only as fast as the integrator allows. The map is also prone to numerical
noise when the curvature is small (``h`` close to zero) or the strength is low, so
integrating this model is not recommended.

``mat-kick-mat``
----------------

*What it is.* The expanded combined-function map ``M`` (``K0``, ``K1``, ``h``,
ported from MAD-X) is the propagation part; only ``K1h`` and multipoles are
kicks.

*Why use it.* The quadrupole is part of the exact transport, so a quadrupole or
a combined-function magnet needs no slicing for the linear optics.

*Issues.* The map is expanded: transverse positions advance paraxially
(``px * L / P``, with ``P`` the total momentum) rather than with the exact
``px / pz``. Slicing does not remove this difference. For small amplitudes
and small momentum offset it is negligible, but it sets an error floor that
**does not decrease with the number of slices**.

With ``h != 0``, K0h is part of the linear map ``M`` and only K1h is a kick.

``mat-kick-mat-exact``
----------------------

*What it is.* The same map ``M`` as ``mat-kick-mat``, plus an exact,
momentum-only path correction applied at the boundaries of every kick. The
correction replaces the paraxial transverse advance by the exact one:

.. math::

    \Delta x = p_x\,\ell\left(\frac{1}{p_z} - \frac{1}{1+\delta}\right),
    \qquad
    p_z = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2},

and likewise for ``y``. It is a function of the momenta only, which are not
changed by the correction, so two consecutive corrections of half length are
equal to one of full length. When a kick has nothing between the two
corrections, a single full-length correction is applied.

*Why use it.* When accuracy at large amplitude or large momentum offset
matters (large ``delta``, large ``px``/``py``), and a thick quadrupole or
combined-function magnet is being modelled. The error converges with the
number of slices instead of stalling at the paraxial floor.

*Issues.*

* It is exact only for straight magnets (``h == 0``). When the magnet has
  curvature the model behaves exactly like ``mat-kick-mat``.
* It costs slightly more than ``mat-kick-mat`` (path correction per kick) and
  more per slice than the drift-kick-drift models, see the timings below.
* ``zeta`` is corrected too, so the exact path length is recovered.

``drift-kick-drift-exact``
--------------------------

*What it is.* Exact drift; every strength, including ``K0`` and ``K1``, is a
kick. This is the classical thin-kick scheme and the reference used in the
convergence studies.

*Why use it.* No approximation other than the splitting itself, and the
simplest map. Very good as a reference - reliable map especially for convergence
studies. Best for straight magnets.

*Issues.* In general this map needs a lot of integration steps to have a reasonable accuracy.
The curvature ``h`` is also a kick, not an exact map.

``drift-kick-drift-expanded``
-----------------------------

Same as ``drift-kick-drift-exact`` with the expanded drift, therefore has an
accuracy floor much higher than the exact drift. The expansion error
is not removed by slicing.

Integration and numerical stability
===================================

* Increasing the integrator order is the efficient way to improve accuracy:
  going up one order is roughly equivalent to increasing ``nslice`` by a
  factor 10 at a lower order, at much lower cost. ``yoshida-6`` (the
  ``adaptive`` integrator) is typically the best choice.
* Every drift or map evaluation adds round-off, so fewer, higher-order
  slices are numerically more stable than many low-order slices.
* ``yoshida-4``, ``yoshida-6`` and ``yoshida-8`` use 3, 7 and 15 kicks per
  slice: the number of kicks is rounded up to a complete slice. The older name
  ``yoshida4`` is deprecated and was equivalent to ``yoshida-6``.
* If the model has no kick part for the element, slicing is pointless: a pure
  dipole with ``bend-kick-bend`` is exact in one step.
* If ``num_multipole_kicks`` is zero and kicks are needed, a heuristic is used:
  one kick in the middle for a straight magnet, about 2 kicks/mrad for a bend.
  Check convergence explicitly when the result matters.
* Models that include an expansion (``mat-kick-mat``,
  ``drift-kick-drift-expanded``) have an error that does not depend on the
  number of kicks. Convergence tests should therefore be made against an exact
  reference, not against a higher slice count of the same model.

Which model for which element
=============================

* Straight quadrupole, large amplitude or large ``delta``, converged result
  needed: ``mat-kick-mat-exact``.
* Straight quadrupole, small amplitude, speed: ``mat-kick-mat``.
* Bending magnet that is purely dipolar: ``bend-kick-bend`` with
  ``nslice=1``.
* Bending magnet with quadrupole and multipole content, or strong curvature:
  ``rot-kick-rot``, or ``rot-kick-rot-high-order`` if the
  dipole and curvature dominate the error. Sextupole and higher content in a
  strongly curved magnet lacks its :math:`O(h)` terms in every model.
* Straight dipole (``h = 0``): ``bend-kick-bend``, since ``rot-kick-rot`` treats
  ``K0`` as a plain kick.
* Reference or cross-check calculation: ``drift-kick-drift-exact`` with a high
  order integrator and many slices.
