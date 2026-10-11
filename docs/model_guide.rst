===========
Model guide
===========

.. contents:: Table of Contents
    :depth: 2

This page describes the magnet *models* available for thick elements
(:class:`xtrack.Magnet`, :class:`xtrack.Bend`, :class:`xtrack.RBend`,
:class:`xtrack.Quadrupole`, ...). For each model it explains what is integrated
exactly, what is left to thin kicks, when to use it, and where it can go wrong.

* The Hamiltonians and the derivation of each map are in the section
  "Symplectic integrators and magnet models" of the
  :doc:`Physics Guide <physicsguide>`.
* How to set models, integrators and the number of kicks is shown in
  :doc:`lattice_model`.

Which model for which element
=============================

Follow the chart from the element type to a recommended model.

.. figure:: figures/model_flowchart.svg
   :align: center
   :width: 100%
   :alt: Flow chart for choosing a magnet model

   Choosing a magnet model. Source: ``figures/model_flowchart.dot``.

Two points the chart does not show:

* For a reference or cross-check calculation, use ``drift-kick-drift-exact``
  with a high-order integrator and many slices.
* For any model with a kick part, check convergence in the number of slices
  when the result matters.

.. note::

    Sextupole and higher components inside a strongly curved magnet lack their
    :math:`O(h)` terms in every model (see `Curvature terms in the kick`_).

Models at a glance
==================

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

Propagation parts:

* ``R``: exact drift (polar if ``h`` is non-zero)
* ``B``: exact bend
* ``M``: expanded combined-function (bend-quadrupole) map
* ``D``, ``De``: exact and expanded drifts

Kick parts: ``K0``, ``K1``, ``K0h``, ``K1h`` and ``Kn`` are the dipole,
quadrupole, weak-focusing, quadrupole-curvature and general multipole kicks.

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

How a model is built
====================

Every model splits the Hamiltonian of the magnet into a *propagation part*,
solved exactly or with a closed-form approximation, and a *kick part*, applied
as thin kicks. The *integrator* (``uniform``, ``teapot``, ``yoshida-4``,
``yoshida-6`` or ``yoshida-8``) distributes the kicks along the element, and the
propagation part is applied between them.

Three rules explain most of the behaviour described below:

* The propagation part is integrated exactly, up to the approximation made in
  that map, whatever the number of kicks.
* The kick part carries the splitting error of the integrator, which decreases
  with the integrator order and with the number of slices.
* If a model has *no active kick part* for the element (e.g. a pure dipole with
  ``bend-kick-bend``), the result is exact for any split, and extra kicks only
  add cost and round-off. For ``h != 0`` the integrator still runs, so use
  ``uniform`` with ``num_multipole_kicks=1``.

Curvature terms in the kick
---------------------------

In a curved frame (``h != 0``) every multipole picks up terms proportional to
``h`` times a higher power of ``x``. Only the first-order terms of the normal
``K0`` and ``K1`` components are implemented, as thin kicks (following eq. 5.15
of the MAD-8 physics manual):

.. math::

    H_{K0h} = \chi h K_0 \frac{x^2}{2}, \qquad
    H_{K1h} = \chi h K_1 \left(\frac{x^3}{3} - \frac{x y^2}{2}\right)

The following are **not** corrected, in any model:

* skew ``K0s`` and ``K1s``,
* all components from ``K2`` up, including ``knl``/``ksl`` entries of order 2
  and above,
* terms of order :math:`h^2`.

A sextupole or octupole inside a bend therefore lacks its :math:`O(h)`
feed-down. The relative error is of order :math:`h x`: negligible for large
rings, but up to the percent level in compact, strongly bending lattices. It
**does not decrease with the number of kicks or the integrator order**.

Model details
=============

``rot-kick-rot`` (default for bends)
------------------------------------

*What it is.* The drift is the exact polar drift, with ``K0`` nested inside it:
each drift step of the integrator becomes four polar drifts interleaved with
three ``K0`` kicks (Yoshida-4). ``K1``, higher multipoles, K0h and K1h are thin
kicks.

*Why use it.* It is the general-purpose model for bends. Curvature is exact in
the drift, ``K0`` is integrated to higher order than the other strengths, and
all multipole components are supported.

*Issues.*

* It is not exact for a pure dipole, because K0h is a kick and not part of the
  nested step. For a strong dipole the design orbit can deviate significantly,
  so check against ``bend-kick-bend``.
* The nested order is 4 whatever the integrator, so with the default
  ``yoshida-6`` the ``K0`` part limits the convergence. See
  ``rot-kick-rot-high-order``.
* Each drift step costs four polar drifts, and large bending angles need many
  slices for accurate integration. Prefer ``bend-kick-bend`` for pure dipoles.
* The nesting needs ``|h| >= 1e-8``. For ``h = 0`` (e.g. a steering dipole) the
  model reduces to ``drift-kick-drift-exact`` and ``K0`` is a plain kick. The
  same applies to the other ``rot-kick-rot`` variants.
* Quadrupole and higher components are kicks, with only the curvature terms
  listed above.

``rot-kick-rot-high-order``
---------------------------

*What it is.* ``rot-kick-rot`` with Yoshida-6 for the nested slices: eight polar
drifts interleaved with seven ``K0`` kicks.

*Why use it.* For a strong dipole or a dipole with large curvature, where the
order-4 nesting would otherwise cap the order of ``yoshida-6``.

*Issues.* Each drift step is about twice as expensive as in ``rot-kick-rot``.

``rot-kick-rot-low-order``
--------------------------

*What it is.* The exact polar drift only. ``K0`` is a kick like all other
strengths, so the dipole carries the splitting error.

*Why use it.* It is cheaper than ``rot-kick-rot``.

*Issues.* Lower accuracy than ``rot-kick-rot``, since the dipole is no longer
integrated to higher order.

``bend-kick-bend``
------------------

*What it is.* The exact bend map ``B`` (``K0`` and ``h``, including the
weak-focusing term) is the propagation part. ``K1``, ``K1h`` and multipoles are
kicks.

*Why use it.* For a pure dipole (no ``K1``, no multipoles) the map is exact for
any number of slices, so ``nslice=1`` is enough and slicing only adds cost and
round-off. It is the natural choice for pure bending magnets and a must for
strong curvature.

*Issues.*

* With a ``K1`` component, the quadrupole is a kick and converges only as fast
  as the integrator allows.
* The map is prone to numerical noise when the curvature is small (``h`` close
  to zero) or the strength is low. Integrating it then introduces round-off
  errors and can make optics from finite differences unstable. In that case use
  ``rot-kick-rot``, or get the optics from TPSA.

``mat-kick-mat``
----------------

*What it is.* The expanded combined-function map ``M`` (``K0``, ``K1``, ``h``,
ported from MAD-X) is the propagation part. Only ``K1h`` and multipoles are
kicks. With ``h != 0``, K0h is part of the linear map ``M``.

*Why use it.* The quadrupole is part of the map, so a quadrupole or a
combined-function magnet needs no slicing for the linear optics. The map is very
fast and reasonably accurate for small amplitudes and small momentum offsets.

*Issues.* The map is expanded: transverse positions advance paraxially
(``px * L / P``, with ``P`` the total momentum) rather than with the exact
``px / pz``. The difference is negligible for small amplitudes and momentum
offsets, but it sets an error floor that **does not decrease with the number of
slices**.

``mat-kick-mat-exact``
----------------------

*What it is.* The same map ``M`` as ``mat-kick-mat``, plus an exact,
momentum-only path correction applied at the boundaries of every kick. The
correction replaces the paraxial transverse advance by the exact one:

.. math::

    \Delta x = p_x\,\ell\left(\frac{1}{p_z} - \frac{1}{1+\delta}\right),
    \qquad
    p_z = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2},

and likewise for ``y``. The correction depends only on the momenta, which it
leaves unchanged, so two consecutive corrections of half length equal one of
full length. When a kick has nothing between its two corrections, a single
full-length correction is applied.

*Why use it.* It removes the paraxial approximation, and with it the error floor
of ``mat-kick-mat``, at a modest extra cost. It is the best choice for straight
quadrupoles and straight combined-function magnets, and it converges to the
exact result with the number of slices.

*Issues.*

* It is exact only for straight magnets (``h == 0``). With curvature it behaves
  exactly like ``mat-kick-mat``.
* It costs slightly more than ``mat-kick-mat`` (a path correction per kick), and
  more per slice than the drift-kick-drift models.

``drift-kick-drift-exact``
--------------------------

*What it is.* The exact drift, with every strength, including ``K0`` and
``K1``, applied as a kick. This is the classical thin-kick scheme and the
reference used in the convergence studies.

*Why use it.* It makes no approximation other than the splitting itself, and it
is the simplest map. It is a reliable reference, especially for convergence
studies, and works best for straight magnets.

*Issues.* It needs many integration steps for reasonable accuracy, and the
curvature ``h`` is a kick rather than an exact map.

``drift-kick-drift-expanded``
-----------------------------

*What it is.* ``drift-kick-drift-exact`` with the expanded drift. It is the
default for sextupoles, octupoles and multipoles.

*Why use it.* It is the cheapest drift-kick-drift model, adequate for small
amplitudes.

*Issues.* The accuracy floor is much higher than with the exact drift, and the
expansion error is not removed by slicing.

Integration and numerical stability
===================================

* Raising the integrator order is the efficient way to improve accuracy: one
  order up is roughly equivalent to multiplying ``nslice`` by 10 at the lower
  order, at much lower cost. ``yoshida-6`` (the ``adaptive`` integrator) is
  typically the best choice.
* Every drift or map evaluation adds round-off, so fewer, higher-order slices
  are numerically more stable than many low-order ones.
* ``yoshida-4``, ``yoshida-6`` and ``yoshida-8`` use 3, 7 and 15 kicks per
  slice, so the number of kicks is rounded up to a complete slice. The older name
  ``yoshida4`` is deprecated and was equivalent to ``yoshida-6``.
* If the model has no kick part for the element, slicing is pointless. A pure
  dipole with ``bend-kick-bend`` is exact in one step, and ``mat-kick-mat``
  gains nothing from slicing without multipole content.
* If ``num_multipole_kicks`` is zero and kicks are needed, a heuristic is used:
  one kick in the middle for a straight magnet, about 2 kicks/mrad for a bend.
  Check convergence explicitly when the result matters.
* Models that include an expansion (``mat-kick-mat``,
  ``drift-kick-drift-expanded``) have an error that does not depend on the
  number of kicks. Test convergence against an exact reference, not against a
  higher slice count of the same model.

References
==========

The model names and splittings follow the kick-code tradition of TEAPOT and PTC.
The reference description of the PTC integrators is

* E. Forest, F. Schmidt and E. McIntosh, *Introduction to the Polymorphic
  Tracking Code: Fibre Bundles, Polymorphic Taylor Types and "Exact Tracking"*,
  KEK and CERN, July 2002
  (`ptc_report_2002.pdf <https://madx.web.cern.ch/doc/ptc_report_2002.pdf>`_).

Where the Xsuite models relate to it (section numbers refer to that report;
the correspondence is of concept, not a line-by-line identity of the splitting):

* Yoshida integrators of order 4 and 6 built from the second-order scheme: K.4.2
  (DKD2). In PTC, the integration method is regarded as part of the model.
* Exact bend in cylindrical geometry, with either the polar drift plus multipole
  kick split (cf. ``rot-kick-rot-low-order``) or the exact bend plus multipole
  split (cf. ``bend-kick-bend``): K.4.9 (TEAPOT). This PTC "Teapot split" is
  unrelated to the Xsuite ``teapot`` *integrator*.
* Matrix-kick-matrix integration and the exact-drift correction term
  :math:`C_\text{exact}`, supported by PTC only for straight magnets (cf.
  ``mat-kick-mat-exact``): K.4.6 (KTK) and K.4.7 (TKTF).

The derivation of the maps is in the section "Symplectic integrators and magnet
models" of the :doc:`Physics Guide <physicsguide>`.
