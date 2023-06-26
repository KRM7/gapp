Multi-objective algorithms
===================================================

.. doxygennamespace:: gapp::algorithm
   :project: gapp
   :desc-only:


class NSGA2
---------------------------------------------------

.. code-block::

   #include <algorithm/nsga2.hpp>

.. doxygenclass:: gapp::algorithm::NSGA2
   :project: gapp
   :members:
   :protected-members:


class NSGA3
---------------------------------------------------

.. code-block::

   #include <algorithm/nsga3.hpp>

.. doxygenclass:: gapp::algorithm::NSGA3
   :project: gapp
   :members:
   :protected-members:


Reference line generation methods
---------------------------------------------------

.. doxygennamespace:: gapp::algorithm::reflines
   :project: gapp
   :desc-only:

.. code-block::

   #include <algorithm/reference_lines.hpp>

.. doxygenfunction:: gapp::algorithm::reflines::quasirandomSimplexPointsMirror
   :project: gapp

.. doxygenfunction:: gapp::algorithm::reflines::quasirandomSimplexPointsSort
   :project: gapp

.. doxygenfunction:: gapp::algorithm::reflines::quasirandomSimplexPointsRoot
   :project: gapp

.. doxygenfunction:: gapp::algorithm::reflines::quasirandomSimplexPointsLog
   :project: gapp

.. doxygenfunction:: gapp::algorithm::reflines::pickSparseSubset
   :project: gapp

