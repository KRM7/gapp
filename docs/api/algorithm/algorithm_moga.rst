Multi-objective algorithms
===================================================

.. doxygennamespace:: genetic_algorithm::algorithm
   :project: gapp
   :desc-only:


class NSGA2
---------------------------------------------------

.. code-block::

   #include <algorithm/nsga2.hpp>

.. doxygenclass:: genetic_algorithm::algorithm::NSGA2
   :project: gapp
   :members:
   :protected-members:


class NSGA3
---------------------------------------------------

.. code-block::

   #include <algorithm/nsga3.hpp>

.. doxygenclass:: genetic_algorithm::algorithm::NSGA3
   :project: gapp
   :members:
   :protected-members:


Reference line generation methods
---------------------------------------------------

.. doxygennamespace:: genetic_algorithm::algorithm::reflines
   :project: gapp
   :desc-only:

.. code-block::

   #include <algorithm/reference_lines.hpp>

.. doxygenfunction:: genetic_algorithm::algorithm::reflines::quasirandomSimplexPointsMirror
   :project: gapp

.. doxygenfunction:: genetic_algorithm::algorithm::reflines::quasirandomSimplexPointsSort
   :project: gapp

.. doxygenfunction:: genetic_algorithm::algorithm::reflines::quasirandomSimplexPointsRoot
   :project: gapp

.. doxygenfunction:: genetic_algorithm::algorithm::reflines::quasirandomSimplexPointsLog
   :project: gapp

.. doxygenfunction:: genetic_algorithm::algorithm::reflines::pickSparseSubset
   :project: gapp

