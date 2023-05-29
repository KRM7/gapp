Single-objective algorithms
===================================================

.. toctree::
   :maxdepth: 1

   algorithm_selection
   algorithm_replacement


class SingleObjective
---------------------------------------------------

.. code-block::

   #include <algorithm/single_objective.hpp>

.. doxygennamespace:: genetic_algorithm::algorithm
   :project: gapp
   :desc-only:

.. doxygenclass:: genetic_algorithm::algorithm::SingleObjective
   :project: gapp
   :members:
   :protected-members:


class Selection
---------------------------------------------------

.. code-block::

   #include <algorithm/selection_base.hpp>

.. doxygennamespace:: genetic_algorithm::selection
   :project: gapp
   :desc-only:

.. doxygenclass:: genetic_algorithm::selection::Selection
   :project: gapp
   :members:
   :protected-members:
   :private-members:


class Replacement
---------------------------------------------------

.. code-block::

   #include <algorithm/replacement_base.hpp>

.. doxygennamespace:: genetic_algorithm::replacement
   :project: gapp
   :desc-only:

.. doxygenclass:: genetic_algorithm::replacement::Replacement
   :project: gapp
   :members:
   :protected-members:
   :private-members:
