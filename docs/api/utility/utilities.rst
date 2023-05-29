Utility classes and functions
===================================================

Random number generation
---------------------------------------------------

.. code-block::

   #include <utility/rng.hpp>

.. doxygennamespace:: genetic_algorithm::rng
   :project: gapp
   :members:
   :undoc-members:


Floating-point comparison tolerances
---------------------------------------------------

.. code-block::

   #include <utility/math.hpp>

.. doxygenclass:: genetic_algorithm::math::ScopedTolerances
   :project: gapp
   :members:
   :protected-members:
   :private-members:


.. code-block::

   #include <utility/math.hpp>

.. doxygenclass:: genetic_algorithm::math::Tolerances
   :project: gapp
   :members:
   :protected-members:
   :private-members:


Bounded value types
---------------------------------------------------

.. code-block::

   #include <utility/bounded_value.hpp>

.. doxygentypedef:: genetic_algorithm::BoundedValue
   :project: gapp

.. doxygentypedef:: genetic_algorithm::NonNegative
   :project: gapp

.. doxygentypedef:: genetic_algorithm::Negative
   :project: gapp

.. doxygentypedef:: genetic_algorithm::Positive
   :project: gapp

.. doxygentypedef:: genetic_algorithm::Probability
   :project: gapp

.. doxygentypedef:: genetic_algorithm::Normalized
   :project: gapp

