Utility classes and functions
===================================================

Random number generation
---------------------------------------------------

.. code-block::

   #include <utility/rng.hpp>

.. doxygennamespace:: gapp::rng
   :project: gapp
   :members:
   :undoc-members:


Floating-point comparison tolerances
---------------------------------------------------

.. code-block::

   #include <utility/math.hpp>

.. doxygenclass:: gapp::math::ScopedTolerances
   :project: gapp
   :members:
   :protected-members:
   :private-members:


.. code-block::

   #include <utility/math.hpp>

.. doxygenclass:: gapp::math::Tolerances
   :project: gapp
   :members:
   :protected-members:
   :private-members:


Bounded value types
---------------------------------------------------

.. code-block::

   #include <utility/bounded_value.hpp>

.. doxygentypedef:: gapp::BoundedValue
   :project: gapp

.. doxygentypedef:: gapp::NonNegative
   :project: gapp

.. doxygentypedef:: gapp::Negative
   :project: gapp

.. doxygentypedef:: gapp::Positive
   :project: gapp

.. doxygentypedef:: gapp::Probability
   :project: gapp

.. doxygentypedef:: gapp::Normalized
   :project: gapp

