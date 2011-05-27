*****************************************************
Resource requirement
*****************************************************


.. highlight:: text

``pmap`` output

.. literalinclude:: examples/map1
  :language: text
  :linenos:

Lines 2-4
  Memory requirement of the actual script. May vary depending on the complexity of the system.

Line 5
  Simulation data (state vector, density operator, etc.). Depends on the dimensionality of the system.

Line 6-85
  Shared libraries. Constant.

.. literalinclude:: examples/map2
  :language: text
  :linenos:
  :lines: 1-5

.. literalinclude:: examples/map3
  :language: text
  :linenos:
  :lines: 1-5

.. literalinclude:: examples/map4
  :language: text
  :linenos:
  :lines: 1-5

A somewhat more complex script takes slightly more memory (Line 2):

.. literalinclude:: examples/map5
  :language: text
  :linenos:
  :lines: 1-6

.. literalinclude:: examples/map6
  :language: text
  :linenos:
  :lines: 1-6

.. literalinclude:: examples/map7
  :language: text
  :linenos:
  :lines: 1-6

One of the most complex scripts used so far:

.. literalinclude:: examples/map8
  :language: text
  :linenos:
  :lines: 1-6


.. highlight:: c++
  :linenothreshold: 10

