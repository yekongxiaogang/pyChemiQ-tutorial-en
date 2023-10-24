:py:class:`pychemiq.PauliOperator`
=========================================

Classes
----------

.. py:class:: PauliOperator({pauli_string: coefficient})

   :param str pauli_string: The Pauli operator in character form.
   :param float coefficient: The coefficient of the Pauli operator for this term.

   :return: Pauli operator class.

   .. py:function:: get_max_index()

   Obtain the maximum index value. If it is an empty Pauli operator item calling the get_max_index() interface, it returns SIZE_MAX (specific value depends on the operating system), otherwise return its maximum index value.

   .. py:function:: data()

   The Pauli operator class provides a data interface that can return data maintained internally by Pauli operators.

.. note::
    For a detailed introduction of this class, please refer to :doc:`../04advanced/fermionpauliop` in the advanced tutorial.
   
