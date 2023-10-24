:py:class:`pychemiq.FermionOperator`
=========================================

Classes
----------

.. py:class:: FermionOperator({fermion_string: coefficient})

   :param str fermion_string: The Fermi operator in character form.
   :param float coefficient: The coefficient of this Fermi operator.

   :return: Fermi operator class.


   .. py:function:: normal_ordered()

   The normal\_ordered interface organizes fermionic operators. In this transformation, it is specified that the orbital encoding applied is sorted from high to low, and the generated operator appears before the annihilation operator.

   .. py:function:: data()

   The fermionic operator class also provides a data interface that can return data maintained internally by the fermionic operator.

.. note::
    For a detailed introduction to this class, please refer to the :doc:`../04advanced/fermionpauliop` in the advanced tutorial.
   
