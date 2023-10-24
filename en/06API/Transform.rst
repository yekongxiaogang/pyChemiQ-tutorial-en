:py:mod:`pychemiq.Transform`
===============================

Module Contents
---------------
- :py:mod:`pychemiq.Transform.mapping`  

Map the fermionic operator to the mapping submodule of the Pauli operator. When constructing parameterized quantum circuits using UCC to prepare experimental wavefunctions, we need to input the mapping type of the unitary coupling cluster operator. At this point, we specify the mapping type of the unitary coupling cluster operator through pychemiq.Transform.mapping.MappingType. The mapping type here needs to be consistent with the mapping method of the Hamiltonian, which means we need to maintain the same set of mapping methods in the calculation.

Classes
~~~~~~~~~~~

.. py:module:: pychemiq.Transform.mapping

   .. py:class:: MappingType

     The MappingType enumeration class has four different values, namely:

      .. py:attribute:: Bravyi_Kitaev

      .. py:attribute:: Jordan_Wigner

      .. py:attribute:: Parity

      .. py:attribute:: SegmentParity

---------


**Interface Example：**

.. code:: 

      from pychemiq import ChemiQ,QMachineType
      from pychemiq.Transform.Mapping import MappingType
      from pychemiq.Circuit.Ansatz import UCC

      chemiq = ChemiQ()
      machine_type = QMachineType.CPU_SINGLE_THREAD
      # Using JW mapping method
      mapping_type = MappingType.Jordan_Wigner
      # Using BK mapping method
      # mapping_type = MappingType.Bravyi_Kitaev
      # Using Parity mapping method
      # mapping_type = MappingType.Parity
      # Using SegmentParity mapping method
      # mapping_type = MappingType.SegmentParity
      chemiq.prepare_vqe(machine_type,mapping_type,2,1,4)
      ansatz = UCC("UCCSD",2,mapping_type,chemiq=chemiq)


Functions
~~~~~~~~~~~

   .. py:function:: jordan_wigner(fermion)

      Map the input fermionic operator into a Pauli operator through the Jordan Wigner transformation.

      :param FermionOperator fermion: Enter the fermionic operator to be mapped.

      :return: The mapped Pauli operator. Pauli operator class.



   .. py:function:: bravyi_kitaev(fermion)

      Map the input fermionic operator into a Pauli operator through the Bravyi Kitaev transformation.

      :param FermionOperator fermion: Enter the fermionic operator to be mapped.

      :return: The mapped Pauli operator. Pauli operator class.



   .. py:function:: parity(fermion)

      Map the input fermionic operator into a Pauli operator through Parity transformation.

      :param FermionOperator fermion: Enter the fermionic operator to be mapped.

      :return:  The mapped Pauli operator. Pauli operator class.



   .. py:function:: segment_parity(fermion)

      Pass the input fermionic operator through segment_party transformation maps to the Pauli operator.

      :param FermionOperator fermion: Enter the fermionic operator to be mapped.

      :return: The mapped Pauli operator. Pauli operator class.


---------


**Interface example：**

In the following example, we use the above four mapping methods to map the Hamiltonian of the hydrogen molecule after second-quantized from the fermionic operator to the form of the Pauli operator. Firstly, initialize the electronic structure parameters of the molecule to obtain the Hamiltonian in fermionic form.

.. code::

      from pychemiq import Molecules

      multiplicity = 1
      charge = 0
      basis =  "sto-3g"
      geom = "H 0 0 0,H 0 0 0.74"

      mol = Molecules(
      geometry = geom,
      basis    = basis,
      multiplicity = multiplicity,
      charge = charge)
      fermion_H2 = mol.get_molecular_hamiltonian()

Obtain the Hamiltonian of hydrogen molecules in Pauli form through JW transformation and print the results.

.. code::

      from pychemiq.Transform.Mapping import jordan_wigner
      pauli_H2 = jordan_wigner(fermion_H2)
      print(pauli_H2)

.. code::

      {
      "" : -0.097066,
      "X0 X1 Y2 Y3" : -0.045303,
      "X0 Y1 Y2 X3" : 0.045303,
      "Y0 X1 X2 Y3" : 0.045303,
      "Y0 Y1 X2 X3" : -0.045303,
      "Z0" : 0.171413,
      "Z0 Z1" : 0.168689,
      "Z0 Z2" : 0.120625,
      "Z0 Z3" : 0.165928,
      "Z1" : 0.171413,
      "Z1 Z2" : 0.165928,
      "Z1 Z3" : 0.120625,
      "Z2" : -0.223432,
      "Z2 Z3" : 0.174413,
      "Z3" : -0.223432
      }

Obtain the Hamiltonian of hydrogen molecules in Pauli form through BK transformation and print the results.

.. code::

      from pychemiq.Transform.Mapping import bravyi_kitaev
      pauli_H2 = bravyi_kitaev(fermion_H2)
      print(pauli_H2)

.. code::

      {
      "" : -0.097066,
      "X0 Z1 X2" : 0.045303,
      "X0 Z1 X2 Z3" : 0.045303,
      "Y0 Z1 Y2" : 0.045303,
      "Y0 Z1 Y2 Z3" : 0.045303,
      "Z0" : 0.171413,
      "Z0 Z1" : 0.171413,
      "Z0 Z1 Z2" : 0.165928,
      "Z0 Z1 Z2 Z3" : 0.165928,
      "Z0 Z2" : 0.120625,
      "Z0 Z2 Z3" : 0.120625,
      "Z1" : 0.168689,
      "Z1 Z2 Z3" : -0.223432,
      "Z1 Z3" : 0.174413,
      "Z2" : -0.223432
      }

Obtain the Hamiltonian of hydrogen molecules in Pauli form through Parity transformation and print the results.

.. code::

      from pychemiq.Transform.Mapping import parity
      pauli_H2 = parity(fermion_H2)
      print(pauli_H2)


.. code::

      {
      "" : -0.097066,
      "X0 Z1 X2" : 0.045303,
      "X0 Z1 X2 Z3" : 0.045303,
      "Y0 Y2" : 0.045303,
      "Y0 Y2 Z3" : 0.045303,
      "Z0" : 0.171413,
      "Z0 Z1" : 0.171413,
      "Z0 Z1 Z2" : 0.120625,
      "Z0 Z1 Z2 Z3" : 0.120625,
      "Z0 Z2" : 0.165928,
      "Z0 Z2 Z3" : 0.165928,
      "Z1" : 0.168689,
      "Z1 Z2" : -0.223432,
      "Z1 Z3" : 0.174413,
      "Z2 Z3" : -0.223432
      }

Obtain the Hamiltonian of hydrogen molecules in Pauli form through SP transformation and print the results.

.. code::

      from pychemiq.Transform.Mapping import segment_parity
      pauli_H2 = segment_parity(fermion_H2)
      print(pauli_H2)

.. code::

      {
      "" : -0.097066,
      "X0 Z1 X2" : 0.045303,
      "X0 Z1 X2 Z3" : 0.045303,
      "Y0 Z1 Y2" : 0.045303,
      "Y0 Z1 Y2 Z3" : 0.045303,
      "Z0" : 0.171413,
      "Z0 Z1" : 0.171413,
      "Z0 Z1 Z2" : 0.165928,
      "Z0 Z1 Z2 Z3" : 0.165928,
      "Z0 Z2" : 0.120625,
      "Z0 Z2 Z3" : 0.120625,
      "Z1" : 0.168689,
      "Z1 Z2 Z3" : -0.223432,
      "Z1 Z3" : 0.174413,
      "Z2" : -0.223432
      }
