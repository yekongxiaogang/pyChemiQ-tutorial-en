:py:mod:`pychemiq.Circuit`
==============================

.. py:module:: pychemiq.Circuit


Module Contents
---------------
- :py:mod:`pychemiq.Circuit.Ansatz`  

Construct Ansatz submodule for quantum circuit ansatz.


Functions
~~~~~~~~~~~

.. py:module:: pychemiq.Circuit.Ansatz

   .. py:function:: UCC(ucc_type, n_electrons, mapping_type, chemiq=None)

         Construct a quantum circuit ansatz using the unitary coupled cluster operator.


        :param str ucc_type: The excitation level of the input unitary coupling cluster. Currently available: UCCS, UCCD, UCCSD.
        :param int n_electrons: Input the number of electrons in the molecular system.
        :param MappingType mapping_type: Input the mapping type of the unitary coupling cluster operator. Please refer to pychemiq.Transform.Mapping for details.
        :param ChemiQ chemiq: Specify the chemiq class. Please refer to pychemiq.ChemiQ for details.

        :return: Output the AbstractAnsatz class with the specified excitation level.


   .. py:function:: HardwareEfficient(n_electrons, chemiq=None)

         Construct quantum circuit ansatz using HardwareEfficient.

        :param int n_electrons: Input the number of electrons in the molecular system.
        :param ChemiQ chemiq: Specify the chemiq class. Please refer to pychemiq.ChemiQ for details.

        :return: Output the proposed AbstractAnsatz class.


   .. py:function:: SymmetryPreserved(n_electrons, chemiq=None)

         Construct quantum circuit ansatz using SymmetryPreserved.


        :param int n_electrons: Input the number of electrons in the molecular system.
        :param ChemiQ chemiq: Specify the chemiq class. Please refer to pychemiq.ChemiQ for details.

        :return: Output the proposed Abstract Ansatz class.



   .. py:function:: UserDefine(n_electrons, circuit=None, fermion=None, chemiq=None)

        Construct quantum circuit ansatz using user-defined methods.


        :param int n_electrons: Input the number of electrons in the molecular system.
        :param str circuit: Construct the originIR string of the quantum circuit.
        :param FermionOperator fermion: The fermionic operator class for constructing quantum circuits.
        :param ChemiQ chemiq: Specify the chemiq class. Please refer to pychemiq.ChemiQ for details.

        :return: Output the customized AbstractAnsatz class.


.. note::
    For detailed call examples of the first three functions of the Ansatz module, please refer to :doc:`../03basis/ansatz` in the basic tutorial. For an example of calling the last function, please refer to :doc:`../04advanced/circuit` in the advanced tutorial




