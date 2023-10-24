:py:mod:`pychemiq.Utils`
============================

.. py:module:: pychemiq.Utils


Module Contents
---------------

Functions
~~~~~~~~~~~

.. py:function:: get_cc_n_term(n_qubits, n_elec, excited_level)

   Obtain the number of coupling cluster operator terms at the specified excitation level. For example, for four quantum bits, the single and double excited coupling cluster operator of a 2-electron system, with spin orbitals 0 and 1 as occupied states, has five coupling cluster terms: 0->2,0->3,1->2,1->3,01->23.

   :param int n_qubits: Enter the number of qubits required for calculation.
   :param int n_elec: Enter the number of electrons in the molecular system.
   :param str excited_level: The excitation level of the input coupling cluster operator. Currently available: CCS, CCD, CCSD.

   :return: Output the number of coupling cluster operator terms at the specified excitation level. Integer type.



.. py:function:: get_cc(n_qubits, n_elec, para, excited_level='SD')

   Obtain the coupling cluster operator for the specified excitation level with parameters. For example, for four qubits, the single and double excited coupling cluster operator of a 2-electron system has spin orbitals 0 and 1 as occupied states, so the excited coupling cluster term is: 0->2,0->3,1->2,1->3,01->23. The output Fermion operators are: {{"2+0": para [0]}, {"3+0": para [1]}, {"2+1": para [2]}, {"3+1": para [3]}, {"3+2+10": para [4]}

   :param int n_qubits: Enter the number of qubits required for calculation.
   :param int n_elec: Enter the number of electrons in the molecular system.
   :param list[float] para: Enter the initial parameter list.
   :param str excited_level: The excitation level of the input coupling cluster operator. Currently available: CCS, CCD, CCSD. The default is the single double excited coupling cluster operator (CCSD).

   :return: Output the coupling cluster operator for the specified excitation level. Fermion operator class.





.. py:function:: transCC2UCC(Pauli)

   Only unitary operators can be placed on quantum circuits for simulation. On the basis of the coupling cluster operator, this function removes the operator that is not a Hermitian matrix and constructs a "unitary operator version" of the coupling cluster operator.

   :param PauliOperator Pauli: Input coupling cluster operator. Pauli operator class.

   :return: Output the unitary coupling cluster operator. Pauli operator class.
   

---------

**Interface example:**

In the following example, we calculate the single and double excited coupling cluster operators for 4 quantum bits and 2 electron systems, and convert them into unitary coupling cluster operators that can directly construct quantum circuits.

.. code::

    from pychemiq.Utils import get_cc_n_term,get_cc,transCC2UCC
    from pychemiq.Transform.Mapping import jordan_wigner
    import numpy as np

    # Calculate the number of terms for the single and double excited coupling cluster operators of 4 quantum bits and 2 electron systems to initialize the parameter list. Here, we first set the initial parameters to a list of all 1
    n_para = get_cc_n_term(4,2,"CCSD")
    para = np.ones(n_para)

    # After printing and setting the initial parameters, the single and double excited coupling cluster operator for a 4-qubit, 2-electron system
    cc_fermion = get_cc(4,2,para,"CCSD")
    print(cc_fermion)

The printed result is：

.. code:: 

    {
    2+ 0 : 1.000000
    3+ 0 : 1.000000
    2+ 1 : 1.000000
    3+ 1 : 1.000000
    3+ 2+ 1 0 : 1.000000
    }


.. code:: 

    # Next, use JW mapping to map the Fermion operator to the Pauli operator
    cc_pauli = jordan_wigner(cc_fermion)
    # Delete the non unitary coupling cluster operator and leave the unitary coupling cluster operator behind
    ucc_pauli = transCC2UCC(cc_pauli)
    print(ucc_pauli)

The printed result is:

.. code:: 

    {
    "X0 X1 X2 Y3" : -0.125000,
    "X0 X1 Y2 X3" : -0.125000,
    "X0 Y1 X2 X3" : 0.125000,
    "X0 Y1 Y2 Y3" : -0.125000,
    "X0 Z1 Y2" : 0.500000,
    "X0 Z1 Z2 Y3" : 0.500000,
    "X1 Y2" : 0.500000,
    "X1 Z2 Y3" : 0.500000,
    "Y0 X1 X2 X3" : 0.125000,
    "Y0 X1 Y2 Y3" : -0.125000,
    "Y0 Y1 X2 Y3" : 0.125000,
    "Y0 Y1 Y2 X3" : 0.125000,
    "Y0 Z1 X2" : -0.500000,
    "Y0 Z1 Z2 X3" : -0.500000,
    "Y1 X2" : -0.500000,
    "Y1 Z2 X3" : -0.500000
    }

