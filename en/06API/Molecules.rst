:py:class:`pychemiq.Molecules`
==================================

Classes
----------

.. py:class:: Molecules(geometry=None, basis=None, multiplicity=None, charge=0, active=None, nfrozen=None)

    Initialize the electronic structure parameters of molecules, including charge, base group, atomic coordinates, spin multiplicity, etc

   :param str geometry: Enter the type and coordinates of atoms in the molecule. It can be a string type or a string list. For example: geometry="H 0 00, H 0 0 0.74" or geometry=["H 0 00", "H 0 0 0.74"]
   :param str basic: Enter the base group level for performing calculations. The currently supported basis set are Gaussian function basis groups such as MINI, sto-3G, sto-6G, 3-21G, 6-31G, etc. Polarization and dispersion basis sets are not supported.
   :param int multiplicity: The spin multiplicity of the input molecular system. The relationship with the total spin quantum number of molecules is M=2S+1. Currently, pyChemiQ only supports RHF singlet calculations, and UHF and ROHF are under development.
   :param int charge: Input the charge of the molecular system.
   :param list [int] active: The active space of a molecular system, formatted as [m, n], where m is the number of active orbitals and n is the number of active electrons. By default, no active space is set.
   :param int nfrozen: The number of frozen orbitals in a molecular system, starting from the molecular orbital with the lowest energy and freezing the electrons in that orbital. The default is not to set frozen tracks.

   :return: Output the result of executing HF calculation.


   **Attributes**

   .. py:attribute:: n_atoms

      Obtaining the number of atoms in a molecular system

   .. py:attribute:: n_electrons

      Obtaining the number of electrons in the molecular system

   .. py:attribute:: n_orbitals

      Obtain the total molecular orbital number of the molecular system

   .. py:attribute:: n_qubits

      Obtain the total number of quantum bits required for calculation (i.e. number of spin orbitals, 2 * number of molecular orbitals)

   .. py:attribute:: hf_energy

      Obtain the energy calculated by HF (unit: Hartree)

   .. py:attribute:: nuclear_repulsion

      Obtaining the inter nuclear repulsion of the molecular system (unit: Hartree)

   .. py:attribute:: canonical_orbitals

      Obtain the canonical orbital coefficients (i.e. molecular orbital coefficients) of the molecular system

   .. py:attribute:: orbital_energies

      Obtaining the energy of each molecular orbital in the system
      
   .. py:attribute:: overlap_integrals

      Obtaining overlapping integrals of molecular systems

   .. py:attribute:: one_body_integrals

      Obtaining single electron integrals for molecular systems

   .. py:attribute:: two_body_integrals

      Obtaining the Double Electron Integral of a Molecular System



   **Methods**

   .. py:method:: get_molecular_hamiltonian()

      Obtaining the Hamiltonian of the initialized molecular system


---------

**Interface example:**

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

Call the following interface to obtain information about the molecular system：

.. code:: 

    print("The number of atoms is", mol.n_atoms)
    print("The number of electrons is", mol.n_electrons)
    print("The number of orbitals is", mol.n_orbitals)
    print("The number of qubits is", mol.n_qubits)
    print("The Hartree-Fock energy is", mol.hf_energy)
    print("The nuclear repulsion is", mol.nuclear_repulsion)


.. parsed-literal::

    The number of atoms is 2
    The number of electrons is 2
    The number of orbitals is 2
    The number of qubits is 4
    The Hartree-Fock energy is -1.1167593072992057
    The nuclear repulsion is 0.7151043390810812


.. code:: 

    print("The canonical orbitals are\n", mol.canonical_orbitals)
    print("The orbital energies are", mol.orbital_energies)
    print("The overlap integrals are\n", mol.overlap_integrals)


.. parsed-literal::

    The canonical orbitals are
     [[-0.54884228  1.21245192]
     [-0.54884228 -1.21245192]]
     
    The orbital energies are [-0.57855386  0.67114349]

    The overlap integrals are
     [[1.         0.65987312]
     [0.65987312 1.        ]]


.. code:: 

    print("The one body integrals are\n", mol.one_body_integrals)
    print("The two body integrals are\n", mol.two_body_integrals)


.. parsed-literal::

    The one body integrals are
     [[-1.25330979e+00  0.00000000e+00]
     [ 4.16333634e-17 -4.75068849e-01]]

    The two body integrals are
     [[[[ 6.74755927e-01 -1.11022302e-16]
       [-8.32667268e-17  6.63711401e-01]]
    
      [[-3.46944695e-17  1.81210462e-01]
       [ 1.81210462e-01  0.00000000e+00]]]
    
    
     [[[-4.85722573e-17  1.81210462e-01]
       [ 1.81210462e-01 -2.22044605e-16]]
    
      [[ 6.63711401e-01 -2.22044605e-16]
       [-1.66533454e-16  6.97651504e-01]]]]

.. code:: 

    print("The molecular hamiltonian is", mol.get_molecular_hamiltonian())


.. parsed-literal::

    The molecular hamiltonian is {
    : 0.715104
    0+ 0 : -1.253310
    1+ 0+ 1 0 : -0.674756
    1+ 0+ 3 2 : -0.181210
    1+ 1 : -1.253310
    2+ 0+ 2 0 : -0.482501
    2+ 1+ 2 1 : -0.663711
    2+ 1+ 3 0 : 0.181210
    2+ 2 : -0.475069
    3+ 0+ 2 1 : 0.181210
    3+ 0+ 3 0 : -0.663711
    3+ 1+ 3 1 : -0.482501
    3+ 2+ 1 0 : -0.181210
    3+ 2+ 3 2 : -0.697652
    3+ 3 : -0.475069
    }
    