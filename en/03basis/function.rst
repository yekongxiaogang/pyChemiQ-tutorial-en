Functional Tutorials
=================================

We all know that the essence of chemical reactions is the process of breaking old chemical bonds and forming new ones. The fracture(breaking) and formation(forming) of chemical bonds are only qualitative descriptions of them at the macro level.
In fact, the reason why we called chemical bonds "forming" is because the atoms between molecules approach each other along the direction of the reaction coordinates, and when they are in a certain relative position, the energy of the entire system is the lowest and the most stable. The "breaking" of chemical bonds is when the atoms between molecules overcome the energy barrier and move away from each other along the direction of the reaction coordinate after external energy input. In computational chemistry, we can use Potential Energy Surface (PES) to describe the energy of molecules at different positions of their atoms, thereby helping us understand the process of chemical reactions.

.. image:: ./picture/PES_diatom.png
   :align: center
   :scale: 40%
.. centered:: Figure 1: Potential energy surface of bimolecular atoms

For example, as a certain bond within a molecule grows, the energy will change. A curve of energy bond length change is called the potential energy curve, as shown in Figure 1; If you make an image of the potential energy of a molecule changing with two coordinate parameters, you will find that this is a surface (because there are three quantities: two coordinate variables plus energy, forming a three-dimensional space), which is called the potential energy surface, as shown in Figure 2; By analogy, the entire molecular potential energy varies with all possible atomic coordinate variables, forming a complex hypersurface in multidimensional space, collectively known as the potential energy surface [1]_.

.. image:: ./picture/PES.png
   :align: center
   :scale: 80%
.. centered::  Figure 2: Water molecule potential energy surface: The lowest point of potential energy corresponds to the optimized water molecule structure, with an O-H bond length of 0.0958 nm and an H-O-H angle of :math:`104.5^{\circ}`. Figure cited from [2]_

In this tutorial, we will demonstrate how to use pyChemiQ to obtain the potential energy surface data of molecules and use matplotlib to draw the potential energy surface.
Here, we take the diatomic molecule :math:`H_2` as an example, using JW transformation for mapping and ansatz to utilizes UCCSD and SLSQP for the classical optimizer. Finally, we compare the potential energy curve of pyChemiQ with that of PySCF:

.. code-block::

    # Import the required package
    from pychemiq import Molecules,ChemiQ,QMachineType
    from pychemiq.Transform.Mapping import jordan_wigner,MappingType
    from pychemiq.Optimizer import vqe_solver
    from pychemiq.Circuit.Ansatz import UCC
    import numpy as np
    from pyscf import gto, scf, fci
    import matplotlib.pyplot as plt

    # Perform potential energy surface scanning: initialize parameters first, then construct molecular systems with different bond lengths, and perform multiple single point energy calculations
    basis = 'sto-3g'
    multiplicity = 1
    charge=0

    ## Define step spacing and number of steps
    bond_length_interval = 0.1
    n_points = 40
    bond_lengths = []
    energies = []
    for point in range(3, n_points + 1):
        bond_length = bond_length_interval * point
        bond_lengths += [bond_length]
        geometry = ["H 0 0 0", f"H 0 0 {bond_length}"]
    
        mol = Molecules(
            geometry = geometry,
            basis    = basis,
            multiplicity = multiplicity,
            charge = charge)
    
        fermion_H2 = mol.get_molecular_hamiltonian()
        pauli_H2 = jordan_wigner(fermion_H2)
    
        chemiq = ChemiQ()
        machine_type = QMachineType.CPU_SINGLE_THREAD
        mapping_type = MappingType.Jordan_Wigner
        pauli_size = len(pauli_H2.data())
        n_qubits = mol.n_qubits
        n_elec = mol.n_electrons
        chemiq.prepare_vqe(machine_type,mapping_type,n_elec,pauli_size,n_qubits)
        ansatz = UCC("UCCSD",n_elec,mapping_type,chemiq=chemiq)
    
        method = "SLSQP"
        init_para = np.zeros(ansatz.get_para_num())
        solver = vqe_solver(
                method = method,
                pauli = pauli_H2,
                chemiq = chemiq,
                ansatz = ansatz,
                init_para=init_para)
        energy = solver.fun_val
        energies += [energy]

    # Using the FCI method of classical computational chemistry software PySCF to calculate the energy of hydrogen molecules at different bond lengths
    pyscf_energies = []
    bond_length_interval = 0.1
    n_points = 40
    for point in range(3, n_points + 1):
        bond_length = bond_length_interval * point
        atom = f'H 0 0 0; H 0 0 {bond_length}'
    
        mol = gto.M(atom=atom,   # in Angstrom
                basis='STO-3G',
                charge=0,
                spin=0)
        myhf = scf.HF(mol).run() 
        cisolver = fci.FCI(myhf) 
        pyscf_energies += [cisolver.kernel()[0]]

    # Using the FCI method of classical computational chemistry software PySCF to calculate the energy of hydrogen molecules at different bond lengths
    plt.figure()
    plt.plot(bond_lengths, energies, '-g',label='pyChemiQ')
    plt.plot(bond_lengths, pyscf_energies, '--r',label='PySCF')
    plt.ylabel('Energy in Hartree')
    plt.xlabel('Bond length in angstrom')
    plt.legend()
    plt.show()

The comparison of the obtained hydrogen molecular potential energy maps is shown in the following figure. Due to the close calculation results of the two, most of the potential energy surfaces are in a state of overlap.

.. image:: ./picture/PES_H2.png
   :align: center
   :scale: 8%
.. centered:: Figure 3: Hydrogen molecular potential energy surface obtained from pyChemiQ and PySCF
















**References**

.. [1]  Baidu. https://baike.baidu.com/item/%E5%8A%BF%E8%83%BD%E9%9D%A2/6295493, last access on 6th January, 2023
.. [2]  Wikipedia. Potential energy surface. https://en.wikipedia.org/wiki/Potential_energy_surface, last access on 6th January, 2023
