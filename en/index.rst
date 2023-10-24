.. pychemiq documentation master file, created by
   sphinx-quickstart on Mon Nov 14 14:20:39 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyChemiQ
====================================
PyChemiQ is a python software library developed by Origin Quantum, which can be used to simulate and calculate the fermionic system on a real Quantum computer or virtual machine. Support user-defined mapping, ansatz, and optimizers, facilitating user secondary development. This software package provides a simple, lightweight and efficient platform for Quantum chemistry calculation and method development. PyChemiQ can be used to solve Quantum chemistry problems by using the mean field and post mean field methods to simulate molecules on a Quantum computer. PyChemiQ simplifies the conversion between molecular structure input and quantum circuits, minimizes the domain expertise required to enter the field, and makes it easier for interested scholars to solve and study electronic structure problems on Quantum computer.

At present, pyChemiQ supports inputting molecular structures to derive the second-quantized Fermion Hamiltonian; In terms of mapping, pyChemiQ supports Jordan-Wigner (JW) transformation, Bravyi-Kitaev (BK) transformation, Parity transformation and Multilayer Segmented Parity (MSP)
transformation methods to map the second-quantized Fermion Hamiltonian operators to Pauli Hamiltonian operators on quantum computers; In terms of ansatz, pyChemiQ also supports different ansatz to construct quantum circuits, such as Unitary Coupled Cluster (UCC), Hardware-Efficient, symmetry-preserved, etc; In terms of optimization, pyChemiQ provides several classical optimizers for optimizing variational parameters: NELDER-MEAD, POWELL, COBYLA, L-BFGS-B, SLSQP, and Gradient Descent. Users can also construct or optimize the Quantum circuit by themselves to obtain the optimal ground state energy solution of the electronic structure problem through user-defined mapping and ansatz.

* Welcome to Fork our project on Github or leave us a message： `Github <https://github.com/OriginQ/pyChemiQ/>`_.
* To customize the pyChemiQ function for specific quantum chemistry problems, please contact Mr. Chen of origin Quantum:00-86-18221003869 or send us an email : dqa@originqc.com.
* Experience all the latest features of the visual teaching tool ChemiQ, please go to the  `official website <https://qcloud.originqc.com.cn/zh/chemistryIntroduce/>`_ to download.
* If you want to cite pyChemiQ or ChemiQ in the literature, please follow the following format: Wang Q, Liu H Y, Li Q S, et al. Chemiq: A chemistry simulator for quantum computer[J]. arXiv preprint arXiv:2106.10162, 2021.





.. toctree::
   :maxdepth: 2
   :caption: Installation Introduction
   
   01install/install.rst


.. toctree::
   :maxdepth: 2
   :caption: Quick Start
   
   02start/quickstart.rst


.. toctree::
   :maxdepth: 2
   :caption: Basic Tutorial
   
   03basis/hamiltonian.rst
   03basis/mapping.rst
   03basis/ansatz.rst
   03basis/algorithm.rst
   03basis/function.rst

.. toctree::
   :maxdepth: 2
   :caption: Advanced Tutorial
   
   04advanced/fermionpauliop.rst
   04advanced/optimizer.rst
   04advanced/circuit.rst


.. toctree::
   :maxdepth: 1
   :caption: API

   06API/index
   06API/Configs.rst


.. toctree::
   :maxdepth: 2
   :caption: Communication and Feedback
   
   07FAQ/feedback.rst
   
