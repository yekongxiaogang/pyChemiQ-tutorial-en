Introduction to Configuration File Parameters
====================================================

In addition to calling the basic interface of pyChemiQ for calculation, you can also set the configuration file for direct calculation, and more advanced functions are available for use in the configuration file. You can shorten the quantum circuit and reduce running time by using built-in optimization methods. There are also interfaces with rich functions such as slice number setting, proposed truncation, and MP2 initial parameter setting that can be called. If you want to try out these advanced features of pyChemiQ, please go to `official website <https://qcloud.originqc.com.cn/chemistryIntroduce>`_ to apply for an authorization code. The licenses for ChemiQ and pyChemiQ are universal. If you already have a license for ChemiQ, simply fill in the globally set license parameters in the configuration file. After setting the configuration file, enter the following command line on the terminal to run the calculation:

.. code-block::

    from pychemiq import directly_run_config
    directly_run_config("test.chemiq")  # Replace test with the name of your configuration file



The configuration file is usually a file with the suffix .chmiq, which is mainly set in five aspects. The detailed parameter introduction and default parameters are as follows:

1. General settings
   -task: Set the calculation type [energy/MD], that is, single point energy calculation/molecular dynamics simulation. The default value is energy.

   -backend: Set the backend type for executing calculations [CPU_SINGLE_THREAD]. Currently, only single threaded CPUs are supported, and more backend types of access are still in progress.

   -print_out: Set whether to print the scf iteration process. The default is F.

   -print_iters: Set whether to print the optimizer iteration process. The default is F.

   -console_level: Set the terminal printing log level, where 0 is output and 6 is not output. The default value is 0.

   -logfile_name: Set the log file name, which is empty by default. For example, if the log name set is chemiq.log, the final output log name is chemiq - date of the day. log. Only when the log file name is set can the log file be output.

   -logfile_level: Set the file output log level, where 0 is output and 6 is not output. The default value is 6.

   -license: Set the authorization serial number. Apply for an authorization code, Please go to the `Official-Website <https://qcloud.originqc.com.cn/chemistryIntroduce>`_.

2. Molecular specification setting
   -geoms: Set molecular coordinates, where the atomic type and atomic coordinates are separated by spaces.

   -bohr: Whether the coordinate unit is set to bohr. Boolean value. The default is F, in angstrom units.

   -charge: Set the number of charges in the system. When the number of charges is positive, there is no positive sign, and when it is negative, a negative sign is written. The default value is 0.

   -spin: Set the spin multiplicity (M=2S+1), with a default value of 1.

   -basis: Set the base group level required for calculation [MINI/to 3G/to 6G/3-21G/6-31G].

   -pure: Use spherical harmonic or Cartesian Gaussian functions. The default is T, which means using the spherical harmonic Gaussian function.

   -local: Whether to localize HF orbitals. The default is F, which means that the HF orbit is not localized.

   -active: Set the active space, separated by commas, with no active space set by default. For example, active=4, where the first 4 represents 4 spatially active orbitals, and the second 4 represents 4 electrons in the active space.

   -nfrozen: Set the number of frozen space tracks, default to 0. Note: Active and NFrozen cannot be used simultaneously.

   -mix_SCF: Setting this parameter can effectively solve the problem of SCF non convergence. The principle of the method is the damping method. Change the density matrix D (n+1) of the Fock matrix in step n+1 to w * D (n-1)+(1-w) * D (n), where w is the set parameter. The density matrix after parameter averaging weakens the difference between the current step and the previous step's density matrix, making the density matrix smoother with iteration and helping convergence. The value of this parameter is [0.0, 1.0], with a default value of 0.5.

3. Ansatz parameter settings (ansatz settings)
Set the proposed type of quantum circuit [UCC/Hardware efficient/Symmetry reserved/User defined]. Select the first three types of schemes, and the quantum circuit will be automatically generated. Select the last scheme and input the originIR format quantum circuit in the parameter circuit.

   -mapping: Set mapping [JW/P/BK/SP]. These mapping methods are Jordan Wigner Transform, Parity Transform, Bravyi Kitaev Transform, and Segment Parity Transform.

   -excluded_Level: Set the excitation level [S/D/SD], which only takes effect when ansatz is UCC.

   -restricted: Restrict the excitation term and prepare superposition states with fewer configuration wave functions to shorten the circuit. The default is T. Only effective when ansatz is UCC.

   -cutoff: Truncate the excitation term proposed for UCC based on the initial parameters of MP2. Only if ansatz is UCC and init_para Effective when type=MP2. The default is F.

   -reorder: Arrange the qubits in order, with the first half of the qubits encoding spin up and the second half encoding spin down. Setting this parameter to T can reduce the use of qubits. It takes effect when mapping is P and BK. The default is F.

   -circuit: Set the circuit through the originIR string, which only takes effect when ansatz is User define.

4. Optimizer settings
Set the classic optimizer type [Nelder Mead/Power/Gradient Descent/COBYLA/L-BFGS-B/SLSQP].

   -init_Para_Type: Set the method for constructing initial parameters [Zero/Random/input/MP2], where Zero represents the initial parameter as all zeros, Random represents the initial parameter as a random number within the [0,1) interval, input represents the custom initial parameter, and MP2 represents the initial parameter result obtained from second-order perturbation. MP2 is only available when intended for UCCD and UCCSD. The initial parameter defaults to Zero.

   -slices: Set the number of slices, i.e. the number of quantum circuit repetitions, with a default value of 1.

   -learning_Rate: Set the learning rate. The default value is 0.1.

   -iters: Set the number of iterations, with a default value of 1000.

   -fcalls: Set the number of function calls, with a default value of 1000.

   -xatol: Set the variable convergence threshold, with a default value of 1e-4.

   -fatol: Set the expected value convergence threshold, with a default value of 1e-4.

5. Molecular dynamics parameter settings

   -HF: Set the correlation sampling method. The default is 1.

   -axis: Set the system to move in a specific direction in the form of a string, in the format of "x y z".

   -save_trajectory: Set the name of the saved molecular coordinate file. The default is' traj. csv '.

   -save_topology: Set the name of the saved molecular topology file. The default is' topology. txt '.

   -velocity: Set the initial velocity of atoms, separated by commas, "0.1 0.2 0.3, -0.1-0.2-0.3 ", in units of A/fs, with default values of all 0.

   -step_size: Set the step size, greater than 0, in fs, with a default of 0.2.

   -step_number: Set the total number of steps, greater than 1, with a default of 100.

   -delta_r: Set the size of the differential coordinate, greater than 0, with a default of 0.001.

Next, we provide a case study of using a configuration file to calculate the single point energy of hydrogen molecules. The base group uses sto-3G, it is planned to use UCCSD, the mapping uses BK, and the optimizer uses NELDER MEAD. The initial reference is MP2.

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0 0 0
            H 0 0 0.74
        }
        bohr    = F
        charge  = 0
        spin    = 1 
        basis   = sto-3G
        pure    = T 
        local   = F 
    }

    ansatz = UCC {
        excited_level = SD
        restricted    = T
        cutoff        = T
        mapping       = BK
        reorder       = F
    }

    optimizer = NELDER-MEAD {
        learning_rate                 = 0.1 
        init_para_type                = MP2
        slices                        = 1 
        iters                         = 1000 
        fcalls                        = 1000 
        xatol                         = 1e-6 
        fatol                         = 1e-6 
    }


In the second example, we calculate the potential energy curve of hydrogen molecules. Here, we take scanning five points as an example, with each point spaced 0.1 angstrom apart. The base group uses sto-3G, the active space uses [2,2], it is planned to use a custom circuit, the mapping uses party, and the optimizer uses SLSQP. The initial parameter is zero.

.. code-block::

    general = {
        task    = energy
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0 0 0
            H 0 0 0.54;
            H 0 0 0
            H 0 0 0.64;
            H 0 0 0
            H 0 0 0.74;
            H 0 0 0
            H 0 0 0.84;
            H 0 0 0
            H 0 0 0.94
        }
        bohr    = F
        charge  = 0
        spin    = 1 
        basis   = sto-3G
        pure    = T 
        local   = F 
        active = 2,2
    }

    ansatz = User-define {
        circuit = {
            QINIT 4
            CREG 4
            CNOT q[1],q[0]
            CNOT q[2],q[1]
            CNOT q[3],q[2]
            H q[1]
            H q[3]
            S q[1]
    }
        mapping       = P
        reorder       = T
    }

    optimizer = SLSQP {
        learning_rate                 = 0.1 
        init_para_type                = Zero
        slices                        = 1  
        iters                         = 1000 
        fcalls                        = 1000 
        xatol                         = 1e-6 
        fatol                         = 1e-6 
    }


In the third example, we calculate the molecular dynamics trajectory of lithium hydride molecules. The base group uses 3-21G, the active space uses [4,4], hardware efficiency is proposed to be used, mapping uses JW, and optimizer uses L-BFGS-B. The initial parameter is a random number.

.. code-block::

    general = {
        task    = MD
        backend = CPU_SINGLE_THREAD
        license = XXXXX
    }

    mole = {
        geoms = {
            H 0 0 0.38
            Li 0 0 -1.13
        }
        bohr    = F
        charge  = 0
        spin    = 1 
        basis   = 3-21G
        pure    = T 
        local   = F 
        active = 4,4
    }

    ansatz = Hardware-efficient {
        mapping       = JW
        reorder       = F
    }

    optimizer = L-BFGS-B {
        learning_rate                 = 0.1 
        init_para_type                = Random
        slices                        = 1  
        iters                         = 1000 
        fcalls                        = 1000 
        xatol                         = 1e-6 
        fatol                         = 1e-6 
    }

    MD = 1 {
        velocity           = 0.0
        step_size          = 0.2
        step_number        = 100 
        delta_r            = 0.001
    }
