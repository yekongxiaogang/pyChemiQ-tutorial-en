:py:mod:`pychemiq.Optimizer`
============================

.. py:module:: pychemiq.Optimizer


Module Contents
---------------


Functions
~~~~~~~~~~~

.. py:function:: vqe_solver(method="NELDER-MEAD", ansatz=None, pauli=None, init_para=None, chemiq=None, Learning_rate=0.1, Xatol=0.0001, Fatol=0.0001, MaxFCalls=200, MaxIter=200)

   This class is a VQE solver, which requires specifying the classical optimizer method, ansatz, molecular Pauli Hamiltonian, initial parameter, and chemiq class in the parameters.

   :param str method: Specify the classic optimizer method. Currently, pyChemiQ supports methods such as NELDER-MEAD, POWERLL, COBYLA, L-BFGS-B, SLSQP, and Gradient-Descent. If not specified, the NELDER-MEAD optimizer is used by default.
   :param AbstractAnsatz ansatz: Specify the proposed class. Please refer to pychemiq.Circuit.Ansatz for details.
   :param PauliOperator pauli: The Pauli Hamiltonian of the specified molecule. Pauli operator class. Please refer to pychemiq.PauliOperator for details.
   :param list[float] init_para: Specify initial parameters.
   :param ChemiQ chemiq: Specify the chemiq class. Please refer to pychemiq.ChemiQ for details.
   :param float Learning_rate: Specify the learning rate. This parameter is required to select the optimizer method related to gradients. The default is 0.1.
   :param float Xatol: The convergence threshold of the variable. The default is 0.0001.
   :param float Fatol: The convergence threshold of the function value. The default is 0.0001.
   :param int MaxFCalls: The maximum number of times a function can be called. The default is 200.
   :param int MaxIter: Maximum number of optimization iterations. The default is 200.

   :return: QOptimizationResult class。for detail: `pyqpanda.QOptimizationResult <https://pyqpanda-toturial.readthedocs.io/zh/latest/autoapi/pyqpanda/index.html#pyqpanda.QOptimizationResult>`_ 。



.. note::
    For more examples of optimizers and calling external scipy. optimize libraries to achieve classic optimization, please refer to :doc:`../04advanced/optimizer` in the advanced tutorial