Optimizer Tutorial
=================================

In the ansatz tutorial, we constructed a ansatz circuit with parameters, and after conducting quantum expectation estimation, we need to continuously iterate and optimize the parameters involved in Ansatz to obtain the lowest energy, and use the superposition state with the lowest energy as the ground state of the current molecular model.
The optimization of these parameters in VQE is handled using classic optimizers. As of now, pyChemiQ provides the following optimizers: NELDER-MEAD, POWER, COBYLA, L-BFGS-B, SLSQP, and Gradient-Descent.
Among them, the derivative free optimization methods are NELDER-MEAD, and COBYLA; The first order methods are L-BFGS-B, SLSQP, and Gradient-Descent. In pyChemiQ, we specify different classical optimizers by using the method in the final VQE solver.

.. code-block::

    # Using method to specify the classical optimizer and initial parameters and iteratively solve
    #method = "NELDER-MEAD"
    #method = "POWELL"
    #method = "COBYLA"
    #method = "L-BFGS-B"
    #method = "Gradient-Descent"
    method = "SLSQP"
    init_para = np.zeros(ansatz.get_para_num())
    solver = vqe_solver(
            method = method,
            pauli = pauli_H2,
            chemiq = chemiq,
            ansatz = ansatz,
            init_para=init_para)
    result = solver.fun_val
    n_calls = solver.fcalls
    print(result,f"Function called {n calls} times in total")

In addition to using the optimizer that comes with pyChemiQ, we can also call external Python libraries to implement the classic optimization part. Here, we take calling scipy. optimize as an example. Firstly, we need to use chemiq.getLossFuncValue() to obtain the loss function:

.. code-block::

    # The optimizer uses SLSQP from an external scipy library, in which case we need to first define a loss function
    def loss(para,grad,iters,fcalls):
        res = chemiq.getLossFuncValue(0,para,grad,iters,fcalls,pauli,chemiq.qvec,ansatz)
        return res[1]

    def optimScipy():
        import scipy.optimize as opt
        options = {"maxiter":2}
        init_grad = np.zeros(ansatz.get_para_num())
        method = "SLSQP"
        res = opt.minimize(loss,init_para,
                args=(init_grad,0,0),
                method = method,
                options = options)
        final_results = {
                "energy":res.fun,
                "fcalls":f"function called {res.nfev} times in total",
        }
        print(final_results)

    # Specify the Pauli Hamiltonian pauli in the loss function and the initial parameter init_para of the optimizer in the main function
    if __name__ == "__main__":
        pauli = pauli_H2
        init_para = np.zeros(ansatz.get_para_num())
        optimScipy()

If we don't want to use the existing optimizer, we can also define the optimizer ourselves. Let's take the gradient-descent method as an example. Firstly, we need to use chemiq.getExpectationValue() to obtain the expected value to define the loss function:

.. code-block::

    # The optimizer uses a custom gradient-descent method, first defining the loss function:
    def loss3(para,grad):
        new_para[:] = para[:]
        global result
        result = chemiq.getExpectationValue(0,0,0,chemiq.qvec,H,new_para,ansatz,False)
        if len(grad) == len(para):
            for i in range(len(para)):
                new_para[i] += delta_p
                result_p = chemiq.getExpectationValue(0,0,0,chemiq.qvec,H,new_para,ansatz,False)
                grad[i] = (result_p - result)/delta_p
                new_para[i] -= delta_p
        return result
    def GradientDescent():
        para_num = ansatz.get_para_num()
        seed = 20221115
        np.random.seed(seed)
        para = np.zeros(para_num)
        grad = np.zeros(para_num)
        lr = 0.1
        result_previous = 0
        result = 0
        threshold = 1e-8
        for i in range(1000):
            result_previous = result
            result = loss3(para,grad)
            para -= lr*grad
            if abs(result - result_previous) < threshold:
                print("final_energy :",result)
                print("iterations :",i+1)
                break
    # Specify the Pauli Hamiltonian H and parameter delta_p in the loss function in the main function, as well as the initial parameter new_para
    if __name__ == "__main__":
        new_para = np.zeros(ansatz.get_para_num())
        delta_p = 1e-3
        H = pauli_H2.to_hamiltonian(True)
        GradientDescent()
