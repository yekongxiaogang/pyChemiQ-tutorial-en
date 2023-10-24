Operator Class Tutorial
=================================

fermionic operator class
----------------------------------

We have already introduced the creation operator :math:`a^{\dagger}_p` and annihilation operator :math:`a_q` of Fermions in the first section of the basic tutorial. If we want to customize the Hamiltonian or ansatz, we can use the fermionic operator class and Pauli operator class in pyChemiQ.
In the fermionic operator class of pyChemiQ,  :math:`P+` represents the creation operator :math:`a^\dagger_p`, and :math:`Q` represents the annihilation operator :math:`a_q`. For example, "1+3 5+1" represents: :math:`a^\dagger_1 a_3 a^\dagger_5 a_1`. The Fermioperator class can also store the coefficients of this term simultaneously, such as constructing Fermi subterms  :math:`2a_0a^\dagger_1` 、  :math:`3a^\dagger_2a_3` and :math:`a^\dagger_1 a_3 a^\dagger_5 a_1` :

.. code-block::

   from pychemiq import FermionOperator
   a = FermionOperator("0 1+", 2)  
   b = FermionOperator("2+ 3", 3)
   c = FermionOperator("1+ 3 5+ 1",1)

To construct the addition and subtraction of multiple Fermi subterms, construct multiple expressions in lexicographic order. For example, constructing  :math:`2a_0a^\dagger_1+3a^\dagger_2a_3` with p0. In addition, we can also construct an empty fermionic operator class p1, which does not contain any production and annihilation operators or constant terms; Alternatively, only construct the repulsive energy constant term between electrons, such as p2; You can also construct a copy of the already constructed fermionic operator, such as p3.

.. code-block::

   p0 = FermionOperator({"0 1+": 2, "2+ 3": 3})
   p1 = FermionOperator()
   p2 = FermionOperator(2)
   p3 = p2

The FermioOperator class also provides basic operations for adding, subtracting, and multiplying fermionic operators. The result returned by the calculation is still a fermionic operator class. This class also supports printing function, which can directly print and output the fermionic operator class to the screen.

.. code-block::

   plus = a + b
   minus = a - b
   multiply = a * b
   print("a + b = {}".format(plus))
   print("a - b = {}".format(minus))
   print("a * b = {}".format(multiply))

By using the above example code, the calculation result of  :math:`a+b` ， :math:`a-b` and :math:`a*b` is as follows:

.. code-block::

   a + b = {
   0 1+ : 2.000000
   2+ 3 : 3.000000
   }
   a - b = {
   0 1+ : 2.000000
   2+ 3  : -3.000000
   }
   a * b = {
   0 1+ 2+ 3 : 6.000000
   }

You can also organize fermionic operators through the normal\_ordered interface. In this transformation, it is specified that the orbital encoding applied is sorted from high to low, and the creation operator appears before the annihilation operator. The organizing rules are as follows: for two identical numbers, exchanging annihilation and creation operators is equivalent to unit 1 minus normal order. If there are two identical creation or annihilation operators at the same time, the expression for this term is equal to 0 (Pauli exclusion principle); For different numbers, to organize them into a normal order, it is necessary to change the coefficients of the normal order to the opposite number (the anticommutative relationship of the fermionic operator). For example, the expression "0 1+2+3" is organized into a normal order of "2+1+3 0", which is equivalent to exchanging different numbers 4 times without changing the coefficients.

.. code-block::

   print(multiply.normal_ordered())
   {
   2+ 1+ 3 0 : 6.000000
   }

In addition, the fermionic operator class also provides a data interface that can return data maintained internally by the fermionic operator.

.. code-block::

   print("data = {}".format(a.data()))
   data = [(([(0, False), (1, True)], '0 1+'), (2+0j))]


Pauli operator class
----------------------------------

The Pauli operator is a set of three :math:`2×2` unitary Hermitian complex matrices, also known as unitary matrices. We usually use the Greek letter :math:`\sigma` to represent, denoted as  :math:`\sigma_x` ， :math:`\sigma_y` ， :math:`\sigma_z`. In quantum circuits, we call them X gate, Y gate, and Z gate. Their corresponding matrix forms are as follows:

Pauli-X gate：

.. math::
   X=\sigma_x=\begin{bmatrix} 0 & 1\\ 1 & 0 \end{bmatrix}
   
Pauli-Y门：

.. math::
   Y=\sigma_y=\begin{bmatrix} 0 & -i\\ i & 0 \end{bmatrix}

Pauli-Z门：

.. math::
   Z=\sigma_z=\begin{bmatrix} 1 & 0\\ 0 & -1 \end{bmatrix}

The above three Pauli matrices are sometimes referred to as spin matrices. It can be observed that the Pauli-X gate is equivalent to the NOT gate, The function of the Pauli-Y gate is equivalent to rotating around the :math:`Y`-axis of the Bloch ball at an angle of :math:`\pi`. The function of the Pauli-Z gate is equivalent to a rotation angle of :math:`\pi` around the :math:`Z`-axis of the Bloch ball.

**The Pauli operator has the following properties:**

1. The Pauli operator multiplied by itself yields the identity matrix

.. math::
    &\sigma_x \sigma_x=I \\
		&\sigma_y \sigma_y=I \\
		&\sigma_z \sigma_z=I

2. The relationship between the two Pauli operators multiplied in sequence and the Pauli operators not involved in the calculation is :math:`i`-fold

.. math::
   &\sigma_x \sigma_y=i \sigma_z \\
    	&\sigma_y \sigma_z=i \sigma_x \\
    	&\sigma_z \sigma_x=i \sigma_y 

3. The relationship between the two Pauli operators multiplied in reverse order and the Pauli operators not involved in the calculation is :math:`-i` times

.. math::
   &\sigma_y \sigma_x=-i \sigma_z \\
			&\sigma_z \sigma_y=-i \sigma_x \\
			&\sigma_x \sigma_z=-i \sigma_y 

The PauliOperator class is implemented in pyChemiQ. We can easily construct Pauli operator classes, such as constructing an empty Pauli operator term, such as p1; Alternatively, construct the direct product term :math:`2\sigma_z^0\sigma_z^1` of the Pauli operator with coefficients, such as p2. The number in the upper right corner of the Pauli operator represents the specific qubit acting on it. This term represents a Pauli Z gate acting on qubit 0 multiplied by a Pauli Z gate acting on qubit 1, with a coefficient of 2; If you want to construct the sum of multiple direct product terms of Pauli operators, you can use the form of dictionary order, such as p3 constructing :math:`2\sigma_z^0\sigma_z^1 + 3\sigma_x^1\sigma_y^2`; Alternatively, construct an identity matrix such as p4 with a coefficient of 5, or construct it in the form of p5, which is equivalent.

.. code-block::

   from pychemiq import PauliOperator
   p1 = PauliOperator()
   p2 = PauliOperator("Z0 Z1", 2)
   p3 = PauliOperator({"Z0 Z1": 2, "X1 Y2": 3})
   p4 = PauliOperator(5)
   p5 = PauliOperator("", 5)

**note:**  *When constructing the Pauli operator class, the characters contained in the string can only be one or more of spaces, X, Y, and Z, and containing other characters will throw an exception. In addition, the bit index of the same Pauli operator in the same string cannot be the same, for example: PauliOperator ("Z0 Z0", 2) will throw an exception* 。

Like the fermionic operator class, Pauli operator classes can perform addition, subtraction, multiplication, and other operations, and the returned result is still a Pauli operator class. And it also supports the printing function, where we can print and output the Pauli operator class to the screen for easy viewing of its values.

.. code-block::

   a = PauliOperator("Z0 Z1", 4)
   b = PauliOperator("X5 Y6", 3)
   plus = a + b
   minus = a - b
   muliply = a * b
   print(plus)

In practical use, we often need to know how many qubits the Pauli operator term operates on, and at this point, we obtain it by calling the get_max_index() interface of the Pauli operator class. If it is an empty Pauli operator term calling the get_max_index() interface, it returns SIZE_MAX (depending on the operating system), otherwise it returns its maximum index value. In the following example, the former outputs a value of 1, while the latter outputs a value of 6.

.. code-block::

   a = PauliOperator("Z0 Z1", 2)
   b = PauliOperator("X5 Y6", 3)
   print(a.get_max_index())
   print(b.get_max_index())


In addition, the Pauli operator class also provides a data interface that can return data maintained internally by the Pauli operator.

.. code-block::

   print("data = {}".format(a.data()))
   data = [(({0: 'Z', 1: 'Z'}, 'Z0 Z1'), (2+0j))]


