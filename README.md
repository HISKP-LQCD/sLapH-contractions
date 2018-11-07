@mainpage

This is the new version of the contraction code. 


[GitHub Page](https://github.com/maowerner/sLapH-contractions)

[Doxygen Documentation](https://hiskp-lqcd.github.io/sLapH-contractions/)

[![Build
Status](https://travis-ci.org/HISKP-LQCD/sLapH-contractions.svg?branch=master)](https://travis-ci.org/maowerner/sLapH-contractions)


A Wick contraction in this context means the analytic result obtained when
"integrating out" fermionic degrees of freedom from the path integral in 
lattice QCD. Therefor it builds the bridge between the results of the numerical 
integration and the physical observables as described in section sec:intro 
by computing the euclidean correlator. 


In the sLapH framework space is approximated by the span of the eigenvectors
of a Laplace operator, below referred to as eigenvector space. Introducing
random vectors and dilution in time, eigenvector and Dirac space yields a 
stochastical estimate for the propagator. 
The correlation function can be expressed as a trace of products of 
propagators, eigenvectors and field operators (@f$\Gamma@f$). 

In the framework of the LapHs-method, the quark lines could be expressed by 
the @f$ V_s @f$ whose dimension was the number of LapH eigenvectors. Thus for the 
LapHs-method it is possible
to only apply the random vectors in the LapH-subspce spanned by spin, time
and LapH eigenvector number. To distinguish these random vectors from the 
usual ones, the are denoted by @f$ \rho @f$. For the Quark lines this results in
@f{align}{
  \mathcal{Q}
    = \sum_b { E \left( V_s V_s^\dagger \Omega^{-1} V_s P^{(b)} \rho
        ( V_s P^{(b)} \rho)^\dagger \right)}.
@f}

To illustrate these rather abstract concepts, the example of  a single charged 
pion shall be discussed. Its quantum numbers are @f$ J^P = 0^- @f$, a 
pseudoscaler particle. The operator with the largest overlap to the
charged pion is
@f{align}{
	O_{\pi^+}(n) = u(n) \gamma_5 \bar{d}(n).
@f}
The Operators in the euclidean correlator are then chosen to be a pion
creation operator @f$ O_{\pi^+}^\dagger(n) @f$ and a pion annihilation 
operator @f$ O_{\pi^+}(m) @f$ at two different lattice points $m$, $n$.
Inserting this into the path integral eq. {eq:correlatorlqcd}
the fermionic part is a Grassmann integration that can be carried out
analytically by performing a Wick-contraction. It will be denoted as the
fermionic expectation value @f$ \langle \rangle_F @f$. 

@f{align}{
	 \bcontraction[1ex]{\langle O_{\pi^+}(\vec{m}, m_4) O_{\pi^+}^\dagger(\vec{n}, n_4) \rangle_F = \quad ~ \
		- \langle}{\bar{d}}{(\vec{m}, m_4)\gamma_5 u(\vec{m}, m_4) \bar{u}(\vec{n}, n_4) 
		\gamma_5}{d}
	\bcontraction[2ex]{\langle O_{\pi^+}(\vec{m}, m_4) \O_{\pi^+}^\dagger(\vec{n}, n_4) \rangle_F =\quad ~ \
		- \langle \bar{d}(\vec{m}, m_4) \gamma_5}{u}{(\vec{m}, m_4)}{\bar{u}}
	\langle O_{\pi^+}(\vec{m}, m_4) O_{\pi^+}^\dagger(\vec{n}, n_4) \rangle_F =&&
		- \langle \bar{d}(\vec{m}, m_4) \gamma_5 u(\vec{m}, m_4)
		\bar{u}(\vec{n}, n_4)
		\gamma_5 d(\vec{n}, n_4) \rangle_F \nonumber \\
	 =&&
		tr(\gamma_5 D_u^{-1}(\vec{m}, m_4 | \vec{n}, n_4) \gamma_5
		 D_d^{-1}(\vec{n}, n_4 | \vec{m}, m_4)^\dagger ).
@f}

With eq. eq:quarkline as stochastical estimate for the propagator, the correlator
can be expressed in terms of the quantities introduced in section sec:laphs
The cyclicity of the trace and the @f$ \gamma_5 @f$-trick

@f{align}{
	D_d^{-1}(m | n) 
		&=& \gamma_5 D_u^{-1}(n | m) \gamma_5
@f}

were used. More details are given in bulava2009}

@f{align}{
	C_\pi(m_4, n_4)
		&=& tr \left( \gamma_5 \left(\sum_{\vec{m}, \vec{n}'} V_s(\vec{m}, m_4)^\dagger
					\Omega_u^{-1}(\vec{m}, \vec{n}; m_4, n_4) \left(V_s(\vec{n}, n_4) P^{(b)} \rho \right) \right) \right. \nonumber \\
		&&	\left( \sum_{\vec{y}} \left( V_s(\vec{y}, m_4) P^{(b)} \rho \right)^\dagger \left( V_s(\vec{y}, m_4) P^{(b)} \rho \right) \right) 
				\gamma_5 \nonumber \nonumber \\
		&&		\gamma_5 \left(\sum_{\vec{m}, \vec{n}} V_s(\vec{n}, m_4)^\dagger
					\Omega_u^{-1}(\vec{n}, \vec{m}; m_4, n_4) \left( V_s(\vec{m}, n_4) P^{(b)} \rho \right) \right)^\dagger \\
		&& \left. \left( \sum_{\vec{y}} V_s(\vec{y}, n_4)^\dagger V_s(\vec{y}, n_4) \right) 
					\gamma_5 \right) \nonumber .
@f}

This is the same expression one would also obtain from eq. (eq:correlatorquarklines).
The 	correlator has now been calculated for every time combination $m_4, n_4$. 
The values with the same time difference
@f$ (m_4 - n_4 + T) \mod T @f$ were added up to the so-called correlation function or correlator

@f{align}{
	C_\pi(t) 
		&=& \sum_{\stackrel{m_4, n_4}{t = (m_4 - n_4+T) \mod T}} C(m_4, n_4).
@f}

The resulting correlator consists of @f$ T @f$ data points and due to periodic boundary
conditions on the lattice it has to be symmetric in @f$ t @f$.

To obtain multiple scattering phases and to obtain @f$ a_\ell @f$ and @f$ r_\text{eff} @f$
from a fit to the effective range expansion, the computation
of @f$ \pi\pi @f$ scattering can be performed multiple times with different momenta. As a first test 
before simulating two pions
with momenta, a single pion was given momentum and the results were checked by
comparing them to the theoretical expectation given by the dispersion relation.	

Momenta on the lattice are introduced by a Fourier transformation to momentum space.
The correlators are given in cite{bulava2009} as for @f$ p^2 = 0 @f$ before. Including the 
momenta is a generalzation of eq. (eq:laphsC2nomom).
defining formula is eq. (eq:correlatorquarklines):

@f{align}{
		C_\pi(\vec{p}; m_4, n_4)
		&=& tr \left( \gamma_5 \left(\sum_{\vec{m}, \vec{n}'} V_s(\vec{m}, m_4)^\dagger
					\Omega_u^{-1}(\vec{m}, \vec{n}; m_4, n_4) \left(V_s(\vec{n}, n_4) P^{(b)} \rho \right) \right) \right. \nonumber \\
		&&	\left( \sum_{\vec{y}} \left( V_s(\vec{y}, m_4) P^{(b)} \rho \right)^\dagger e^{-i \vec{p} \vec{y}} 
				\left( V_s(\vec{y}, m_4) P^{(b)} \rho \right) \right) 
				\gamma_5 \nonumber \nonumber \\
		&&		\gamma_5 \left(\sum_{\vec{m}, \vec{n}} V_s(\vec{n}, m_4)^\dagger
					\Omega_u^{-1}(\vec{n}, \vec{m}; m_4, n_4) \left( V_s(\vec{m}, n_4) P^{(b)} \rho \right) \right)^\dagger \\
		&& \left. \left( \sum_{\vec{y}} V_s(\vec{y}, n_4)^\dagger e^{i \vec{p} \vec{y}} V_s(\vec{y}, n_4) \right) 
					\gamma_5 \right) \nonumber . .
@f}

Indeed, in the limit case @f$ \vec{p} = 0 @f$ eq. (eq:laphsC2) assumes the form given in eq. (eq:laphsC2nomom).	
At this point it is important
to note, that the @f$ \gamma_5 @f$-trick does not change the sequential arrangement of the last two terms. 
An intuitive idea to save computational effort is multipying the first two brackets to obtain a @f$ u @f$-quark propagator
and daggering the resulting expression to obtain a @f$ d @f$-quark propagator. For general momenta this is something
different than the expression derived in eq (eq:laphsC2). The reason for this is, that the cyclic property of the
trace has been used. The first two factors together are not to be mistaken with the propagator.

For @f$ \vec{p} = \vec{0} @f$ the sums over @f$ \vec{y} @f$ and @f$ \vec{y_0} @f$ assume the unit matrix as values and do not
contibute.

The momentum on th lattice is quantized: @f$ \vec{p} = \frac{2 \pi}{L} \vec{n} @f$.
@f$ \vec{n} \in \mathbb{Z}^3 @f$ denotes the a three-vector in terms of lattice units.
Via @f$ \vec{n} @f$, direction and absolute value of the desired momentum can be chosen.


The stochastic estimates can be brought into a form which allows saving them
to hard disk: The perambulator. This allows to reuse inverted propagators to
"build" multiple correlation functions. This code calculates correlation 
functions in exactly this manner. It featurus

- IO functions for perambulators, random vectors, eigenvectors
- An Infile to specify exactly the physical process in question
  - Different quark flavors
  - Flexible Dilution schemes
  - Every linearly independent Dirac structure
  - Different momenta
- 13 different diagrams for one- and two-meson operators pre-implemented

As a starting point the main function is in @ref contract.cpp

The issue of reusing certain parts for multiple diagrams achieved by calculating 
and saving quarklines. These hold particular combinations of perambulators and 
operators that appear in multiple diagrams.

 In the work group codes existed codes to calculate the perambulators and 
 eigensystems as introduced in ch. sec:laphs already before the
 beginning of this master thesis. 

For the sake of notation, the arguments of functions will
be stylized by parentheses and pointers to data structures
are notated with attached square brackets.

The perambulator   
@f{align}{
  \mathbf{P} =  V_s^\dagger \Omega^{-1} V_s P^{(b)} \rho
@f}
and the term
@f{align}{
  \mathbf{B}
    &=&  V_s^\dagger \exp{(i \vec{p} \vec{x})} V_s
@f}
appear repeatedly in the euclidean correlator eq. (eq:laphsC2nomom), (eq:laphsC2). 

The perambulator can be used as it was read in. The term @f$ \mathbf{B} @f$ is saved into
OperatorForMesons In the infile only the squared absolute value of the
momentum is specified. The routine 

The desired
momenta and displacements are specified in the function call. If both are 
trivial @f$ \mathbf{B} = \mathbb{1} @f$ and it will be explicitly set to unity.

For full dilution of the dirac-index @f$ \mathbf{B} @f$ becomes blockdiagonal 
in dirac space. Writing the euclidean correlator in matrix-vector notation
the product @f$ \mathbf{P} \cdot \mathbf{B} @f$ occurs. At this point the 
blockdiagonal structure can be exploited. Rather than a multiplication
in which twelve of the sixteen blocks of @f$ \mathbf{B} @f$ are zero, it suffices
to perform matrix multiplications for the four nonzero blocks. 

The code I implemented exploits this feature to save the effort of the
mulitplication @f$ D^{-1} \cdot \Gamma @f$. To illustrate the idea, I want to
perform an example calculation. 

Let w.l.o.g 
@f{align}{
  \mathbf{P} \cdot \mathbf{B} &=& \begin{pmatrix}  0 & 1 & 2 & 3 \\ 4 & 5 & 6 & 7 
      \\ 8 & 9 & 10 & 11 \\ 12 & 13 & 14 & 15 \end{pmatrix} \cdot 
      \begin{pmatrix} A & 0 & 0 & 0 \\ 0 & B & 0 & 0 \\ 0 & 0 & C & 0 \\
      0 & 0 & 0 & D \end{pmatrix} \\
    &=& \begin{pmatrix} 0A & 1B & 2C & 3D \\ 4A & 5B & 6C & 7D \\ 
      8A & 9B & 10C & 11D \\ 12A & 13B & 14C & 15D \end{pmatrix} 
@f}

Here, @f$ \{0, \dots, 15\} @f$ denote blocks of size number_of_eigen_vec
@f$ \cdot @f$ number_of_dilution_E and @f$ \{A, \dots, D\} @f$ blocks of size 
number_of_dilution_E @f$ \cdot @f$ number_of_eigen_vec. Performing
@f$ D^{-1} \cdot \Gamma @f$ culminates in reordering columns of eqn. 
(eq:ptimesb) and multiplying with scalars. For example:
@f{align}{
  \mathbf{P} \cdot \mathbf{B} \cdot \gamma_5 &=& \begin{pmatrix} 
      0A & 1B & -2C & -3D \\ 4A & 5B & -6C & -7D \\ 
      8A & 9B & -10C & -11D \\ 12A & 13B & -14C & -15D \end{pmatrix}
@f}

The contraction code supports positve and negative OS quarks. The @f$ u @f$-quark propagator is 
given by

@f{align}{
  D_u^{-1} 
    &=& V_s V_s^\dagger \Omega_u^{-1} V_s V_s^\dagger P^{(b)} \rho (P^{(b)} 
        \rho)^\dagger \\
    &=& V_s \mathbf{P} (P^{(b)} \rho)^\dagger
@f}

Due to the @f$ \gamma_5 @f$-trick 
the @f$ d @f$-quark propagator can be obtained by
@f{align}{
  D_d^{-1} 
    &=& \gamma_5 (D_u^{-1})\dagger V_s^\dagger \gamma_5 \\
    &=& \gamma_5 V_s (P^{(b)} \rho) \mathbf{P}^\dagger V_s^\dagger \gamma_5
@f}
without the need to perform additional inversions for the Dirac matrix.

@todo generalalize that for s and c quarks

For general momenta @f$ \mathbf{P} @f$ and @f$ \mathbf{B} @f$
do not commute. The multiplications required are performed in 

The multiplication with the @f$ \Gamma @f$-structure 
is performed seperately to reuse results as often as possible. 

@todo link to OperatorsForMesons here.
Due to the blockdiagonal form of @f$ \mathbf{B} @f$ multiplication with a gamma-matrix
is equivalent to reordering columns in dirac-space. Details on this reordering
as well as an example are presented in appendix sec:gamma.
The u-quark term is additionaly
multiplied with @f$ P^{(b)} \rho @f$. The reason is, that in the end the trace of
the whole product has to be taken for which only the diagonal elements need
to be calculated. For general gamma structure the random vectors complicate that slightly.
Rather than @f$ 4 N_t N_{\nu} \times 4 N_{\nu} @f$ the dimension of
the matrices becomes @f$ 4 N_t N_{\nu} \times 4 N_{\nu}^{dil} @f$
Thus, the trace can be taken with a factor @f$ \frac{N_{\nu}}{N_{\nu}^{dil}} @f$
less floating point operations at the cost of an additional factor of 

The baseline of the program flow is
1. Read infile
2. Determine the data needed to calculate the specified diagrams: The class 
    GlobalData processes them to vectors of certain structs for Quarklines, 
    Correlators and Operators.

    Internally all physical loops are replaced with auto-loops over these 
    structs. This allows to decouple the computation from the physical context 
    and keep the code very clean where the indices are very complicated.
3. Read all necessary data. The perambulators and the random vectors are saved 
    as attributes of LapH::Perambulator and LapH::RandomVector. 
4. Build parts that are reusable. The operators are constructed 
    in LapH::OperatorsForMesons. This is the only step where the eigenvectors 
    are needed, thus they are read on the fly and discarded afterwards. The 
    quarklines are constructed in LapH::Quarklines.
    @note The operators referred to above live solely in time and eigenvector 
          space. In order to save memory if multiple Dirac structures are to be 
          calculated at once it has been factored out and is performed when 
          building a quarkline or a correlator. They are in particular not 
          identical to the quantum field operators
5. Construct the correlators. This happens in LapH::Correlators. The class has
    a private member function corresponding to and building the diagrams 
    admissible in the infiles.
