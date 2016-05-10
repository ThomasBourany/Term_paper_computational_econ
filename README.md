# Term_paper_computational_econ
Term paper for the course of computational economics.

In this program, I implement a (very) simple growth model with exogenous fiscal policy, with a method called "shooting algorithm"

The model in question is a growth-model, non-stochastic, with an government. The government is purchasing exogenously a stream of goods, financing itself with distortive tax and lump-sum tax to equalize its budget constraint. The particularities of this economy are: inelastic labor supply, Cobb-Douglas production function, CRRA utility function, no saving but perfect Ricardian equivalence.  


Key equations of the model:

- A law of motion of capital determined by the market clearing condition
(11.6.1)
$ k_t+1 = f(k_t)+ (1-delta) k-t - g_t - c_t $

- The euler equation
(11.6.3)
$ u'(c_t) = u'(c_t+1) ((1+tau_ct)/(1+tauct+1)) [(1-tau_kt+1)(f'(k_t+1) - delta) +1]   $ 

- The steady state relation, called the Augmented Golden Rule
(11.6.7)
$ delta + rho/(1-tau_k) = f'(k_bar) $

- A set of 6 equations that compute the other endogenous equilibrium quantities 

From Ljungqvist & Sargent, the method to solve this model is the following : 

1. Solve (11.6.4) for the terminal steady-state k that is associated with the permanent policy vector z=(g, tau_c, tau_k) i.e., find the solution of (11.6.7). \
2. Select a large time index S >> T and guess an initial consumption rate c0 . Compute u′(c0) and solve (11.6.1) for k1 . \
3. For t = 0, use (11.6.3) to solve for u′(ct+1). Then invert u' and compute ct+1 . Use (11.6.1) to compute kt+2 . \
4. Iterate on step 3 to compute candidate values kˆt,t = 1,...,S. \
5. Compute kˆS −k. \
6. If kˆS > k, raise c0 and compute a new kˆt,t = 1,...,S. \
7 . I f kˆ S < k , lower c0 \
8. In this way, search for a value of c0 that makes kˆS ≈ k \
9. Compute the lump-sum taxes that satisfies the government budget constraint at equality.\


- This is only a draft, more documentation will be provided in the following days, about the following points: 

- Why the "shooting algorithm?"
- A description of the convergence of the algorithm (shown in the graphs 1, 2 and thoses accompanying the policy experiment) and its accuracy
- A description of the implementation of "policy experiments" concerning spending and taxes, and computation of the path and displayed in Impulse-Response function. 


- The code need to be improved:
 
- Tests to be added
- More policy experiments



 

