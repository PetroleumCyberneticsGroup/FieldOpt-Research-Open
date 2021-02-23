# FieldOpt-Research-Open
FieldOpt [Open Research Version] is a C++ programming framework
for efficient prototyping and testing of optimization methodologies
for problems involving large-scale numerical simulations.

The FieldOpt framework serves as a multi-disciplinary knowledge
repository of optimization and petroleum domain expertise.
Technology development is based on efficient knowledge
retention and integration of project work and synergy
through mining of implemented iterative procedures and
domain parametrizations.
FieldOpt facilitates research and innovation through up-scaling of
prototype methodology to realistic cases, coupling, integration and
hybridization of optimization methodology and problem solutions,
and cross-application of existing methods to new domains.

[//]: # (
Open-source software framework for efficient prototyping of
optimization methodologies for complex and large-scale systems
)

## Target problems
### Petroleum Field Development
- [x] Well placement optimization [[1]](#Bellout2012JntWplcCntrl)
- [x] Production optimization
- [x] Optimization of inflow-control valve settings
- [x] Well completion optimization and model-update while drilling
- [ ] Minimization of C02 emissions

## Optimization methodologies

### Deterministic
- [x] Compass Search (CS)
- [x] Asynchronous Paralell Pattern Search (APPS)
- [x] Derivative-Free Trust-Region Algorithm (DFTR) [[2]](#Silva2020DfTrAlgWcntrlOpt)

### Stochastic / probabilistic
- [x] Genetic Algorithm (GA)
- [x] Particle Swarm Optimization (PSO)
- [x] Covariance Matrix Adaption Evolutionary Strategy (CMA-ES)
- [x] Bayesian Optimization (EGO)
- [x] Simultaneous Perturbation Stochastic Approximation (SPSA)

### Hybrid approaches
- [ ] mPSO
- [ ] APPS/PSO + data-driven meta-optimization
- [ ] Joint optimization using embedded reduced-order sub-routines

### Problem structure
- [ ] Multi-level joint optimization (concurrent, sequential, embedded)
- [ ] Automatic variable segregation for multi-level optimization
- [x] Variable scaling

### Objective terms
- [x] Weighted function, Net Present Value
- [x] Well cost
- [x] Augmented terms: Geology & geophysics-based (SWCT)

### Thirdparty solvers/libraries
- [x] SNOPT [[3]](#Gill2002SNOPTSIAMRev)
- [x] Ensemble based Reservoir Tool (ERT) [[4]](#EquinorERT)
- [ ] TensorFlow

## Functionalities
### Interfaces subsurface flow simulators
- [x] Schlumberger's E100/E300/IX
- [x] Open Porous Media Flow
- [x] Stanford's AD-GPRS
- [x] Pre-/Post-processing
  - [x] E300 adjoint-gradient read-in

### Well trajectory development
- [ ] Automatic well planner (AWP) [[5]](#Kristoffersen2020AWPGeoUncer)
- [x] State-of-the-art well connection transmissibility factor calculation [[6]](#ResInsightv2020.04.1)
- [x] Variable mapping onto multi-segmented well model (WELSEGS/COMPSEGS/WSEGVALV) [[7]](#SLB2012EclipseTD)

### Well placement constraint-handling
- [ ] Method of Alternating Projections (MAP) [[8]](#Bellout2018EffConstrHandlWplcOpt)
- [ ] Length, inter-well distance, user-defined convex-polytope reservoir-boundary

### Network/facility modeling
- [ ] Topside facility model for CO2 emission calculation

### Uncertainty-handling
- [x] Expected cost function evaluation over realization set
- [ ] Reduced random sampling strategy[[9]](#Jesmani2020RedRanSamStr)

### Parallelization
- [x] Algorithm-level parallelization of cost function
evaluations (simulations) through MPI runtime library
[[10]](#Baumann2020FieProFrmwrk)






## References

[//]: # (== 1 ==)
<a id="Bellout2012JntWplcCntrl">[1]</a>
Bellout, M.C.; Echeverria Ciaurri, D.; Durlofsky, L.J.; Foss, B.; Kleppe, J.
(2012).
Joint optimization of oil well placement and controls.
Computational Geosciences, 16(4), pp.1061-1079.
https://doi.org/10.1007/s10596-012-9303-5

[//]: # (== 2 ==)
<a id="Silva2020DfTrAlgWcntrlOpt">[2]</a>
Silva, T.L.; Bellout, M.C.; Giuliani, C.; Camponogara, E.; Pavlov, A.
(2020).
A Derivative-Free Trust-Region Algorithm for Well Control Optimization.
17th European Conference on the Mathematics of Oil
Recovery, 14th-17th September, Online Event.
https://doi.org/10.3997/2214-4609.202035086

[//]: # (== 3 ==)
<a id="Gill2002SNOPTSIAMRev">[3]</a>
Gill, P.E.; Murray, W.; Saunders, M.A.
(2005).
SNOPT: An SQP Algorithm for Large-Scale Constrained Optimization.
SIAM Review, 47(1), pp.99-131.
http://dx.doi.org/10.1137/S0036144504446096

[//]: # (== 4 ==)
<a id="EquinorERT">[4]</a>
Equinor.
(2021).
Ensemble based Reservoir Tool.
https://github.com/equinor/ert

[//]: # (== 5 ==)
<a id="Kristoffersen2020AWPGeoUncer">[5]</a>
Kristoffersen, B.S.; Silva, T.L.; Bellout, M.C.; Berg, C.F.
(2020).
An Automatic Well Planner for Efficient Well Placement
Optimization Under Geological Uncertainty.
17th European Conference on the Mathematics of Oil
Recovery, 14th-17th September, Online Event.
https://doi.org/10.3997/2214-4609.202035211

[//]: # (== 6 ==)
<a id="ResInsightv2020.04.1">[6]</a>
Ceetron Solutions AS; Equinor ASA.
(2020).
ResInsight.
http://resinsight.org

[//]: # (== 7 ==)
<a id="SLB2012EclipseTD">[7]</a>
Schlumberger AS.
(2012).
Eclipse technical description.
Chp.44: Multi-segment Wells. pp.683-703.
https://www.software.slb.com/products/eclipse/simulators

[//]: # (== 8 ==)
<a id="Bellout2018EffConstrHandlWplcOpt">[8]</a>
Bellout, M.C.; Volkov, O.
(2018).
Development of efficient constraint-handling approaches
for well placement optimization.
16th European Conference on the Mathematics of Oil
Recovery, 3rd-6th September, Barcelona, Spain.
https://doi.org/10.3997/2214-4609.201802247

[//]: # (== 9 ==)
<a id="Jesmani2020RedRanSamStr">[9]</a>
Jesmani, M.; Jafarpour, B.; Bellout, M.C.; Foss, B.
(2020).
A reduced random sampling strategy
for fast robust well placement optimization.
Journal of Petroleum Science and Engineering, 184, pp.106414.
https://doi.org/10.1016/j.petrol.2019.106414

[//]: # (== 10 ==)
<a id="Baumann2020FieProFrmwrk">[10]</a>
Baumann, E.J.M.; Dale, S.I.; Bellout, M.C.
(2020).
FieldOpt: A powerful and effective programming
framework tailored for field development optimization.
Computers & Geosciences, 135, pp.104379.
https://doi.org/10.1016/j.cageo.2019.104379


















