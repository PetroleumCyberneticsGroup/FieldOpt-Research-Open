# FieldOpt-Research-Open
FieldOpt [Open Research Version] is a C++ programming framework
for rapid prototyping and testing of optimization methodologies
for problems involving large-scale numerical simulations.

## Target problems
### Petroleum Field Development
- [x] Well placement optimization
- [x] Production optimization
- [x] Optimization of inflow-control valve settings

## Optimization methodologies

### Deterministic
- [x] CS
- [x] APPS

### Stochastic
- [ ] GA
- [ ] PSO

### Hybrid approaches
- [ ] APPS/PSO-ML1
- [ ] APPS/PSO-ML2

### Thirdparty solvers/libraries
- [ ] SNOPT
- [ ] TensorFlow

## Functionalities
### Reservoir simulator interfaces
- [x] Schlumberger's E100/E300/IX
- [x] Open Porous Media Flow
- [x] Stanford's AD-GPRS

- [ ] Pre-/Post-processing
  - [x] E300 adjoint-gradient read-in

### Parallelization
- [x] Algorithm-level parallelization of cost function
evaluations (simulations) through MPI runtime library
[[1]](#Baumann2020FieProFrmwrk)

### Uncertainty-handling
- [x] Expected cost function evaluation over realization set
- [ ] Reduced random sampling strategy
[[2]](#Jesmani2020RedRanSamStr)

### Well placement constraint-handling
- [ ] Method of Alternating Projections (MAP)
[[3]](#Bellout2018EffConstrHandlWplcOpt)
- [ ] Extended MAP

### Well trajectory development
- [ ] Automatic well planner (AWP)

## References
<a id="Baumann2020FieProFrmwrk">[1]</a>
Baumann, E.J.M.; Dale, S.I.; Bellout, M.C.
(2020).
FieldOpt: A powerful and effective programming
framework tailored for field development optimization.
Computers & Geosciences, 135, 104379.
https://doi.org/10.1016/j.cageo.2019.104379

<a id="Jesmani2020RedRanSamStr">[2]</a>
Jesmani, M.; Jafarpour, B.; Bellout, M.C.; Foss, B.
(2020).
A reduced random sampling strategy
for fast robust well placement optimization.
Journal of Petroleum Science and Engineering, 184, 106414.
https://doi.org/10.1016/j.petrol.2019.106414

<a id="Bellout2018EffConstrHandlWplcOpt">[3]</a>
Bellout, M.C.; Volkov, O.
(2018).
Development of efficient constraint-handling approaches
for well placement optimization.
16th European Conference on the Mathematics of Oil
Recovery, 3rd-6th September, Barcelona, Spain.
https://doi.org/10.3997/2214-4609.201802247
