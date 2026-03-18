# vertical-dynamics

Julia testbed for non-wetting impact models with:

- structured Dirichlet-to-Neumann operators on a 2D surface grid
- a prescribed-contact rigid-plane solver for non-axisymmetric bath dynamics
- dense reference and matrix-free Krylov solve paths

Start with [docs/project-map.md](/Users/eaguerov/Documents/Github/vertical-dynamics/docs/project-map.md).

Verification:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```
