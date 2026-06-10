# ternary-chemistry

**Chemical reaction networks where concentrations are {-1, 0, +1}. Catalysts, equilibrium, and mass conservation — in discrete algebra, not differential equations.**

In classical chemistry, you model reactions with systems of ordinary differential equations. Concentrations are real-valued, reaction rates are continuous, and equilibrium is a fixed point of the ODE. This crate asks: what happens when every concentration is clamped to one of three values — negative, neutral, or positive?

The answer turns out to be surprisingly rich. You get discrete reaction dynamics that always converge (the state space is finite), conservation laws that you can check in O(n), and catalysis that works by boosting reaction rates from 0 to ±1. The entire system runs in `no_std` with zero allocations beyond the initial `Vec` construction.

## The Insight

Continuous chemistry and ternary chemistry share the same skeleton: reactants combine, products form, rates govern speed, and the system evolves toward equilibrium. The difference is that ternary chemistry replaces the continuum with three discrete levels, and every operation becomes integer arithmetic clamped to [-1, 1].

This means:
- **No floating-point errors** — everything is exact integer math
- **Guaranteed termination** — finite state space means equilibrium is always reachable
- **Trivial conservation check** — sum concentrations before and after a step
- **`no_std` compatible** — no `f64`, no heap allocation during simulation

The tradeoff is resolution: you can't distinguish between "slightly positive" and "very positive." But for discrete decision systems (accept/reject/abstain, buy/hold/sell), ternary concentrations are the natural representation.

## Quick Start

```toml
[dependencies]
ternary-chemistry = "0.1.0"
```

```rust
use ternary_chemistry::*;

// Build a simple reaction network: A + B → C
let mut net = ReactionNetwork::new();
net.add_species(0, 1);   // A at concentration +1
net.add_species(1, 1);   // B at concentration +1
net.add_species(2, 0);   // C at concentration 0
net.add_reaction(vec![(0, 1)], vec![(2, 1)], 1);  // A → C at rate 1

// Apply the reaction
apply_reaction(&mut net, 0);
// A consumed, C produced (subject to ternary clamping)

// Run to equilibrium
let steps = equilibrium_concentrations(&mut net, 100);
println!("Equilibrium reached in {} steps", steps);
```

## Architecture

```
  ReactionNetwork
  ├── species: Vec<Species>
  │   └── Species { name: usize, concentration: i8 }  ∈ {-1, 0, +1}
  └── reactions: Vec<Reaction>
      └── Reaction { reactants, products, rate }

  Functions:
  ┌─────────────────────┐
  │ apply_reaction()    │  Single reaction step
  │ step()             │  All reactions, one pass
  │ equilibrium()      │  step() until stable or max_steps
  │ catalysis()        │  Boost a reaction rate
  │ conservation_check │  Sum concentrations before/after
  │ reaction_quotient  │  Q = Π(products) / Π(reactants)
  └─────────────────────┘
```

Species are identified by `usize` IDs. Concentrations are `i8` clamped to [-1, 1]. Reactions consume reactants and produce products, both modulated by the rate.

## API Reference

### Core Types

```rust
pub struct Species {
    pub name: usize,
    pub concentration: i8,  // clamped to [-1, 1]
}

pub struct Reaction {
    pub reactants: Vec<(usize, i8)>,  // (species_id, stoichiometry)
    pub products: Vec<(usize, i8)>,
    pub rate: i8,  // reaction rate ∈ [-1, 0, 1]
}

pub struct ReactionNetwork {
    pub species: Vec<Species>,
    pub reactions: Vec<Reaction>,
}
```

### Network Construction

```rust
ReactionNetwork::new() -> ReactionNetwork
net.add_species(name: usize, concentration: i8)
net.add_reaction(reactants: Vec<(usize, i8)>, products: Vec<(usize, i8)>, rate: i8)
```

### Simulation

```rust
apply_reaction(network: &mut ReactionNetwork, reaction_idx: usize)
step(network: &mut ReactionNetwork)  // apply all reactions once
equilibrium_concentrations(network: &mut ReactionNetwork, max_steps: usize) -> usize
```

- **`apply_reaction`** — checks reactant availability, consumes reactants, produces products. All concentrations clamped to [-1, 1]. If rate is 0 or reactants are insufficient, the reaction is skipped.
- **`step`** — applies all reactions in order, one pass.
- **`equilibrium_concentrations`** — calls `step` repeatedly until concentrations stop changing or `max_steps` is reached. Returns the number of steps taken.

### Analysis

```rust
catalysis(network: &ReactionNetwork, reaction_idx: usize, boost: i8) -> ReactionNetwork
conservation_check(network: &mut ReactionNetwork) -> (i16, bool)
reaction_quotient(network: &ReactionNetwork, reaction_idx: usize) -> i8
```

- **`catalysis`** — returns a new network with the specified reaction's rate boosted (clamped to [-1, 1]). The original network is not modified.
- **`conservation_check`** — steps the network once and compares total mass before/after. Returns `(total_mass_before, is_conserved)`. Ternary clamping may cause non-conservation.
- **`reaction_quotient`** — simplified Q = product of product concentrations / product of reactant concentrations, clamped to ternary. Returns 0 if any product concentration is 0, 1 if any reactant concentration is 0 (forward-favored).

## Real-World Example: Three-Way Decision Network

```rust
use ternary_chemistry::*;

// Model accept/reject/abstain as a reaction network
// Species: 0=accept_signal, 1=reject_signal, 2=neutral
let mut net = ReactionNetwork::new();
net.add_species(0, 1);   // accept signal present
net.add_species(1, -1);  // reject signal negative
net.add_species(2, 0);   // neutral starting state

// Accept signal catalyzes positive state
net.add_reaction(vec![(0, 1)], vec![(2, 1)], 1);

// Check if system reaches stable state
let steps = equilibrium_concentrations(&mut net, 50);
let final_state = net.species[2].concentration;
// +1 → accept, 0 → abstain, -1 → reject
```

## Reaction Semantics

A reaction fires when:
1. Its rate is nonzero
2. Every reactant has sufficient concentration (sign matches stoichiometry direction)

When it fires:
1. Each reactant's concentration decreases by `rate × sign(stoichiometry)`
2. Each product's concentration increases by `rate × stoichiometry`
3. All results are clamped to [-1, 1]

The clamping means that a product already at +1 can't go higher. This is the ternary analog of saturation — the system can't "accumulate" beyond its discrete levels.

## `no_std` Support

This crate is `#![no_std]` with `extern crate alloc`. The only allocation is the `Vec` used for species and reactions. Once constructed, `step` and `apply_reaction` don't allocate.

## Ecosystem

Part of the **ternary fleet** — 200+ crates, 4,300+ tests, one pattern:

- **ternary-core** — shared traits and Z₃ arithmetic
- **ternary-grid** — spatial grid with ternary cells
- **ternary-automata** — three-state cellular automata
- **ternary-compiler** — expression compiler and optimizer

## Open Questions

- **Stochastic reactions**: Current model is deterministic. Adding probabilistic firing (based on concentration levels) would model real chemical kinetics more faithfully.
- **Reversible reactions**: All reactions currently fire in one direction. Supporting equilibrium reactions (forward + backward rates) would enable true thermodynamic modeling.
- **Conservation under clamping**: Ternary clamping can break mass conservation. Characterizing *when* conservation holds (which stoichiometries, which rates) is an open problem.
- **Spatial chemistry**: Coupling with `ternary-grid` for reaction-diffusion systems.

## Stats

| Metric | Value |
|--------|-------|
| Lines of Rust | 335 |
| Tests | 15 |
| Public types | 3 |
| Public functions | 9 |
| `no_std` | Yes |
| `forbid(unsafe_code)` | Yes |

## License

MIT
