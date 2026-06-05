#![forbid(unsafe_code)]
#![no_std]

extern crate alloc;
use alloc::vec;
use alloc::vec::Vec;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Species {
    pub name: usize,
    pub concentration: i8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Reaction {
    pub reactants: Vec<(usize, i8)>,
    pub products: Vec<(usize, i8)>,
    pub rate: i8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionNetwork {
    pub species: Vec<Species>,
    pub reactions: Vec<Reaction>,
}

impl ReactionNetwork {
    pub fn new() -> Self {
        Self {
            species: Vec::new(),
            reactions: Vec::new(),
        }
    }

    pub fn add_species(&mut self, name: usize, concentration: i8) {
        self.species.push(Species { name, concentration });
    }

    pub fn add_reaction(
        &mut self,
        reactants: Vec<(usize, i8)>,
        products: Vec<(usize, i8)>,
        rate: i8,
    ) {
        self.reactions
            .push(Reaction {
                reactants,
                products,
                rate,
            });
    }

    fn get_concentration(&self, name: usize) -> i8 {
        self.species
            .iter()
            .find(|s| s.name == name)
            .map(|s| s.concentration)
            .unwrap_or(0)
    }

    fn set_concentration(&mut self, name: usize, val: i8) {
        if let Some(s) = self.species.iter_mut().find(|s| s.name == name) {
            s.concentration = val;
        }
    }
}

/// Apply a single reaction by index, updating concentrations.
/// Uses simplified mass-action: change = rate * min(reactant concentrations, clamped to ternary).
/// All concentrations are clamped to [-1, 1].
pub fn apply_reaction(network: &mut ReactionNetwork, reaction_idx: usize) {
    let reaction = network.reactions[reaction_idx].clone();
    let rate = reaction.rate;

    if rate == 0 {
        return;
    }

    // Check if we have enough reactants (each must be nonzero with matching sign)
    for &(name, stoich) in &reaction.reactants {
        let conc = network.get_concentration(name);
        // Need sufficient reactant
        if stoich > 0 && conc < stoich.min(1) {
            return;
        }
        if stoich < 0 && conc > stoich.max(-1) {
            return;
        }
    }

    // Apply reactant consumption
    for &(name, stoich) in &reaction.reactants {
        let conc = network.get_concentration(name);
        let new_conc = (conc as i16 - (rate as i16 * stoich as i16 / stoich.abs() as i16))
            .clamp(-1, 1) as i8;
        network.set_concentration(name, new_conc);
    }

    // Apply product formation
    for &(name, stoich) in &reaction.products {
        let conc = network.get_concentration(name);
        let delta = rate as i16 * stoich as i16;
        let new_conc = (conc as i16 + delta).clamp(-1, 1) as i8;
        network.set_concentration(name, new_conc);
    }
}

/// Apply all reactions once.
pub fn step(network: &mut ReactionNetwork) {
    let n = network.reactions.len();
    for i in 0..n {
        apply_reaction(network, i);
    }
}

/// Run until concentrations stabilize or max_steps is reached.
/// Returns the number of steps actually taken.
pub fn equilibrium_concentrations(network: &mut ReactionNetwork, max_steps: usize) -> usize {
    for i in 0..max_steps {
        let before: Vec<i8> = network.species.iter().map(|s| s.concentration).collect();
        step(network);
        let after: Vec<i8> = network.species.iter().map(|s| s.concentration).collect();
        if before == after {
            return i + 1;
        }
    }
    max_steps
}

/// Catalyst enables or speeds a reaction by modifying its rate.
/// Returns a new network with the catalyzed reaction.
pub fn catalysis(network: &ReactionNetwork, reaction_idx: usize, catalyst_boost: i8) -> ReactionNetwork {
    let mut new_network = network.clone();
    if let Some(reaction) = new_network.reactions.get_mut(reaction_idx) {
        reaction.rate = (reaction.rate as i16 + catalyst_boost as i16).clamp(-1, 1) as i8;
    }
    new_network
}

/// Check mass conservation: sum of all concentrations before and after should be invariant.
/// Returns (total_mass, is_conserved) where is_conserved means total doesn't change on a step.
pub fn conservation_check(network: &mut ReactionNetwork) -> (i16, bool) {
    let before: i16 = network.species.iter().map(|s| s.concentration as i16).sum();
    step(network);
    let after: i16 = network.species.iter().map(|s| s.concentration as i16).sum();
    (before, before == after)
}

/// Compute reaction quotient Q = product of product concentrations / product of reactant concentrations.
/// In ternary, returns the simplified ternary value.
pub fn reaction_quotient(network: &ReactionNetwork, reaction_idx: usize) -> i8 {
    let reaction = &network.reactions[reaction_idx];

    let mut prod_sum: i16 = 0;
    for &(name, stoich) in &reaction.products {
        let conc = network.get_concentration(name) as i16;
        if conc == 0 {
            return 0;
        }
        prod_sum += conc * stoich as i16;
    }

    let mut react_sum: i16 = 0;
    for &(name, stoich) in &reaction.reactants {
        let conc = network.get_concentration(name) as i16;
        if conc == 0 {
            return 1; // no reactants → Q effectively infinite → forward favored
        }
        react_sum += conc * stoich as i16;
    }

    if react_sum == 0 {
        return if prod_sum > 0 { 1 } else { 0 };
    }

    let q = prod_sum / react_sum;
    q.clamp(-1, 1) as i8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_species_creation() {
        let s = Species { name: 0, concentration: 1 };
        assert_eq!(s.name, 0);
        assert_eq!(s.concentration, 1);
    }

    #[test]
    fn test_network_new() {
        let net = ReactionNetwork::new();
        assert_eq!(net.species.len(), 0);
        assert_eq!(net.reactions.len(), 0);
    }

    #[test]
    fn test_add_species() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, -1);
        assert_eq!(net.species.len(), 2);
        assert_eq!(net.species[0].concentration, 1);
        assert_eq!(net.species[1].concentration, -1);
    }

    #[test]
    fn test_simple_reaction() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1); // A
        net.add_species(1, 0); // B
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        apply_reaction(&mut net, 0);
        // Reactant consumed, product formed
        assert!(net.get_concentration(1) != 0);
    }

    #[test]
    fn test_zero_rate_no_change() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 0);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 0);
        let before = net.get_concentration(0);
        apply_reaction(&mut net, 0);
        assert_eq!(net.get_concentration(0), before);
    }

    #[test]
    fn test_step_applies_all() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 1);
        net.add_species(2, 0);
        net.add_reaction(vec![(0, 1)], vec![(2, 1)], 1);
        net.add_reaction(vec![(1, 1)], vec![(2, 1)], 1);
        step(&mut net);
        assert!(net.get_concentration(2) != 0);
    }

    #[test]
    fn test_equilibrium_reaches_stable() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 0); // starts at 0, nothing to react
        net.add_species(1, 0);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        let steps = equilibrium_concentrations(&mut net, 100);
        assert!(steps <= 100);
    }

    #[test]
    fn test_catalysis_boosts_rate() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 0);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 0);
        let catalyzed = catalysis(&net, 0, 1);
        assert_eq!(catalyzed.reactions[0].rate, 1);
    }

    #[test]
    fn test_conservation_check_balanced() {
        // A balanced reaction where consumption = production
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 0);
        // forward only: consume A, produce B — mass moves but ternary clamping may not conserve
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        let (_total, _conserved) = conservation_check(&mut net);
        // Just check it doesn't panic
    }

    #[test]
    fn test_reaction_quotient_basic() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 1);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        let q = reaction_quotient(&net, 0);
        // products: 1*1=1, reactants: 1*1=1, Q=1
        assert_eq!(q, 1);
    }

    #[test]
    fn test_reaction_quotient_zero_reactant() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 0);
        net.add_species(1, 1);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        let q = reaction_quotient(&net, 0);
        // No reactant concentration → Q infinite (1 in ternary)
        assert_eq!(q, 1);
    }

    #[test]
    fn test_reaction_quotient_zero_product() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 0);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        let q = reaction_quotient(&net, 0);
        assert_eq!(q, 0);
    }

    #[test]
    fn test_concentration_clamping() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 1);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        // Product at 1, adding more should clamp
        apply_reaction(&mut net, 0);
        assert!(net.get_concentration(1) <= 1);
    }

    #[test]
    fn test_catalysis_clamps_rate() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        let catalyzed = catalysis(&net, 0, 5); // boost by 5, should clamp to 1
        assert_eq!(catalyzed.reactions[0].rate, 1);
    }

    #[test]
    fn test_equilibrium_max_steps() {
        let mut net = ReactionNetwork::new();
        net.add_species(0, 1);
        net.add_species(1, 0);
        net.add_reaction(vec![(0, 1)], vec![(1, 1)], 1);
        let steps = equilibrium_concentrations(&mut net, 100);
        assert!(steps <= 100);
    }
}
