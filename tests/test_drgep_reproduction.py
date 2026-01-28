
import unittest
import sys

# Simplified extract to avoid copy-paste errors
def _reduce_drgep(species, reactions, target_species, error_threshold, min_species, preserve_radicals, preserve_soa_precursors):
    # Build species interaction graph
    graph = {s: {} for s in species}
    for rxn in reactions:
        # FIX: Directional Edges (Product depends on Reactant)
        reactants = rxn["reactants"]
        products = rxn["products"]
        for r in reactants:
            if r in graph:
                for p in products:
                    if p in graph and p != r:
                        # R leads to P (Causality)
                        graph[r][p] = graph[r].get(p, 0) + 1

    # Normalize interaction coefficients
    for s1, connections in graph.items():
        total = sum(connections.values()) if connections else 1
        for s2 in connections:
            graph[s1][s2] /= total
            
    # Essential species
    essential = set(target_species)
    
    # DRGEP path-based analysis
    memo = {}
    def drgep_coefficient(source, target, visited):
        state = (source, target, tuple(sorted(list(visited))))
        if state in memo: return memo[state] # optimization
        
        if source == target: return 1.0
        
        visited.add(source)
        max_coeff = 0.0
        
        for neighbor, weight in graph.get(source, {}).items():
            if neighbor not in visited:
                path_coeff = weight * drgep_coefficient(neighbor, target, visited.copy())
                max_coeff = max(max_coeff, path_coeff)
                
        visited.remove(source)
        return max_coeff

    # Calculate importance
    importance = {}
    for s in species:
        max_importance = 0.0
        for target in target_species:
            if target in species:
                # Need visited set! The original code used a set instance
                coeff = drgep_coefficient(s, target, set())
                max_importance = max(max_importance, coeff)
        importance[s] = max_importance

    print(f"Graph: {graph}")
    print(f"Importance: {importance}")
    
    kept_species = set()
    for s, imp in sorted(importance.items(), key=lambda x: x[1], reverse=True):
        if s in essential or imp >= error_threshold:
            kept_species.add(s)
            
    return sorted(list(kept_species)), []

class TestDRGEP(unittest.TestCase):
    def test_downstream_dependency(self):
        """
        Test that downstream products are NOT kept if they don't affect target.
        A -> B -> C. Target = A.
        B and C should be removed.
        """
        species = ['A', 'B', 'C']
        reactions = [
            {'reactants': ['A'], 'products': ['B']},
            {'reactants': ['B'], 'products': ['C']},
        ]
        kept, _ = _reduce_drgep(species, reactions, ['A'], 0.01, 0, False, False)
        print(f"Target A, Kept: {kept}")
        self.assertNotIn('B', kept)
        self.assertNotIn('C', kept)

    def test_linear_dependency(self):
        """
        Test a simple linear chain: A -> B -> C -> D
        """
        species = ['A', 'B', 'C', 'D']
        reactions = [
            {'reactants': ['A'], 'products': ['B']},
            {'reactants': ['B'], 'products': ['C']},
            {'reactants': ['C'], 'products': ['D']},
        ]
        kept, _ = _reduce_drgep(species, reactions, ['D'], 0.01, 0, False, False)
        print(f"Target D, Kept: {kept}")
        self.assertIn('A', kept)
        self.assertIn('B', kept)
        self.assertIn('C', kept)

if __name__ == '__main__':
    unittest.main()
