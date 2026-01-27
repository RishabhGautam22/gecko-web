import json

def inspect_edges():
    path = "/Users/rgpro/Documents/GEKO/data/output/boxmodel-1765036269/reaction_tree.json"
    with open(path) as f:
        data = json.load(f)
        
    print(f"Nodes: {len(data['nodes'])}")
    print(f"Edges: {len(data['edges'])}")
    if data['edges']:
        print(f"Sample edge: {data['edges'][0]}")

if __name__ == "__main__":
    inspect_edges()
