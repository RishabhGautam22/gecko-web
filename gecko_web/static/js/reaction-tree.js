/**
 * GECKO-A Web Interface - Reaction Tree Visualization Module
 * Handles Cytoscape.js visualization of reaction pathways
 *
 * Author: Deeksha Sharma
 */

const ReactionTreeViz = {
    instance: null,
    drawer: null,

    /**
     * Initialize SmilesDrawer for molecule rendering
     */
    initSmilesDrawer() {
        if (typeof SmilesDrawer !== 'undefined' && !this.drawer) {
            try {
                this.drawer = new SmilesDrawer.Drawer({
                    width: 300,
                    height: 200,
                    compactDrawing: false,
                    atomVisualization: 'default'
                });
                console.log('[ReactionTree] SmilesDrawer initialized');
            } catch (e) {
                console.warn('[ReactionTree] SmilesDrawer init failed:', e);
            }
        }
    },

    /**
     * Generate molecule image from SMILES
     * @param {string} smiles - SMILES string
     * @returns {Promise<string|null>} Data URL of molecule image or null
     */
    getMoleculeImage(smiles) {
        return new Promise(resolve => {
            if (!smiles || !this.drawer) {
                return resolve(null);
            }
            const canvas = document.createElement('canvas');
            canvas.width = 300;
            canvas.height = 200;
            try {
                SmilesDrawer.parse(smiles, (tree) => {
                    this.drawer.draw(tree, canvas, 'light', false);
                    resolve(canvas.toDataURL());
                }, (err) => {
                    resolve(null);
                });
            } catch (e) {
                resolve(null);
            }
        });
    },

    /**
     * Destroy existing Cytoscape instance
     */
    cleanup() {
        if (this.instance) {
            try {
                if (typeof this.instance.destroy === 'function') {
                    this.instance.destroy();
                    console.log('[ReactionTree] Previous instance destroyed');
                }
            } catch (err) {
                console.warn("[ReactionTree] Failed to destroy previous cy instance:", err);
            }
            this.instance = null;
        }
    },

    /**
     * Run layout on the graph
     * @param {string} layoutName - Layout algorithm name
     * @returns {Promise<void>}
     */
    runLayout(layoutName = 'dagre') {
        return new Promise((resolve, reject) => {
            if (!this.instance) {
                reject(new Error('No Cytoscape instance'));
                return;
            }

            let layoutOptions;

            if (layoutName === 'dagre' && typeof dagre === 'undefined') {
                console.warn('[ReactionTree] dagre not available, falling back to breadthfirst');
                layoutName = 'breadthfirst';
            }

            switch (layoutName) {
                case 'dagre':
                    layoutOptions = {
                        name: 'dagre',
                        rankDir: 'TB',
                        padding: 50,
                        spacingFactor: 1.2,
                        nodeSep: 50,
                        rankSep: 100,
                        animate: false
                    };
                    break;
                case 'breadthfirst':
                    layoutOptions = {
                        name: 'breadthfirst',
                        directed: true,
                        padding: 50,
                        spacingFactor: 1.5,
                        animate: false
                    };
                    break;
                case 'circle':
                    layoutOptions = {
                        name: 'circle',
                        padding: 50,
                        animate: false
                    };
                    break;
                case 'grid':
                    layoutOptions = {
                        name: 'grid',
                        padding: 50,
                        animate: false
                    };
                    break;
                default:
                    layoutOptions = {
                        name: 'breadthfirst',
                        directed: true,
                        padding: 50,
                        animate: false
                    };
            }

            console.log('[ReactionTree] Running layout:', layoutName);

            try {
                const layout = this.instance.layout(layoutOptions);
                layout.on('layoutstop', () => {
                    console.log('[ReactionTree] Layout complete:', layoutName);
                    resolve();
                });
                layout.run();
            } catch (layoutErr) {
                console.error('[ReactionTree] Layout error:', layoutErr);
                if (layoutName !== 'grid') {
                    console.log('[ReactionTree] Trying fallback to grid layout');
                    const fallback = this.instance.layout({ name: 'grid', padding: 50, animate: false });
                    fallback.run();
                }
                resolve();
            }
        });
    },

    /**
     * Add tooltips to nodes
     */
    addTooltips() {
        if (!this.instance) return;

        this.instance.on('mouseover', 'node', function(evt) {
            const node = evt.target;
            const nodeData = node.data();
            const tooltip = document.createElement('div');
            tooltip.id = 'cy-tooltip';
            tooltip.style.cssText = `
                position: fixed;
                background-color: #333;
                color: #fff;
                padding: 8px 12px;
                border-radius: 4px;
                font-size: 12px;
                z-index: 10000;
                pointer-events: none;
                max-width: 300px;
            `;
            tooltip.innerHTML = `<strong>${nodeData.label}</strong><br>${nodeData.smiles || ''}<br><em>${nodeData.raw_formula || ''}</em>`;
            document.body.appendChild(tooltip);

            const updateTooltipPosition = (e) => {
                tooltip.style.left = (e.clientX + 15) + 'px';
                tooltip.style.top = (e.clientY + 15) + 'px';
            };
            document.addEventListener('mousemove', updateTooltipPosition);
            node.data('tooltipMoveHandler', updateTooltipPosition);
        });

        this.instance.on('mouseout', 'node', function(evt) {
            const tooltip = document.getElementById('cy-tooltip');
            if (tooltip) tooltip.remove();
            const handler = evt.target.data('tooltipMoveHandler');
            if (handler) document.removeEventListener('mousemove', handler);
        });
    },

    /**
     * Render the reaction tree
     * @param {HTMLElement} container - Container element
     * @param {Object} treeData - Tree data with nodes and edges
     */
    async render(container, treeData) {
        console.log('[ReactionTree] Starting reaction tree rendering...');
        console.log('[ReactionTree] Data received:', {
            hasTree: !!treeData,
            nodeCount: treeData?.nodes?.length || 0,
            edgeCount: treeData?.edges?.length || 0
        });

        // Cleanup previous instance
        this.cleanup();
        container.innerHTML = '';

        // Cleanup previous controls and warnings
        const oldControls = document.getElementById('cy-controls');
        if (oldControls) oldControls.remove();
        const oldWarning = document.getElementById('cy-warning');
        if (oldWarning) oldWarning.remove();

        if (typeof cytoscape === 'undefined') {
            throw new Error('Cytoscape library not loaded');
        }

        if (!treeData || !treeData.nodes || treeData.nodes.length === 0) {
            console.log('[ReactionTree] No data available');
            container.innerHTML = '<p>No reaction tree data available.</p>';
            return;
        }

        // Warning for large graphs
        if (treeData.nodes.length > 150) {
            const warning = document.createElement('div');
            warning.id = 'cy-warning';
            warning.style.cssText = `
                padding: 10px;
                background-color: #fff3cd;
                border: 1px solid #ffc107;
                border-radius: 4px;
                margin-bottom: 10px;
            `;
            warning.innerHTML = `Warning: Large graph detected (${treeData.nodes.length} nodes). Showing first 200 nodes for performance.`;
            container.parentNode.insertBefore(warning, container);
        }

        // Add controls
        const controls = document.createElement('div');
        controls.id = 'cy-controls';
        controls.style.marginBottom = '10px';
        controls.innerHTML = `
            <span style="margin-right: 10px; font-weight: bold;">Nodes: ${Math.min(treeData.nodes.length, 200)}</span>
            <button onclick="ReactionTreeViz.fitView()" style="width: auto; padding: 5px 10px; font-size: 0.8em;">Fit View</button>
            <select onchange="ReactionTreeViz.changeLayout(this.value)" style="width: auto; padding: 5px; font-size: 0.8em;">
                <option value="dagre">Hierarchical (Dagre)</option>
                <option value="breadthfirst">Breadth First</option>
                <option value="circle">Circle</option>
                <option value="grid">Grid</option>
            </select>
        `;
        container.parentNode.insertBefore(controls, container);

        // Initialize SmilesDrawer
        this.initSmilesDrawer();

        // Limit nodes for performance
        const maxNodes = 200;
        const nodesToRender = treeData.nodes.slice(0, maxNodes);
        const nodeIds = new Set(nodesToRender.map(n => n.id));

        // Filter edges to only include those between rendered nodes
        const edgesToRender = treeData.edges.filter(e =>
            nodeIds.has(e.from) && nodeIds.has(e.to)
        );

        console.log('[ReactionTree] Rendering:', {
            nodes: nodesToRender.length,
            edges: edgesToRender.length
        });

        // Prepare nodes with images
        const processedNodes = await Promise.all(nodesToRender.map(async n => {
            let img = null;
            if (n.smiles) {
                img = await this.getMoleculeImage(n.smiles);
            }
            return {
                data: {
                    id: n.id,
                    smiles: n.smiles,
                    imageUrl: img,
                    label: n.code || n.id,
                    raw_formula: n.raw_formula,
                    title: n.title
                }
            };
        }));

        // Prepare edges
        const processedEdges = edgesToRender.map(e => ({
            data: {
                source: e.from,
                target: e.to,
                yield: e.yield,
                label: e.label || '',
                reaction: e.reaction || ''
            }
        }));

        // Ensure container has dimensions
        if (container.offsetWidth === 0 || container.offsetHeight === 0) {
            console.warn('[ReactionTree] Container has zero dimensions, forcing size');
            container.style.width = '100%';
            container.style.height = '500px';
            container.style.minHeight = '500px';
        }

        // Initialize Cytoscape
        this.instance = cytoscape({
            container: container,
            elements: {
                nodes: processedNodes,
                edges: processedEdges
            },
            style: [
                {
                    selector: 'node',
                    style: {
                        'background-color': '#e8f4f8',
                        'border-width': 2,
                        'border-color': '#4a90a4',
                        'label': 'data(label)',
                        'color': '#333',
                        'font-size': '9px',
                        'font-weight': 'bold',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'text-wrap': 'wrap',
                        'text-max-width': '120px',
                        'width': 130,
                        'height': 90,
                        'shape': 'roundrectangle',
                        'padding': '10px'
                    }
                },
                {
                    selector: 'node[imageUrl]',
                    style: {
                        'background-image': 'data(imageUrl)',
                        'background-fit': 'contain',
                        'background-opacity': 1,
                        'background-color': '#ffffff',
                        'text-valign': 'bottom',
                        'text-margin-y': 5
                    }
                },
                {
                    selector: 'edge',
                    style: {
                        'width': 2,
                        'line-color': '#888',
                        'target-arrow-color': '#666',
                        'target-arrow-shape': 'triangle',
                        'curve-style': 'bezier',
                        'arrow-scale': 1.0
                    }
                }
            ],
            layout: { name: 'preset' }
        });

        console.log('[ReactionTree] Cytoscape instance created:', {
            nodeCount: this.instance.nodes().length,
            edgeCount: this.instance.edges().length
        });

        // Make instance available globally for backward compatibility
        window.cyInstance = this.instance;

        // Run initial layout after a short delay
        setTimeout(async () => {
            try {
                if (this.instance) {
                    this.instance.resize();
                    await this.runLayout('dagre');
                    this.instance.fit();
                    this.instance.center();
                    this.addTooltips();
                    console.log('[ReactionTree] Initial layout and fit complete');
                }
            } catch (e) {
                console.error('[ReactionTree] Post-init error:', e);
            }
        }, 100);
    },

    /**
     * Fit view to show all elements
     */
    fitView() {
        if (this.instance) {
            this.instance.fit();
        }
    },

    /**
     * Change layout algorithm
     * @param {string} layoutName - Layout name
     */
    async changeLayout(layoutName) {
        try {
            await this.runLayout(layoutName);
            if (this.instance) {
                this.instance.fit();
                this.instance.center();
            }
        } catch (e) {
            console.error('[ReactionTree] Layout change failed:', e);
        }
    }
};

// Global alias for backward compatibility
window.changeLayout = function(layoutName) {
    ReactionTreeViz.changeLayout(layoutName);
};

// Export for module usage
if (typeof module !== 'undefined' && module.exports) {
    module.exports = ReactionTreeViz;
}
