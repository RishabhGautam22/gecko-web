/**
 * GECKO-A Web Interface - 3D Structure Visualization Module
 *
 * Uses 3Dmol.js for interactive molecular visualization
 * Implements Recommendation #13: Add interactive 3D structures
 *
 * Author: Deeksha Sharma
 */

const Structure3D = {
    viewers: {},  // Cache of 3Dmol viewers by container ID

    /**
     * Initialize 3Dmol.js library if not loaded
     * @returns {Promise<void>}
     */
    async init() {
        if (window.$3Dmol) {
            return;
        }

        return new Promise((resolve, reject) => {
            const script = document.createElement('script');
            script.src = 'https://3dmol.org/build/3Dmol-min.js';
            script.async = true;
            script.onload = () => {
                console.log('3Dmol.js loaded successfully');
                resolve();
            };
            script.onerror = () => {
                reject(new Error('Failed to load 3Dmol.js'));
            };
            document.head.appendChild(script);
        });
    },

    /**
     * Create a 3D viewer in a container
     * @param {string|HTMLElement} container - Container element or ID
     * @param {Object} options - Viewer options
     * @returns {Object} 3Dmol viewer instance
     */
    createViewer(container, options = {}) {
        const element = typeof container === 'string'
            ? document.getElementById(container)
            : container;

        if (!element) {
            throw new Error(`Container not found: ${container}`);
        }

        const defaultOptions = {
            backgroundColor: 'white',
            antialias: true,
            cartoonQuality: 10
        };

        const viewerOptions = { ...defaultOptions, ...options };
        const viewer = $3Dmol.createViewer(element, viewerOptions);

        // Cache viewer
        const containerId = element.id || `viewer-${Date.now()}`;
        this.viewers[containerId] = viewer;

        return viewer;
    },

    /**
     * Load and display a molecule from MOL/SDF format
     * @param {Object} viewer - 3Dmol viewer instance
     * @param {string} molData - MOL/SDF format string
     * @param {Object} styleOptions - Visualization style options
     */
    loadMolecule(viewer, molData, styleOptions = {}) {
        viewer.removeAllModels();

        const model = viewer.addModel(molData, 'mol');

        // Default stick + ball style
        const defaultStyle = {
            stick: { radius: 0.15, colorscheme: 'Jmol' },
            sphere: { scale: 0.25, colorscheme: 'Jmol' }
        };

        const style = { ...defaultStyle, ...styleOptions };
        viewer.setStyle({}, style);

        viewer.zoomTo();
        viewer.render();

        return model;
    },

    /**
     * Load molecule from SMILES string (requires backend conversion)
     * @param {Object} viewer - 3Dmol viewer instance
     * @param {string} smiles - SMILES string
     * @param {string} compoundName - Compound name for API lookup
     * @returns {Promise<Object>} Model
     */
    async loadFromSMILES(viewer, smiles, compoundName) {
        try {
            const response = await fetch(`/api/structure/${encodeURIComponent(compoundName)}/3d?format=mol`);
            if (!response.ok) {
                throw new Error(`Failed to get 3D structure: ${response.status}`);
            }

            const data = await response.json();
            return this.loadMolecule(viewer, data.data);
        } catch (error) {
            console.error('Failed to load 3D structure:', error);
            throw error;
        }
    },

    /**
     * Set visualization style
     * @param {Object} viewer - 3Dmol viewer instance
     * @param {string} style - Style name: 'stick', 'sphere', 'line', 'cartoon'
     */
    setStyle(viewer, style) {
        const styles = {
            stick: { stick: { radius: 0.15, colorscheme: 'Jmol' } },
            sphere: { sphere: { scale: 0.3, colorscheme: 'Jmol' } },
            ballAndStick: {
                stick: { radius: 0.15, colorscheme: 'Jmol' },
                sphere: { scale: 0.25, colorscheme: 'Jmol' }
            },
            line: { line: { colorscheme: 'Jmol' } },
            surface: {
                stick: { radius: 0.1, colorscheme: 'Jmol' },
                surface: { opacity: 0.7, colorscheme: 'electrostatic' }
            }
        };

        const selectedStyle = styles[style] || styles.ballAndStick;
        viewer.setStyle({}, selectedStyle);
        viewer.render();
    },

    /**
     * Add surface representation
     * @param {Object} viewer - 3Dmol viewer instance
     * @param {string} surfaceType - 'VDW', 'MS', 'SAS', 'SES'
     * @param {Object} options - Surface options
     */
    addSurface(viewer, surfaceType = 'VDW', options = {}) {
        const defaultOptions = {
            opacity: 0.7,
            colorscheme: 'whiteCarbon'
        };

        const surfaceOptions = { ...defaultOptions, ...options };

        const surfaceTypes = {
            'VDW': $3Dmol.SurfaceType.VDW,
            'MS': $3Dmol.SurfaceType.MS,
            'SAS': $3Dmol.SurfaceType.SAS,
            'SES': $3Dmol.SurfaceType.SES
        };

        viewer.addSurface(surfaceTypes[surfaceType] || $3Dmol.SurfaceType.VDW, surfaceOptions);
        viewer.render();
    },

    /**
     * Remove all surfaces
     * @param {Object} viewer - 3Dmol viewer instance
     */
    removeSurfaces(viewer) {
        viewer.removeAllSurfaces();
        viewer.render();
    },

    /**
     * Enable/disable rotation animation
     * @param {Object} viewer - 3Dmol viewer instance
     * @param {boolean} enable - Enable or disable
     * @param {string} axis - Rotation axis: 'y', 'x', 'z'
     * @param {number} speed - Rotation speed (degrees per frame)
     */
    setRotation(viewer, enable, axis = 'y', speed = 1) {
        if (enable) {
            viewer.spin(axis, speed);
        } else {
            viewer.spin(false);
        }
    },

    /**
     * Take a screenshot
     * @param {Object} viewer - 3Dmol viewer instance
     * @param {string} filename - Output filename
     * @returns {string} Data URL of the image
     */
    screenshot(viewer, filename = 'molecule.png') {
        const canvas = viewer.pngURI();

        // Create download link
        const link = document.createElement('a');
        link.download = filename;
        link.href = canvas;
        link.click();

        return canvas;
    },

    /**
     * Create a 3D structure viewer component
     * @param {string} containerId - Container element ID
     * @param {string} compoundName - Compound name
     * @returns {Promise<Object>} Component with viewer and controls
     */
    async createStructureViewer(containerId, compoundName) {
        await this.init();

        const container = document.getElementById(containerId);
        if (!container) {
            throw new Error(`Container not found: ${containerId}`);
        }

        // Create viewer container and controls
        container.innerHTML = `
            <div class="structure-viewer-wrapper">
                <div id="${containerId}-viewer" class="structure-viewer-canvas" style="width: 100%; height: 400px; position: relative;"></div>
                <div class="structure-viewer-controls" style="padding: 10px; background: #f5f5f5; border-radius: 0 0 8px 8px;">
                    <div style="display: flex; gap: 10px; flex-wrap: wrap; align-items: center;">
                        <label>Style:</label>
                        <select id="${containerId}-style" style="padding: 5px;">
                            <option value="ballAndStick">Ball & Stick</option>
                            <option value="stick">Stick</option>
                            <option value="sphere">Space Fill</option>
                            <option value="line">Line</option>
                            <option value="surface">Surface</option>
                        </select>
                        <button id="${containerId}-rotate" style="padding: 5px 10px;">Toggle Rotation</button>
                        <button id="${containerId}-reset" style="padding: 5px 10px;">Reset View</button>
                        <button id="${containerId}-screenshot" style="padding: 5px 10px;">Screenshot</button>
                    </div>
                    <div style="margin-top: 8px; font-size: 12px; color: #666;">
                        <span>Drag to rotate | Scroll to zoom | Shift+drag to translate</span>
                    </div>
                </div>
            </div>
        `;

        // Create viewer
        const viewer = this.createViewer(`${containerId}-viewer`, {
            backgroundColor: 'white'
        });

        // Load molecule
        try {
            await this.loadFromSMILES(viewer, null, compoundName);
        } catch (error) {
            document.getElementById(`${containerId}-viewer`).innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #666;">
                    <p>Could not load 3D structure for ${compoundName}</p>
                </div>
            `;
            return null;
        }

        // Setup controls
        let rotating = false;

        document.getElementById(`${containerId}-style`).addEventListener('change', (e) => {
            this.setStyle(viewer, e.target.value);
        });

        document.getElementById(`${containerId}-rotate`).addEventListener('click', () => {
            rotating = !rotating;
            this.setRotation(viewer, rotating);
        });

        document.getElementById(`${containerId}-reset`).addEventListener('click', () => {
            viewer.zoomTo();
            viewer.render();
        });

        document.getElementById(`${containerId}-screenshot`).addEventListener('click', () => {
            this.screenshot(viewer, `${compoundName}_3d.png`);
        });

        return { viewer, container };
    },

    /**
     * Create a grid of 3D viewers for multiple compounds
     * @param {string} containerId - Container element ID
     * @param {Array} compounds - Array of compound objects with name and smiles
     * @param {number} columns - Number of columns in grid
     */
    async createStructureGrid(containerId, compounds, columns = 3) {
        await this.init();

        const container = document.getElementById(containerId);
        if (!container) {
            throw new Error(`Container not found: ${containerId}`);
        }

        const gridHTML = `
            <div class="structure-grid" style="display: grid; grid-template-columns: repeat(${columns}, 1fr); gap: 15px;">
                ${compounds.map((c, i) => `
                    <div class="structure-card" style="border: 1px solid #ddd; border-radius: 8px; overflow: hidden;">
                        <div id="${containerId}-mol-${i}" style="height: 250px; background: white;"></div>
                        <div style="padding: 10px; background: #f9f9f9; text-align: center;">
                            <strong>${c.label || c.name || c.id}</strong>
                            ${c.smiles ? `<br><code style="font-size: 10px; color: #666;">${c.smiles.substring(0, 40)}${c.smiles.length > 40 ? '...' : ''}</code>` : ''}
                        </div>
                    </div>
                `).join('')}
            </div>
        `;

        container.innerHTML = gridHTML;

        // Create viewers for each compound
        const viewers = [];
        for (let i = 0; i < compounds.length; i++) {
            const compound = compounds[i];
            const viewerContainer = document.getElementById(`${containerId}-mol-${i}`);

            try {
                const viewer = this.createViewer(viewerContainer, {
                    backgroundColor: 'white'
                });

                // Try to load from API first, then from mol_block if available
                if (compound.mol_block) {
                    this.loadMolecule(viewer, compound.mol_block);
                } else if (compound.name || compound.id) {
                    await this.loadFromSMILES(viewer, compound.smiles, compound.name || compound.id);
                }

                viewers.push(viewer);
            } catch (error) {
                console.warn(`Failed to create viewer for ${compound.name || compound.id}:`, error);
                viewerContainer.innerHTML = `<div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #999; font-size: 12px;">No 3D structure</div>`;
            }
        }

        return viewers;
    },

    /**
     * Load 3D structures for a completed job
     * @param {string} containerId - Container element ID
     * @param {string} jobId - Job ID
     * @param {number} maxStructures - Maximum number of structures to load
     */
    async loadJobStructures(containerId, jobId, maxStructures = 12) {
        try {
            const response = await fetch(`/api/jobs/${jobId}/structures/3d?max_structures=${maxStructures}`);
            if (!response.ok) {
                throw new Error(`Failed to load structures: ${response.status}`);
            }

            const data = await response.json();

            if (data.structures && data.structures.length > 0) {
                return await this.createStructureGrid(containerId, data.structures, 3);
            } else {
                document.getElementById(containerId).innerHTML = `
                    <p style="text-align: center; color: #666; padding: 20px;">
                        No 3D structures available for this job.
                    </p>
                `;
                return [];
            }
        } catch (error) {
            console.error('Failed to load job structures:', error);
            document.getElementById(containerId).innerHTML = `
                <p style="text-align: center; color: #c00; padding: 20px;">
                    Error loading 3D structures: ${error.message}
                </p>
            `;
            return [];
        }
    }
};

// Export for module usage
if (typeof module !== 'undefined' && module.exports) {
    module.exports = Structure3D;
}
