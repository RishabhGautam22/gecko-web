/**
 * GECKO-A Web Interface - API Module
 * Handles all API communications with the backend
 *
 * Extended with new endpoints for:
 * - Compound database (100+ compounds)
 * - Reaction kinetics
 * - Combined workflow
 * - Mass balance verification
 * - VOC comparison
 * - Mechanism reduction
 * - 3D structures
 * - Enhanced visualization
 *
 * Author: Deeksha Sharma
 */

const GeckoAPI = {
    // ==========================================================================
    // Core Job Management
    // ==========================================================================

    /**
     * Create a new simulation job
     * @param {string} vocName - VOC name or SMILES
     * @param {string} jobType - 'generator', 'boxmodel', 'combined', 'comparison'
     * @param {Object} scenario - Basic scenario options
     * @param {Object} extendedGenerator - Extended generator options
     * @param {Object} extendedBoxmodel - Extended box model options
     * @returns {Promise<Object>} Job creation response
     */
    async createJob(vocName, jobType, scenario = null, extendedGenerator = null, extendedBoxmodel = null) {
        const response = await fetch('/api/jobs', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                voc_name: vocName,
                job_type: jobType,
                scenario: scenario,
                extended_generator: extendedGenerator,
                extended_boxmodel: extendedBoxmodel
            })
        });
        return response.json();
    },

    /**
     * Get all jobs
     * @returns {Promise<Array>} List of jobs
     */
    async getJobs() {
        const response = await fetch('/api/jobs');
        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }
        return response.json();
    },

    /**
     * Get job details
     * @param {string} jobId - Job ID
     * @returns {Promise<Object>} Job details including logs
     */
    async getJob(jobId) {
        const response = await fetch(`/api/jobs/${jobId}`);
        return response.json();
    },

    /**
     * Get job results
     * @param {string} jobId - Job ID
     * @returns {Promise<Object>} Job results
     */
    async getResults(jobId) {
        const response = await fetch(`/api/jobs/${jobId}/results`);
        if (!response.ok) {
            throw new Error(`Failed to load results: ${response.status}`);
        }
        return response.json();
    },

    /**
     * Archive a completed job
     * @param {string} jobId - Job ID
     * @returns {Promise<Object>} Archive response
     */
    async archiveJob(jobId) {
        const response = await fetch(`/api/jobs/${jobId}/archive`, { method: 'POST' });
        return { ok: response.ok, data: await response.json() };
    },

    /**
     * Delete a job
     * @param {string} jobId - Job ID
     * @returns {Promise<Object>} Delete response
     */
    async deleteJob(jobId) {
        const response = await fetch(`/api/jobs/${jobId}`, { method: 'DELETE' });
        return { ok: response.ok, data: await response.json() };
    },

    // ==========================================================================
    // Compound Database API
    // ==========================================================================

    /**
     * Get all available compounds
     * @returns {Promise<Object>} Compounds list with count
     */
    async getCompounds() {
        const response = await fetch('/api/compounds');
        if (!response.ok) {
            throw new Error(`Failed to load compounds: ${response.status}`);
        }
        return response.json();
    },

    /**
     * Get compound categories
     * @returns {Promise<Object>} Categories hierarchy
     */
    async getCompoundCategories() {
        const response = await fetch('/api/compounds/categories');
        if (!response.ok) {
            throw new Error(`Failed to load categories: ${response.status}`);
        }
        return response.json();
    },

    /**
     * Get compounds for dropdown (fast, optimized endpoint)
     * @returns {Promise<Object>} Simple category->compounds structure
     */
    async getCompoundsForDropdown() {
        const response = await fetch('/api/compounds/dropdown');
        if (!response.ok) {
            throw new Error(`Failed to load compounds: ${response.status}`);
        }
        return response.json();
    },

    /**
     * Get specific compound details
     * @param {string} compoundName - Compound name
     * @returns {Promise<Object>} Compound details
     */
    async getCompound(compoundName) {
        const response = await fetch(`/api/compounds/${encodeURIComponent(compoundName)}`);
        if (!response.ok) {
            throw new Error(`Compound not found: ${compoundName}`);
        }
        return response.json();
    },

    /**
     * Search compounds
     * @param {string} query - Search query
     * @returns {Promise<Object>} Search results
     */
    async searchCompounds(query) {
        const response = await fetch(`/api/compounds/search/${encodeURIComponent(query)}`);
        return response.json();
    },

    /**
     * Get compounds by category
     * @param {string} category - Category name
     * @returns {Promise<Object>} Compounds in category
     */
    async getCompoundsByCategory(category) {
        const response = await fetch(`/api/compounds/category/${encodeURIComponent(category)}`);
        if (!response.ok) {
            throw new Error(`Category not found: ${category}`);
        }
        return response.json();
    },

    // ==========================================================================
    // Reaction Kinetics API
    // ==========================================================================

    /**
     * Get reaction kinetics for a compound
     * @param {string} compoundName - Compound name
     * @param {number} temperatureK - Temperature in Kelvin (default 298)
     * @returns {Promise<Object>} Rate constants
     */
    async getKinetics(compoundName, temperatureK = 298.0) {
        const response = await fetch(
            `/api/kinetics/${encodeURIComponent(compoundName)}?temperature_k=${temperatureK}`
        );
        if (!response.ok) {
            throw new Error(`No kinetics data for: ${compoundName}`);
        }
        return response.json();
    },

    /**
     * Calculate atmospheric lifetime
     * @param {string} compoundName - Compound name
     * @param {Object} conditions - Oxidant concentrations
     * @returns {Promise<Object>} Lifetime data
     */
    async getAtmosphericLifetime(compoundName, conditions = {}) {
        const params = new URLSearchParams({
            oh_concentration: conditions.oh || 1e6,
            o3_concentration: conditions.o3 || 7e11,
            no3_concentration: conditions.no3 || 5e8,
            temperature_k: conditions.temperature || 298.0
        });
        const response = await fetch(
            `/api/kinetics/${encodeURIComponent(compoundName)}/lifetime?${params}`
        );
        if (!response.ok) {
            throw new Error(`No kinetics data for: ${compoundName}`);
        }
        return response.json();
    },

    // ==========================================================================
    // Combined Workflow API
    // ==========================================================================

    /**
     * Create a combined Generator + Box Model workflow
     * @param {Object} request - Combined workflow request
     * @returns {Promise<Object>} Job creation response
     */
    async createCombinedWorkflow(request) {
        const response = await fetch('/api/workflow/combined', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(request)
        });
        if (!response.ok) {
            const error = await response.json();
            throw new Error(error.detail || 'Failed to create workflow');
        }
        return response.json();
    },

    // ==========================================================================
    // Mass Balance API
    // ==========================================================================

    /**
     * Get mass balance verification for a job
     * @param {string} jobId - Job ID
     * @returns {Promise<Object>} Mass balance results
     */
    async getMassBalance(jobId) {
        const response = await fetch(`/api/jobs/${jobId}/mass-balance`);
        if (!response.ok) {
            const error = await response.json();
            throw new Error(error.detail || 'Failed to get mass balance');
        }
        return response.json();
    },

    // ==========================================================================
    // VOC Comparison API
    // ==========================================================================

    /**
     * Create a VOC comparison workflow
     * @param {Array<string>} vocNames - List of VOC names to compare
     * @param {Object} options - Comparison options
     * @returns {Promise<Object>} Job creation response
     */
    async createComparisonWorkflow(vocNames, options = {}) {
        const request = {
            voc_names: vocNames,
            generator_options: options.generatorOptions || null,
            compare_mechanism_size: options.compareMechanismSize !== false,
            compare_product_distribution: options.compareProductDistribution !== false,
            compare_soa_yield: options.compareSoaYield !== false,
            compare_radical_budget: options.compareRadicalBudget !== false
        };

        const response = await fetch('/api/workflow/comparison', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(request)
        });
        if (!response.ok) {
            const error = await response.json();
            throw new Error(error.detail || 'Failed to create comparison');
        }
        return response.json();
    },

    // ==========================================================================
    // Mechanism Reduction API
    // ==========================================================================

    /**
     * Create a mechanism reduction workflow
     * @param {string} sourceJobId - Source job ID with full mechanism
     * @param {Object} options - Reduction options
     * @returns {Promise<Object>} Job creation response
     */
    async createReductionWorkflow(sourceJobId, options = {}) {
        const request = {
            job_id: sourceJobId,
            reduction_method: options.method || 'drgep',
            target_species: options.targetSpecies || [],
            error_threshold: options.errorThreshold || 0.1,
            min_species: options.minSpecies || 20,
            preserve_radicals: options.preserveRadicals !== false,
            preserve_soa_precursors: options.preserveSoaPrecursors !== false
        };

        const response = await fetch('/api/mechanism/reduce', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(request)
        });
        if (!response.ok) {
            const error = await response.json();
            throw new Error(error.detail || 'Failed to create reduction');
        }
        return response.json();
    },

    // ==========================================================================
    // 3D Structure API
    // ==========================================================================

    /**
     * Get 3D structure for a compound
     * @param {string} compoundName - Compound name
     * @param {string} format - Output format: 'mol', 'sdf', 'xyz', 'pdb'
     * @returns {Promise<Object>} 3D structure data
     */
    async get3DStructure(compoundName, format = 'mol') {
        const response = await fetch(
            `/api/structure/${encodeURIComponent(compoundName)}/3d?format=${format}`
        );
        if (!response.ok) {
            throw new Error(`No 3D structure for: ${compoundName}`);
        }
        return response.json();
    },

    /**
     * Get 3D structures for all species in a job
     * @param {string} jobId - Job ID
     * @param {number} maxStructures - Maximum number of structures
     * @returns {Promise<Object>} 3D structures data
     */
    async getJob3DStructures(jobId, maxStructures = 20) {
        const response = await fetch(
            `/api/jobs/${jobId}/structures/3d?max_structures=${maxStructures}`
        );
        if (!response.ok) {
            throw new Error(`Failed to get 3D structures for job: ${jobId}`);
        }
        return response.json();
    },

    // ==========================================================================
    // Enhanced Visualization API
    // ==========================================================================

    /**
     * Generate enhanced pathway visualization
     * @param {string} jobId - Job ID
     * @param {Object} options - Visualization options
     * @returns {Promise<Object>} Visualization result
     */
    async generateVisualization(jobId, options = {}) {
        const params = new URLSearchParams({
            layout: options.layout || 'hierarchical',
            show_branching: options.showBranching !== false,
            color_by_reaction: options.colorByReaction !== false,
            format: options.format || 'png'
        });

        const response = await fetch(`/api/jobs/${jobId}/visualize?${params}`, {
            method: 'POST'
        });
        if (!response.ok) {
            const error = await response.json();
            throw new Error(error.detail || 'Visualization failed');
        }
        return response.json();
    },

    // ==========================================================================
    // Utility Functions
    // ==========================================================================

    /**
     * Poll job status until completion
     * @param {string} jobId - Job ID
     * @param {Function} onProgress - Progress callback(job)
     * @param {number} interval - Poll interval in ms
     * @param {number} timeout - Timeout in ms
     * @returns {Promise<Object>} Completed job
     */
    async pollJobStatus(jobId, onProgress = null, interval = 2000, timeout = 3600000) {
        const startTime = Date.now();

        while (Date.now() - startTime < timeout) {
            const job = await this.getJob(jobId);

            if (onProgress) {
                onProgress(job);
            }

            if (job.status === 'completed' || job.status === 'failed') {
                return job;
            }

            await new Promise(resolve => setTimeout(resolve, interval));
        }

        throw new Error(`Job ${jobId} timed out after ${timeout}ms`);
    },

    /**
     * Get environment status
     * @returns {Promise<Object>} Environment check results
     */
    async checkEnvironment() {
        const response = await fetch('/api/environment');
        return response.json();
    }
};

// Export for module usage
if (typeof module !== 'undefined' && module.exports) {
    module.exports = GeckoAPI;
}
