/**
 * GECKO-A Web Interface - Main Application Module
 * Handles UI state management and user interactions
 *
 * Author: Deeksha Sharma
 */

const GeckoApp = {
    currentJobId: null,
    refreshInterval: null,

    /**
     * Initialize the application
     */
    init() {
        console.log('[GeckoApp] Initializing v3.0.0...');
        this.checkLibraries();
        this.loadCompounds();  // Load compounds from API
        this.loadJobs();
        this.startAutoRefresh();
        console.log('[GeckoApp] Initialization complete');
    },

    /**
     * Load compounds from API and populate dropdown
     * Uses the fast /api/compounds/dropdown endpoint
     */
    async loadCompounds() {
        const select = document.getElementById('voc_select');
        const countEl = document.getElementById('compound-count');

        try {
            // Use the fast dropdown endpoint
            const data = await GeckoAPI.getCompoundsForDropdown();

            // Clear loading message
            select.innerHTML = '<option value="" disabled selected>Select a Chemical...</option>';

            // Add compounds by category (simple structure: { category_name: [compounds] })
            if (data.categories) {
                for (const [categoryName, compounds] of Object.entries(data.categories)) {
                    if (compounds && compounds.length > 0) {
                        const optgroup = document.createElement('optgroup');
                        optgroup.label = categoryName;

                        compounds.forEach(compound => {
                            const option = document.createElement('option');
                            option.value = compound;
                            // Format display name (capitalize words)
                            const displayName = compound
                                .replace(/_/g, ' ')
                                .replace(/-/g, '-')
                                .split(' ')
                                .map(word => word.charAt(0).toUpperCase() + word.slice(1))
                                .join(' ');
                            option.textContent = displayName;
                            optgroup.appendChild(option);
                        });

                        select.appendChild(optgroup);
                    }
                }
            }

            // Add custom option at the end
            const customOption = document.createElement('option');
            customOption.value = 'custom';
            customOption.textContent = 'Other (Custom SMILES)...';
            select.appendChild(customOption);

            // Update count display
            if (countEl && data.total_compounds) {
                countEl.textContent = `${data.total_compounds} compounds available`;
            }

            console.log(`[GeckoApp] Loaded ${data.total_compounds || 0} compounds from ${Object.keys(data.categories || {}).length} categories`);

        } catch (e) {
            console.error('[GeckoApp] Failed to load compounds:', e);
            // Show error state
            select.innerHTML = '<option value="" disabled selected>Error loading compounds</option><option value="custom">Enter Custom SMILES...</option>';
            if (countEl) {
                countEl.textContent = 'Could not load compound list. Use custom SMILES input.';
                countEl.style.color = '#c00';
            }
        }
    },

    /**
     * Check library availability
     */
    checkLibraries() {
        console.log('[LibraryCheck] cytoscape:', typeof cytoscape !== 'undefined' ? 'loaded' : 'MISSING');
        console.log('[LibraryCheck] dagre:', typeof dagre !== 'undefined' ? 'loaded' : 'MISSING');
        console.log('[LibraryCheck] SmilesDrawer:', typeof SmilesDrawer !== 'undefined' ? 'loaded' : 'MISSING');

        if (typeof cytoscape !== 'undefined') {
            const testCy = cytoscape({ headless: true });
            try {
                testCy.layout({ name: 'dagre' });
                console.log('[LibraryCheck] dagre layout: available');
                testCy.destroy();
            } catch (e) {
                console.error('[LibraryCheck] dagre layout: NOT AVAILABLE -', e.message);
            }
        }
    },

    /**
     * Toggle custom VOC input visibility
     */
    toggleCustomVoc() {
        const select = document.getElementById('voc_select');
        const customInput = document.getElementById('voc_custom');
        customInput.style.display = select.value === 'custom' ? 'block' : 'none';
    },

    /**
     * Toggle advanced options panel
     */
    toggleAdvanced() {
        const adv = document.getElementById('advanced-options');
        const btn = document.querySelector('button[onclick="GeckoApp.toggleAdvanced()"]');
        if (adv.style.display === 'none') {
            adv.style.display = 'block';
            btn.innerText = '▼ Advanced Options (Generator & Box Model)';
        } else {
            adv.style.display = 'none';
            btn.innerText = '▶ Advanced Options (Generator & Box Model)';
        }
    },

    /**
     * Get scenario parameters from form
     * @returns {Object} Scenario parameters
     */
    getScenarioParams() {
        return {
            // Generator settings
            max_generations: parseInt(document.getElementById('max_generations').value),
            vapor_pressure_threshold: parseFloat(document.getElementById('vp_threshold').value),

            // Box model environment
            temperature_k: parseFloat(document.getElementById('temp_k').value),
            rh_percent: parseFloat(document.getElementById('rh_pct').value),
            latitude_degrees: parseFloat(document.getElementById('latitude').value),
            date_day: parseInt(document.getElementById('date_day').value),
            date_month: parseInt(document.getElementById('date_month').value),
            date_year: parseInt(document.getElementById('date_year').value),

            // Initial conditions
            initial_o3_ppb: parseFloat(document.getElementById('init_o3').value),
            initial_nox_ppb: parseFloat(document.getElementById('init_nox').value),
            seed_aerosol_ug_m3: parseFloat(document.getElementById('seed_aerosol').value),
            dilution_rate_s1: parseFloat(document.getElementById('dilution').value)
        };
    },

    /**
     * Get extended generator options from form
     * @returns {Object} Extended generator options
     */
    getExtendedGeneratorOptions() {
        return {
            vapor_pressure_method: document.getElementById('vp_method')?.value || 'nannoolal',
            critvp: parseFloat(document.getElementById('vp_threshold').value),
            max_generations: parseInt(document.getElementById('max_generations').value),
            nitrate_yield_method: document.getElementById('nitrate_method')?.value || 'carter',
            enable_oh_reactions: document.getElementById('enable_oh')?.checked ?? true,
            enable_o3_reactions: document.getElementById('enable_o3')?.checked ?? true,
            enable_no3_reactions: document.getElementById('enable_no3')?.checked ?? true,
            enable_photolysis: document.getElementById('enable_photolysis')?.checked ?? true,
            enable_isomerization: document.getElementById('enable_isomerization')?.checked ?? true,
            enable_ro2_permutations: document.getElementById('enable_ro2_permutations')?.checked ?? true,
            enable_pan_decomposition: document.getElementById('enable_pan_decomp')?.checked ?? true,
            enable_criegee_chemistry: document.getElementById('enable_criegee')?.checked ?? false,
            enable_autoxidation: document.getElementById('enable_autoxidation')?.checked ?? false
        };
    },

    /**
     * Get extended box model options from form
     * @returns {Object} Extended box model options
     */
    getExtendedBoxModelOptions() {
        return {
            simulation_hours: parseFloat(document.getElementById('sim_hours')?.value || 24.0),
            output_interval_minutes: parseFloat(document.getElementById('output_interval')?.value || 5.0),
            temperature_k: parseFloat(document.getElementById('temp_k').value),
            pressure_hpa: parseFloat(document.getElementById('pressure_hpa')?.value || 1013.25),
            rh_percent: parseFloat(document.getElementById('rh_pct').value),
            latitude_degrees: parseFloat(document.getElementById('latitude').value),
            longitude_degrees: parseFloat(document.getElementById('longitude')?.value || 0.0),
            date_day: parseInt(document.getElementById('date_day').value),
            date_month: parseInt(document.getElementById('date_month').value),
            date_year: parseInt(document.getElementById('date_year').value),
            local_hour_start: parseFloat(document.getElementById('start_hour')?.value || 6.0),
            initial_voc_ppb: parseFloat(document.getElementById('init_voc')?.value || 10.0),
            initial_o3_ppb: parseFloat(document.getElementById('init_o3').value),
            initial_no_ppb: parseFloat(document.getElementById('init_no')?.value || 5.0),
            initial_no2_ppb: parseFloat(document.getElementById('init_no2')?.value || 5.0),
            initial_hono_ppb: parseFloat(document.getElementById('init_hono')?.value || 0.1),
            initial_h2o2_ppb: parseFloat(document.getElementById('init_h2o2')?.value || 1.0),
            initial_hcho_ppb: parseFloat(document.getElementById('init_hcho')?.value || 2.0),
            initial_co_ppm: parseFloat(document.getElementById('init_co')?.value || 0.1),
            initial_ch4_ppm: parseFloat(document.getElementById('init_ch4')?.value || 1.8),
            seed_aerosol_ug_m3: parseFloat(document.getElementById('seed_aerosol').value),
            seed_aerosol_density_g_cm3: parseFloat(document.getElementById('seed_density')?.value || 1.4),
            organic_aerosol_activity_coefficient: parseFloat(document.getElementById('activity_coeff')?.value || 1.0),
            dilution_rate_s1: parseFloat(document.getElementById('dilution').value),
            deposition_velocity_cm_s: parseFloat(document.getElementById('dep_velocity')?.value || 0.0),
            enable_heterogeneous_reactions: document.getElementById('enable_heterogeneous')?.checked ?? true,
            enable_wall_loss: document.getElementById('enable_wall_loss')?.checked ?? false,
            photolysis_scaling_factor: parseFloat(document.getElementById('photolysis_scaling')?.value || 1.0)
        };
    },

    /**
     * Start a new job
     */
    async startJob() {
        const select = document.getElementById('voc_select');
        let voc = select.value;

        if (voc === 'custom') {
            voc = document.getElementById('voc_custom').value;
        }

        const jobType = document.getElementById('job_type').value;
        if (!voc) {
            alert("Please select or enter a VOC name");
            return;
        }

        const scenario = this.getScenarioParams();
        const generatorOptions = this.getExtendedGeneratorOptions();
        const boxmodelOptions = this.getExtendedBoxModelOptions();

        try {
            let data;

            if (jobType === 'comparison') {
                // Use the comparison workflow API
                return this.startComparisonJob();
            } else {
                // Use the standard job API for generator or boxmodel
                data = await GeckoAPI.createJob(voc, jobType, scenario, generatorOptions, boxmodelOptions);
            }

            this.loadJobs();
            this.viewJob(data.job_id);
        } catch (e) {
            alert("Error starting job: " + e);
        }
    },

    /**
     * Load and display all jobs
     */
    async loadJobs() {
        const list = document.getElementById('job-list');

        try {
            const jobs = await GeckoAPI.getJobs();

            list.innerHTML = '';

            if (jobs.length === 0) {
                list.innerHTML = '<p style="color:#777;">No jobs found. Start a new simulation to create a job.</p>';
                return;
            }

            jobs.reverse().forEach(job => {
                const div = document.createElement('div');
                div.className = 'job-item';
                div.innerHTML = `
                    <div class="job-info">
                        <span class="job-name">${job.voc_name}</span>
                        <span class="job-type">(${job.type})</span>
                        <div class="job-id">ID: ${job.id}</div>
                        ${job.archived ? '<div class="job-archived">✓ Archived</div>' : ''}
                    </div>
                    <div class="job-actions">
                        <span class="status-badge status-${job.status}">${job.status}</span>
                        <button onclick="GeckoApp.viewJob('${job.id}')" class="small">View</button>
                        <button onclick="GeckoApp.deleteJob('${job.id}')" class="small danger">Delete</button>
                    </div>
                `;
                list.appendChild(div);
            });
        } catch (e) {
            console.error("Error loading jobs:", e);
            list.innerHTML = `<p style="color:red">Error loading jobs: ${e.message}. Is the server running?</p>`;
        }
    },

    /**
     * View job details
     * @param {string} jobId - Job ID
     */
    async viewJob(jobId) {
        this.currentJobId = jobId;

        // Add loading state to button
        const buttons = document.querySelectorAll('button[onclick*="viewJob"]');
        buttons.forEach(btn => {
            if (btn.onclick && btn.onclick.toString().includes(jobId)) {
                btn.classList.add('loading');
                btn.disabled = true;
                btn.innerText = 'Loading...';
            }
        });

        document.getElementById('details-container').style.display = 'block';
        document.getElementById('details-title').innerText = `Job Details: ${jobId}`;

        const job = await GeckoAPI.getJob(jobId);

        // Remove loading state
        buttons.forEach(btn => {
            btn.classList.remove('loading');
            btn.disabled = false;
            btn.innerText = 'View';
        });

        const logsDiv = document.getElementById('logs');
        logsDiv.innerHTML = job.logs.join('<br>') || 'No logs yet...';
        logsDiv.scrollTop = logsDiv.scrollHeight;

        const resultsContainer = document.getElementById('results-container');
        if (job.status === 'completed') {
            resultsContainer.style.display = 'block';
            this.loadResults(jobId);
        } else {
            resultsContainer.style.display = 'none';
        }
    },

    /**
     * Archive a job
     * @param {string} jobId - Job ID
     */
    async archiveJob(jobId) {
        if (!jobId) return;
        if (!confirm("Are you sure you want to archive this session?")) return;

        try {
            const { ok, data } = await GeckoAPI.archiveJob(jobId);
            if (ok) {
                alert(`Session archived successfully!\nLocation: ${data.path}`);
                this.loadJobs();
            } else {
                alert(`Failed to archive: ${data.message}`);
            }
        } catch (e) {
            alert(`Error: ${e}`);
        }
    },

    /**
     * Delete a job
     * @param {string} jobId - Job ID
     */
    async deleteJob(jobId) {
        if (!jobId) return;
        if (!confirm(`Are you sure you want to delete this job?\n\nJob ID: ${jobId}\n\nThis will remove it from the list. Archived data will be preserved.`)) return;

        try {
            const { ok, data } = await GeckoAPI.deleteJob(jobId);
            if (ok) {
                alert('Job deleted successfully!');
                this.loadJobs();
                if (this.currentJobId === jobId) {
                    document.getElementById('details-container').style.display = 'none';
                    this.currentJobId = null;
                }
            } else {
                alert(`Failed to delete: ${data.message}`);
            }
        } catch (e) {
            alert(`Error: ${e}`);
        }
    },

    /**
     * Load and display job results
     * @param {string} jobId - Job ID
     */
    async loadResults(jobId) {
        try {
            const data = await GeckoAPI.getResults(jobId);

            // Mechanism Summary
            this.renderMechanismSummary(data.mechanism_summary);

            // Pathway Diagram
            this.renderPathwayDiagram(data.pathway_diagram);

            // Static Diagram (Legacy)
            this.renderStaticDiagram(data.static_diagram);

            // Mechanism Exports
            this.renderMechanismExports(data.mechanism_exports);

            // Reactions
            document.getElementById('reactions-content').innerText = data.reactions || "No reactions generated.";

            // SMILES Data
            this.renderSmilesData(data.smiles_data);

            // Hide structures section (rendered in tree nodes now)
            const structuresContainer = document.getElementById('structures-container');
            if (structuresContainer) {
                structuresContainer.innerHTML = '';
                document.getElementById('structures-section').style.display = 'none';
            }

            // Reaction Tree
            const treeContainer = document.getElementById('cy');
            if (data.reaction_tree && data.reaction_tree.nodes && data.reaction_tree.nodes.length > 0) {
                await ReactionTreeViz.render(treeContainer, data.reaction_tree);
            } else {
                treeContainer.innerHTML = '<p>No reaction tree data available.</p>';
            }

            // Plots
            this.renderPlots(data.plots);

            // CSV Data
            this.renderCsvData(data.csv_data, jobId);

            // NEW: Mass Balance Verification
            this.loadMassBalance(jobId);

            // NEW: Populate 3D structure dropdown
            if (data.smiles_data && data.smiles_data.length > 0) {
                this.populate3DSpeciesDropdown(data.smiles_data);
            }

            // NEW: Show report download links
            this.showReportDownloads(jobId, data.reports);

            // NEW: Show mechanism reduction section for completed jobs
            this.showReductionSection(jobId);

            // NEW: Show comparison results if this is a comparison job
            if (data.comparison_results) {
                this.renderComparisonResults(data.comparison_results);
            }

        } catch (e) {
            console.error("Error loading results:", e);
            const container = document.getElementById('results-container');
            if (!container.innerHTML.includes(e.message)) {
                container.innerHTML += `<p style="color:red">Error loading results: ${e.message}</p>`;
            }
        }
    },

    /**
     * Render VOC comparison results
     */
    renderComparisonResults(results) {
        const section = document.getElementById('comparison-results-section');
        if (!results) {
            section.style.display = 'none';
            return;
        }

        section.style.display = 'block';
        const chartsContainer = document.getElementById('comparison-charts');
        const tableContainer = document.getElementById('comparison-table');

        // Render comparison charts
        if (results.charts) {
            chartsContainer.innerHTML = results.charts.map(chart => `
                <div style="text-align: center;">
                    <h5>${chart.title}</h5>
                    <img src="${chart.url}" alt="${chart.title}" style="max-width: 100%; border: 1px solid #ddd; border-radius: 4px;">
                </div>
            `).join('');
        }

        // Render comparison table
        if (results.summary_table) {
            let tableHtml = '<table class="data-table"><thead><tr><th>VOC</th>';
            const metrics = Object.keys(results.summary_table[Object.keys(results.summary_table)[0]] || {});
            metrics.forEach(m => tableHtml += `<th>${m}</th>`);
            tableHtml += '</tr></thead><tbody>';

            for (const [voc, values] of Object.entries(results.summary_table)) {
                tableHtml += `<tr><td><strong>${voc}</strong></td>`;
                metrics.forEach(m => {
                    const val = values[m];
                    tableHtml += `<td>${typeof val === 'number' ? val.toFixed(2) : val}</td>`;
                });
                tableHtml += '</tr>';
            }
            tableHtml += '</tbody></table>';
            tableContainer.innerHTML = tableHtml;
        }
    },

    /**
     * Render mechanism summary section
     */
    renderMechanismSummary(summary) {
        const section = document.getElementById('mechanism-summary-section');
        if (!summary) {
            section.style.display = 'none';
            return;
        }

        section.style.display = 'block';
        document.getElementById('summary-root').textContent = summary.root_species || '';
        document.getElementById('summary-formula').textContent = summary.root_formula || '';
        document.getElementById('summary-species').textContent = summary.total_species || '';
        document.getElementById('summary-reactions').textContent = summary.total_reactions || '';

        const pathwaysDiv = document.getElementById('summary-pathways');
        if (summary.primary_pathways && summary.primary_pathways.length > 0) {
            const oxidantColors = {
                'HO': '#2196F3', 'OH': '#2196F3',
                'O3': '#4CAF50',
                'NO3': '#9C27B0',
                'NO': '#FF9800',
                'HO2': '#00BCD4',
                'RO2': '#795548'
            };

            let html = '<div style="display: flex; flex-wrap: wrap; gap: 10px;">';
            summary.primary_pathways.forEach(p => {
                const color = oxidantColors[p.oxidant] || '#607D8B';
                html += `
                    <div style="background: ${color}22; border: 1px solid ${color}; border-radius: 4px; padding: 8px; min-width: 150px;">
                        <strong style="color: ${color};">+${p.oxidant}</strong><br>
                        <span style="font-family: monospace; font-size: 0.85em;">${p.products.join(', ')}</span><br>
                        <small>Yield: ${(p.branching_ratio * 100).toFixed(1)}%</small>
                    </div>
                `;
            });
            html += '</div>';
            pathwaysDiv.innerHTML = html;
        }
    },

    /**
     * Render pathway diagram section
     */
    renderPathwayDiagram(diagram) {
        const section = document.getElementById('pathway-diagram-section');
        if (!diagram || (!diagram.png && !diagram.svg)) {
            section.style.display = 'none';
            return;
        }

        section.style.display = 'block';
        const img = document.getElementById('pathway-diagram-img');
        const pngLink = document.getElementById('download-pathway-png');
        const svgLink = document.getElementById('download-pathway-svg');

        if (diagram.png) {
            img.src = diagram.png;
            pngLink.href = diagram.png;
            pngLink.style.display = 'inline-block';
        } else {
            pngLink.style.display = 'none';
        }

        if (diagram.svg) {
            if (!diagram.png) img.src = diagram.svg;
            svgLink.href = diagram.svg;
            svgLink.style.display = 'inline-block';
        } else {
            svgLink.style.display = 'none';
        }
    },

    /**
     * Render static diagram section
     */
    renderStaticDiagram(diagram) {
        const section = document.getElementById('static-diagram-section');
        if (!diagram || !diagram.png) {
            section.style.display = 'none';
            return;
        }

        section.style.display = 'block';
        document.getElementById('static-diagram-img').src = diagram.png;
        document.getElementById('download-diagram-png').href = diagram.png;

        const svgLink = document.getElementById('download-diagram-svg');
        if (diagram.svg) {
            svgLink.href = diagram.svg;
            svgLink.style.display = 'inline-block';
        } else {
            svgLink.style.display = 'none';
        }
    },

    /**
     * Render mechanism exports section
     */
    renderMechanismExports(exports) {
        const section = document.getElementById('mechanism-exports-section');
        if (!exports || exports.length === 0) {
            section.style.display = 'none';
            return;
        }

        section.style.display = 'block';
        const container = document.getElementById('mechanism-exports');
        const formatColors = { 'KPP': '#28a745', 'MCM': '#17a2b8', 'FACSIMILE': '#6c757d' };

        container.innerHTML = exports.map(exp => {
            const color = formatColors[exp.format] || '#007bff';
            return `
                <a href="${exp.url}" target="_blank"
                   style="display: inline-block; padding: 10px 20px; background-color: ${color};
                          color: white; text-decoration: none; border-radius: 4px;">
                    ${exp.format} Format<br>
                    <small style="opacity: 0.8;">${exp.filename}</small>
                </a>
            `;
        }).join('');
    },

    /**
     * Render SMILES data table
     */
    renderSmilesData(smilesData) {
        const section = document.getElementById('smiles-section');
        const container = document.getElementById('smiles-container');

        if (!smilesData || smilesData.length === 0) {
            section.style.display = 'none';
            return;
        }

        section.style.display = 'block';
        let html = '<table class="data-table"><thead><tr><th>Code</th><th>SMILES</th><th>MW</th><th>Functional Groups</th></tr></thead><tbody>';
        smilesData.forEach(item => {
            const fgs = item.functional_groups ? item.functional_groups.join(', ') : '';
            html += `<tr>
                <td>${item.code}</td>
                <td style="font-family: monospace;">${item.smiles}</td>
                <td>${item.molecular_weight ? item.molecular_weight.toFixed(1) : ''}</td>
                <td><small>${fgs}</small></td>
            </tr>`;
        });
        html += '</tbody></table>';
        container.innerHTML = html;
    },

    /**
     * Render plots
     */
    renderPlots(plots) {
        const container = document.getElementById('plots-container');
        container.innerHTML = '';

        if (!plots || plots.length === 0) {
            container.innerHTML = '<p>No plots generated.</p>';
            return;
        }

        plots.forEach(plot => {
            const div = document.createElement('div');
            div.className = 'plot-item';
            div.innerHTML = `
                <h4>${plot.title}</h4>
                <img src="${plot.url}" class="plot-img" alt="${plot.title}">
            `;
            container.appendChild(div);
        });
    },

    /**
     * Render CSV data table
     */
    renderCsvData(csvData, jobId) {
        const tableBody = document.querySelector('#data-table tbody');
        const tableHead = document.querySelector('#data-table thead');
        const downloadBtn = document.getElementById('download-csv');

        tableBody.innerHTML = '';
        tableHead.innerHTML = '';

        if (!csvData || csvData.length === 0) {
            tableBody.innerHTML = '<tr><td colspan="5">No aerosol data available (Generator job or no output).</td></tr>';
            downloadBtn.style.display = 'none';
            return;
        }

        // Create headers
        const headers = Object.keys(csvData[0]);
        const trHead = document.createElement('tr');
        headers.forEach(h => {
            const th = document.createElement('th');
            th.innerText = h;
            trHead.appendChild(th);
        });
        tableHead.appendChild(trHead);

        // Create rows
        csvData.forEach(row => {
            const tr = document.createElement('tr');
            headers.forEach(h => {
                const td = document.createElement('td');
                let val = row[h];
                if (typeof val === 'number') {
                    val = val.toExponential(3);
                }
                td.innerText = val;
                tr.appendChild(td);
            });
            tableBody.appendChild(tr);
        });

        downloadBtn.style.display = 'inline-block';
        downloadBtn.href = `/data/output/${jobId}/aerosol_data.csv`;
    },

    /**
     * Start auto-refresh for job list
     */
    startAutoRefresh() {
        this.refreshInterval = setInterval(() => {
            this.loadJobs();
            if (this.currentJobId) {
                GeckoAPI.getJob(this.currentJobId).then(job => {
                    if (job.status === 'running' || job.status === 'queued') {
                        this.viewJob(this.currentJobId);
                    }
                });
            }
        }, 3000);
    },

    // =========================================================================
    // NEW FEATURES - VOC Comparison, Mechanism Reduction, 3D Structures, etc.
    // =========================================================================

    /**
     * Toggle job type specific options (comparison mode, etc.)
     */
    toggleJobTypeOptions() {
        const jobType = document.getElementById('job_type').value;
        const comparisonOptions = document.getElementById('comparison-options');
        const vocSelect = document.getElementById('voc_select');
        const vocLabel = vocSelect.previousElementSibling;

        if (jobType === 'comparison') {
            comparisonOptions.style.display = 'block';
            vocSelect.parentElement.style.display = 'none';
            // Populate comparison VOC list
            this.populateComparisonVOCs();
        } else {
            comparisonOptions.style.display = 'none';
            vocSelect.parentElement.style.display = 'block';
        }
    },

    /**
     * Populate comparison VOC multi-select with same compounds as main dropdown
     */
    async populateComparisonVOCs() {
        const select = document.getElementById('comparison_vocs');
        try {
            const data = await GeckoAPI.getCompoundsForDropdown();
            select.innerHTML = '';

            if (data.categories) {
                for (const [categoryName, compounds] of Object.entries(data.categories)) {
                    const optgroup = document.createElement('optgroup');
                    optgroup.label = categoryName;
                    compounds.forEach(compound => {
                        const option = document.createElement('option');
                        option.value = compound;
                        option.textContent = compound.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
                        optgroup.appendChild(option);
                    });
                    select.appendChild(optgroup);
                }
            }
        } catch (e) {
            console.error('Failed to populate comparison VOCs:', e);
        }
    },

    /**
     * Start a VOC comparison workflow
     */
    async startComparisonJob() {
        const select = document.getElementById('comparison_vocs');
        const selectedVOCs = Array.from(select.selectedOptions).map(opt => opt.value);

        if (selectedVOCs.length < 2) {
            alert('Please select at least 2 VOCs to compare');
            return;
        }
        if (selectedVOCs.length > 10) {
            alert('Please select no more than 10 VOCs');
            return;
        }

        const options = {
            compareMechanismSize: document.getElementById('compare_mechanism_size').checked,
            compareProductDistribution: document.getElementById('compare_products').checked,
            compareSoaYield: document.getElementById('compare_soa').checked,
            compareRadicalBudget: document.getElementById('compare_radicals').checked
        };

        try {
            const data = await GeckoAPI.createComparisonWorkflow(selectedVOCs, options);
            this.loadJobs();
            this.viewJob(data.job_id);
        } catch (e) {
            alert('Error starting comparison: ' + e.message);
        }
    },

    /**
     * Run mechanism reduction on current job
     */
    async reduceMechanism() {
        if (!this.currentJobId) {
            alert('No job selected');
            return;
        }

        const options = {
            method: document.getElementById('reduction_method').value,
            errorThreshold: parseFloat(document.getElementById('reduction_error').value),
            minSpecies: parseInt(document.getElementById('reduction_min_species').value),
            preserveRadicals: document.getElementById('preserve_radicals').checked,
            preserveSoaPrecursors: document.getElementById('preserve_soa').checked
        };

        try {
            const data = await GeckoAPI.createReductionWorkflow(this.currentJobId, options);

            // Show results
            const resultsDiv = document.getElementById('reduction-results');
            resultsDiv.style.display = 'block';

            document.getElementById('red-orig-species').textContent = data.original_species || '-';
            document.getElementById('red-new-species').textContent = data.reduced_species || '-';
            document.getElementById('red-orig-reactions').textContent = data.original_reactions || '-';
            document.getElementById('red-new-reactions').textContent = data.reduced_reactions || '-';
            document.getElementById('red-species-pct').textContent = data.species_reduction_percent ? data.species_reduction_percent.toFixed(1) + '%' : '-';
            document.getElementById('red-reactions-pct').textContent = data.reactions_reduction_percent ? data.reactions_reduction_percent.toFixed(1) + '%' : '-';

            if (data.reduced_mechanism_url) {
                document.getElementById('download-reduced-mechanism').href = data.reduced_mechanism_url;
            }
        } catch (e) {
            alert('Reduction failed: ' + e.message);
        }
    },

    /**
     * Regenerate pathway diagram with current options
     */
    async regenerateDiagram() {
        if (!this.currentJobId) return;

        const layout = document.getElementById('diagram_layout').value;
        const showBranching = document.getElementById('diagram_branching').checked;
        const colorByReaction = document.getElementById('diagram_colors').checked;
        const showLegend = document.getElementById('diagram_legend').checked;

        try {
            const result = await GeckoAPI.generateVisualization(this.currentJobId, {
                layout: layout,
                showBranching: showBranching,
                colorByReaction: colorByReaction,
                showLegend: showLegend,
                format: 'png'
            });

            if (result.diagram_url) {
                const url = result.diagram_url + '?t=' + Date.now();
                document.getElementById('pathway-diagram-img').src = url;
                
                // Update download links
                const pngLink = document.getElementById('download-pathway-png');
                if (pngLink) {
                    pngLink.href = url;
                    // Ensure the filename is correct when saving
                    pngLink.download = `pathway_${layout}.png`;
                }
            }
        } catch (e) {
            console.error('Failed to regenerate diagram:', e);
            alert('Failed to regenerate diagram: ' + e.message);
        }
    },

    /**
     * Load 3D structure for selected species
     */
    async load3DStructure() {
        const select = document.getElementById('structure3d-species');
        const species = select.value;
        if (!species) return;

        // Get display name from selected option
        const selectedOption = select.options[select.selectedIndex];
        const speciesCode = selectedOption.dataset.code || selectedOption.text;
        const speciesSmiles = selectedOption.dataset.smiles || species;
        const functionalGroups = selectedOption.dataset.functionalGroups || '';

        // Update name display with relevant chemical info
        const nameDisplay = document.getElementById('structure3d-c-name');
        if (nameDisplay) {
            // If it's the root VOC, the code usually matches the name.
            // If it's a generated species (e.g. 1U5000), "Code" is opaque.
            // We'll show: "Code: [Functional Groups]" or "Code (SMILES)"
            if (functionalGroups) {
                 nameDisplay.textContent = `Compound ${speciesCode} (${functionalGroups})`;
            } else if (speciesSmiles && speciesSmiles !== speciesCode) {
                 // Truncate long SMILES
                 const displaySmiles = speciesSmiles.length > 50 ? speciesSmiles.substring(0, 47) + '...' : speciesSmiles;
                 nameDisplay.textContent = `Compound ${speciesCode} (${displaySmiles})`;
            } else {
                 nameDisplay.textContent = `Compound ${speciesCode}`;
            }
        }

        const viewerContainer = document.getElementById('viewer-3d');
        if (!viewerContainer) {
            console.error('3D viewer container not found');
            return;
        }

        // Clear existing content and show loading message
        viewerContainer.innerHTML = '<p style="color: #666; padding: 20px;">Loading 3D structure...</p>';

        try {
            const response = await GeckoAPI.get3DStructure(species, 'mol');

            // Use 3Dmol.js if available
            if (typeof $3Dmol !== 'undefined') {
                // Clear the container completely
                viewerContainer.innerHTML = '';

                // Ensure the container has proper dimensions before creating viewer
                viewerContainer.style.position = 'relative';
                viewerContainer.style.width = '100%';
                viewerContainer.style.height = '400px';

                // Create the viewer with explicit configuration
                const viewerConfig = {
                    backgroundColor: '#f8f9fa', // Very light gray/white background
                    antialias: true,
                    cartoonQuality: 10
                };

                const view = $3Dmol.createViewer(viewerContainer, viewerConfig);

                // Get the MOL data from the response
                const molData = response.data || response.structure;

                if (!molData) {
                    throw new Error('No structure data received from server');
                }

                // Add the model
                view.addModel(molData, 'mol');
                
                // Set style with proper CPK colors
                // Carbon: Dark Gray (#909090 or #333333) vs Black. "Jmol" uses various grays.
                // We define specific colors to meet the "general consensus" requirement explicitly.
                view.setStyle({}, {
                    stick: { 
                        radius: 0.15, 
                        colorscheme: {
                            'C': '#1A1A1A', // Dark Gray/Black for Carbon (avoid pure black for shading)
                            'H': '#FFFFFF',
                            'O': '#FF0D0D',
                            'N': '#3050F8',
                            'S': '#FFFF30',
                            'Cl': '#1FF01F',
                            'F': '#90E050',
                            'Br': '#A62929',
                            'I': '#940094'
                        }
                    },
                    sphere: { 
                        scale: 0.25, 
                        colorscheme: {
                            'C': '#1A1A1A', // Dark Gray/Black for Carbon
                            'H': '#FFFFFF', // White for Hydrogen
                            'O': '#FF0D0D', // Red for Oxygen
                            'N': '#3050F8', // Blue for Nitrogen
                            'S': '#FFFF30', // Yellow for Sulfur
                            'Cl': '#1FF01F',
                            'F': '#90E050',
                            'Br': '#A62929',
                            'I': '#940094'
                        }
                    }
                });

                // Add atom labels (Element symbol)
                // We iterate through atoms to add labels
                const model = view.getModel();
                const atoms = model.selectedAtoms({});
                
                for (let i = 0; i < atoms.length; i++) {
                    const atom = atoms[i];
                    // Label atoms (usually all for small molecules, but skip H if too crowded)
                    if (atom.elem !== 'H' || atoms.length < 20) {
                        view.addLabel(atom.elem, {
                            position: {x: atom.x, y: atom.y, z: atom.z},
                            fontColor: 'black',
                            backgroundColor: 'white',
                            backgroundOpacity: 0.7,
                            borderColor: '#ccc', // Light gray border
                            borderThickness: 1,
                            fontSize: 14,
                            showBackground: true,
                            alignment: 'center'
                        });
                    }
                }

                // Center and render
                view.zoomTo();
                view.render();

                // Store for later use
                this.current3DViewer = view;
                this.current3DStructure = response;

                console.log('[GeckoApp] 3D structure loaded successfully for:', species);
            } else {
                viewerContainer.innerHTML = '<p style="color: #666; padding: 20px;">3Dmol.js not loaded. Structure data available for download.</p>';
                this.current3DStructure = response;
            }
        } catch (e) {
            console.error('[GeckoApp] Failed to load 3D structure:', e);
            viewerContainer.innerHTML = `<p style="color: #f00; padding: 20px;">Failed to load 3D structure: ${e.message}<br><small>Note: 3D generation requires RDKit on the server.</small></p>`;
        }
    },

    /**
     * Download 3D structure in specified format
     */
    async download3DStructure(format) {
        const species = document.getElementById('structure3d-species').value;
        if (!species) {
            alert('Please select a species first');
            return;
        }

        try {
            const response = await GeckoAPI.get3DStructure(species, format);
            const structureData = response.data || response.structure;
            const blob = new Blob([structureData], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `${species}.${format}`;
            a.click();
            URL.revokeObjectURL(url);
        } catch (e) {
            alert('Failed to download structure: ' + e.message + '\nNote: 3D generation requires RDKit on the server.');
        }
    },

    /**
     * Load mass balance data for current job
     */
    async loadMassBalance(jobId) {
        try {
            const data = await GeckoAPI.getMassBalance(jobId);
            const section = document.getElementById('mass-balance-section');
            section.style.display = 'block';

            document.getElementById('mb-status').textContent = data.status || 'Unknown';
            document.getElementById('mb-status').className = 'status-badge status-' + (data.status === 'passed' ? 'completed' : 'failed');
            document.getElementById('mb-carbon').textContent = data.carbon_balance ? (data.carbon_balance * 100).toFixed(2) + '%' : '-';
            document.getElementById('mb-hydrogen').textContent = data.hydrogen_balance ? (data.hydrogen_balance * 100).toFixed(2) + '%' : '-';
            document.getElementById('mb-oxygen').textContent = data.oxygen_balance ? (data.oxygen_balance * 100).toFixed(2) + '%' : '-';
            document.getElementById('mb-nitrogen').textContent = data.nitrogen_balance ? (data.nitrogen_balance * 100).toFixed(2) + '%' : '-';

            if (data.warnings && data.warnings.length > 0) {
                document.getElementById('mb-warnings').innerHTML = '<strong>Warnings:</strong><br>' + data.warnings.join('<br>');
            }
        } catch (e) {
            console.log('Mass balance not available:', e);
        }
    },

    /**
     * Populate 3D structure species dropdown from job results
     */
    populate3DSpeciesDropdown(species) {
        const select = document.getElementById('structure3d-species');
        select.innerHTML = '<option value="">Select a species...</option>';

        // Store species data for later use (code -> smiles mapping)
        this.speciesSmilesMap = {};

        if (species && species.length > 0) {
            species.forEach(sp => {
                const code = sp.code || sp;
                const smiles = sp.smiles || '';
                // Try to get a more descriptive name if available
                const functionalGroups = sp.functional_groups ? sp.functional_groups.join(', ') : '';
                
                const option = document.createElement('option');
                // Use SMILES as value if available (more reliable for 3D generation)
                option.value = smiles || code;
                
                // Construct a more descriptive label for the dropdown
                let label = code;
                if (functionalGroups) {
                    label += ` [${functionalGroups}]`;
                } else if (smiles && smiles.length < 20) {
                     label += ` (${smiles})`;
                }
                option.textContent = label;
                
                option.dataset.code = code;
                option.dataset.smiles = smiles;
                if(functionalGroups) option.dataset.functionalGroups = functionalGroups;
                
                select.appendChild(option);
                this.speciesSmilesMap[code] = smiles;
            });
            document.getElementById('structures-3d-section').style.display = 'block';
        }
    },

    /**
     * Show PDF report download links
     */
    showReportDownloads(jobId, reports) {
        const section = document.getElementById('pdf-report-section');
        if (!reports) {
            section.style.display = 'none';
            return;
        }

        section.style.display = 'block';

        if (reports.pdf_url) {
            document.getElementById('download-pdf-report').href = reports.pdf_url;
            document.getElementById('download-pdf-report').style.display = 'inline-block';
        }
        if (reports.mass_balance_url) {
            document.getElementById('download-mass-balance-report').href = reports.mass_balance_url;
            document.getElementById('download-mass-balance-report').style.display = 'inline-block';
        }
        if (reports.summary_url) {
            document.getElementById('download-mechanism-summary').href = reports.summary_url;
            document.getElementById('download-mechanism-summary').style.display = 'inline-block';
        }
    },

    /**
     * Show mechanism reduction section for completed jobs
     */
    showReductionSection(jobId) {
        const section = document.getElementById('reduction-section');
        section.style.display = 'block';
        document.getElementById('reduction-results').style.display = 'none';
    }
};

// Initialize on DOM ready
document.addEventListener('DOMContentLoaded', () => {
    GeckoApp.init();
});

// Global function aliases for HTML onclick handlers
function toggleCustomVoc() { GeckoApp.toggleCustomVoc(); }
function toggleAdvanced() { GeckoApp.toggleAdvanced(); }
function startJob() { GeckoApp.startJob(); }
function viewJob(id) { GeckoApp.viewJob(id); }
function deleteJob(id) { GeckoApp.deleteJob(id); }
function archiveJob(id) { GeckoApp.archiveJob(id); }

// Export for module usage
if (typeof module !== 'undefined' && module.exports) {
    module.exports = GeckoApp;
}
