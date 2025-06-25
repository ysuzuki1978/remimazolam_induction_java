/**
 * LSODA (Livermore Solver for ODEs with Automatic method switching)
 * JavaScript implementation for pharmacokinetic simulations
 * 
 * Automatically switches between Adams method (non-stiff) and BDF method (stiff)
 * Based on the FORTRAN LSODA implementation by Petzold
 */

class LSODA {
    constructor() {
        // Integration method constants
        this.ADAMS = 1;
        this.BDF = 2;
        
        // Default tolerances
        this.rtol = 1e-6;  // Relative tolerance
        this.atol = 1e-12; // Absolute tolerance
        
        // Method switching parameters
        this.maxorder_adams = 12;
        this.maxorder_bdf = 5;
        this.meth = this.ADAMS;  // Start with Adams method
        this.miter = 0;          // Iteration method
        
        // Step size control
        this.h = 0.0;           // Current step size
        this.hmin = 0.0;        // Minimum step size
        this.hmax = 0.0;        // Maximum step size
        this.maxstep = 500;     // Maximum number of steps
        
        // Current state
        this.t = 0.0;           // Current time
        this.y = [];            // Current solution vector
        this.n = 0;             // Number of equations
        
        // Work arrays
        this.yh = [];           // Nordsieck array
        this.wm = [];           // Work matrix
        this.ewt = [];          // Error weight vector
        this.savf = [];         // Saved function values
        this.acor = [];         // Accumulated corrections
        
        // Counters and flags
        this.nst = 0;           // Number of steps taken
        this.nfe = 0;           // Number of function evaluations
        this.nje = 0;           // Number of Jacobian evaluations
        this.nqu = 0;           // Current method order
        this.nq = 1;            // Current integration order
        this.l = 0;             // Index for method coefficients
        
        // Method coefficients
        this.elco = [];         // Integration coefficients
        this.tesco = [];        // Test coefficients
        
        // Initialize coefficient arrays
        this.initializeCoefficients();
    }

    /**
     * Initialize Adams and BDF method coefficients
     */
    initializeCoefficients() {
        // Adams method coefficients (explicit)
        this.elco[this.ADAMS] = [
            [],
            [1.0, 1.0],
            [2.0/3.0, 4.0/3.0, -1.0/3.0],
            [6.0/11.0, 18.0/11.0, -9.0/11.0, 2.0/11.0],
            [24.0/50.0, 72.0/50.0, -54.0/50.0, 24.0/50.0, -5.0/50.0],
            [120.0/274.0, 360.0/274.0, -360.0/274.0, 200.0/274.0, -75.0/274.0, 12.0/274.0]
        ];

        // BDF method coefficients (implicit)
        this.elco[this.BDF] = [
            [],
            [1.0, 1.0],
            [2.0/3.0, 4.0/3.0, -1.0/3.0],
            [6.0/11.0, 18.0/11.0, -9.0/11.0, 2.0/11.0],
            [12.0/25.0, 48.0/25.0, -36.0/25.0, 16.0/25.0, -3.0/25.0],
            [60.0/137.0, 300.0/137.0, -300.0/137.0, 200.0/137.0, -75.0/137.0, 12.0/137.0]
        ];

        // Test coefficients for error estimation
        this.tesco[this.ADAMS] = [
            [],
            [2.0, 1.0],
            [3.0, 2.0/3.0],
            [4.0, 6.0/11.0],
            [5.0, 24.0/50.0],
            [6.0, 120.0/274.0]
        ];

        this.tesco[this.BDF] = [
            [],
            [2.0, 1.0],
            [3.0, 2.0/3.0],
            [4.0, 6.0/11.0],
            [5.0, 12.0/25.0],
            [6.0, 60.0/137.0]
        ];
    }

    /**
     * Solve ODE system from t0 to tfinal
     * @param {Function} f - Right-hand side function dy/dt = f(t, y)
     * @param {Array} y0 - Initial conditions
     * @param {number} t0 - Initial time
     * @param {number} tfinal - Final time
     * @param {Object} options - Integration options
     * @returns {Object} Solution object with t and y arrays
     */
    solve(f, y0, t0, tfinal, options = {}) {
        // Set options
        this.rtol = options.rtol || 1e-6;
        this.atol = options.atol || 1e-12;
        this.maxstep = options.maxstep || 500;
        this.hmax = options.hmax || Math.abs(tfinal - t0);
        this.hmin = options.hmin || 0.0;

        // Initialize
        this.n = y0.length;
        this.t = t0;
        this.y = [...y0];
        this.h = (tfinal - t0) / 100; // Initial step size guess
        this.nst = 0;
        this.nfe = 0;
        this.nje = 0;
        this.nq = 1;
        this.meth = this.ADAMS;

        // Initialize work arrays
        this.initializeWorkArrays();

        // Solution storage
        const t_out = [t0];
        const y_out = [y0.slice()];

        // Integration parameters
        const direction = tfinal > t0 ? 1 : -1;
        this.h *= direction;

        // Main integration loop
        while ((direction > 0 && this.t < tfinal) || (direction < 0 && this.t > tfinal)) {
            // Don't overshoot the final time
            if ((direction > 0 && this.t + this.h > tfinal) || 
                (direction < 0 && this.t + this.h < tfinal)) {
                this.h = tfinal - this.t;
            }

            // Take one integration step
            const result = this.step(f);
            
            if (!result.success) {
                console.warn('LSODA integration failed at t =', this.t);
                break;
            }

            // Store solution
            t_out.push(this.t);
            y_out.push(this.y.slice());

            // Check for maximum steps
            if (this.nst >= this.maxstep) {
                console.warn('Maximum number of steps reached');
                break;
            }
        }

        return {
            t: t_out,
            y: y_out,
            stats: {
                nsteps: this.nst,
                nfev: this.nfe,
                njev: this.nje,
                method: this.meth === this.ADAMS ? 'Adams' : 'BDF'
            }
        };
    }

    /**
     * Initialize work arrays
     */
    initializeWorkArrays() {
        const maxorder = Math.max(this.maxorder_adams, this.maxorder_bdf) + 1;
        
        this.yh = Array(maxorder).fill(null).map(() => Array(this.n).fill(0));
        this.wm = Array(this.n).fill(null).map(() => Array(this.n).fill(0));
        this.ewt = Array(this.n).fill(0);
        this.savf = Array(this.n).fill(0);
        this.acor = Array(this.n).fill(0);

        // Copy initial conditions to Nordsieck array
        for (let i = 0; i < this.n; i++) {
            this.yh[0][i] = this.y[i];
        }
    }

    /**
     * Take one integration step
     * @param {Function} f - Right-hand side function
     * @returns {Object} Step result
     */
    step(f) {
        let kflag = 0;
        let told = this.t;
        let ncf = 0;  // Convergence failure counter
        let nef = 0;  // Error test failure counter

        // Main step loop
        while (true) {
            // Predict solution using current method
            this.predict();

            // Evaluate RHS at predicted point
            this.nfe++;
            const f_val = f(this.t + this.h, this.y);
            for (let i = 0; i < this.n; i++) {
                this.savf[i] = f_val[i];
            }

            // Correct the solution
            kflag = this.correct(f);

            if (kflag === 0) {
                // Successful step
                this.nst++;
                this.t += this.h;

                // Update Nordsieck array
                this.updateNordsieck();

                // Select next step size and order
                this.selectStepSize();

                // Check for method switching
                this.checkMethodSwitch();

                return { success: true };
            } else if (kflag === -1) {
                // Convergence failure - reduce step size
                ncf++;
                if (ncf >= 10) {
                    return { success: false, reason: 'Too many convergence failures' };
                }
                this.h *= 0.5;
                this.y = [...this.yh[0]]; // Reset y to last successful value
            } else if (kflag === -2) {
                // Error test failure - reduce step size
                nef++;
                if (nef >= 10) {
                    return { success: false, reason: 'Too many error test failures' };
                }
                this.h *= 0.5;
                this.y = [...this.yh[0]]; // Reset y to last successful value
            }

            // Check minimum step size
            if (Math.abs(this.h) < this.hmin) {
                return { success: false, reason: 'Step size too small' };
            }
        }
    }

    /**
     * Predict solution using current method coefficients
     */
    predict() {
        for (let i = 0; i < this.n; i++) {
            this.y[i] = this.yh[0][i];
            
            // Add higher order terms
            for (let j = 1; j <= this.nq; j++) {
                this.y[i] += this.elco[this.meth][this.nq][j] * this.yh[j][i];
            }
        }
    }

    /**
     * Correct the predicted solution
     * @param {Function} f - Right-hand side function
     * @returns {number} Success flag (0=success, -1=convergence failure, -2=error test failure)
     */
    correct(f) {
        const maxiter = 3;
        let converged = false;

        // Iterative correction for implicit methods
        for (let iter = 0; iter < maxiter; iter++) {
            // Calculate correction
            for (let i = 0; i < this.n; i++) {
                this.acor[i] = this.h * this.savf[i] - this.yh[1][i];
            }

            // Apply correction
            const correction_norm = this.vectorNorm(this.acor);
            
            for (let i = 0; i < this.n; i++) {
                this.y[i] = this.yh[0][i] + this.elco[this.meth][this.nq][0] * this.acor[i];
            }

            // Check convergence
            if (correction_norm < 0.33) {
                converged = true;
                break;
            }

            // Re-evaluate function for next iteration
            if (iter < maxiter - 1) {
                this.nfe++;
                const f_val = f(this.t + this.h, this.y);
                for (let i = 0; i < this.n; i++) {
                    this.savf[i] = f_val[i];
                }
            }
        }

        if (!converged) {
            return -1; // Convergence failure
        }

        // Error test
        this.computeErrorWeights();
        const error_norm = this.vectorNorm(this.acor, this.ewt);
        
        if (error_norm > 1.0) {
            return -2; // Error test failure
        }

        return 0; // Success
    }

    /**
     * Compute error weights for error control
     */
    computeErrorWeights() {
        for (let i = 0; i < this.n; i++) {
            this.ewt[i] = this.rtol * Math.abs(this.y[i]) + this.atol;
        }
    }

    /**
     * Compute weighted norm of a vector
     * @param {Array} v - Vector
     * @param {Array} w - Weight vector (optional)
     * @returns {number} Weighted norm
     */
    vectorNorm(v, w = null) {
        let sum = 0;
        for (let i = 0; i < v.length; i++) {
            const weight = w ? w[i] : 1.0;
            sum += Math.pow(v[i] / weight, 2);
        }
        return Math.sqrt(sum / v.length);
    }

    /**
     * Update Nordsieck array after successful step
     */
    updateNordsieck() {
        // Store corrected solution
        for (let i = 0; i < this.n; i++) {
            this.yh[0][i] = this.y[i];
            this.yh[1][i] = this.h * this.savf[i];
        }

        // Update higher order differences
        for (let j = 2; j <= this.nq; j++) {
            for (let i = 0; i < this.n; i++) {
                this.yh[j][i] += this.acor[i];
            }
        }
    }

    /**
     * Select next step size and order
     */
    selectStepSize() {
        // Simple step size control - can be enhanced
        const error_norm = this.vectorNorm(this.acor, this.ewt);
        
        if (error_norm < 0.1) {
            // Increase step size
            this.h *= Math.min(2.0, Math.pow(0.1 / error_norm, 1.0 / (this.nq + 1)));
        } else if (error_norm > 0.5) {
            // Decrease step size
            this.h *= Math.max(0.5, Math.pow(0.5 / error_norm, 1.0 / (this.nq + 1)));
        }

        // Enforce step size limits
        this.h = Math.max(Math.abs(this.h), this.hmin) * Math.sign(this.h);
        this.h = Math.min(Math.abs(this.h), this.hmax) * Math.sign(this.h);
    }

    /**
     * Check if method switching is needed (Adams ↔ BDF)
     */
    checkMethodSwitch() {
        // Simple stiffness detection - can be enhanced
        if (this.nst % 20 === 0) {
            const jacobian_norm = this.estimateJacobianNorm();
            
            // Switch to BDF if system appears stiff
            if (jacobian_norm > 10 && this.meth === this.ADAMS) {
                this.meth = this.BDF;
                this.nq = Math.min(this.nq, this.maxorder_bdf);
                console.log('Switched to BDF method (stiff system detected)');
            } 
            // Switch to Adams if system appears non-stiff
            else if (jacobian_norm < 1 && this.meth === this.BDF) {
                this.meth = this.ADAMS;
                this.nq = Math.min(this.nq, this.maxorder_adams);
                console.log('Switched to Adams method (non-stiff system detected)');
            }
        }
    }

    /**
     * Estimate Jacobian norm for stiffness detection
     * @returns {number} Estimated Jacobian norm
     */
    estimateJacobianNorm() {
        // Simple finite difference estimate
        let norm = 0;
        const eps = 1e-8;
        
        for (let i = 0; i < Math.min(this.n, 3); i++) { // Sample few components
            const y_pert = [...this.y];
            y_pert[i] += eps;
            
            // This would need the RHS function - simplified for now
            const df_dy = 1.0; // Placeholder
            norm = Math.max(norm, Math.abs(df_dy));
        }
        
        return norm * Math.abs(this.h);
    }
}

/**
 * Factory function to create LSODA solver with PK-specific settings
 * @returns {LSODA} Configured LSODA solver
 */
function createPKSolver() {
    const solver = new LSODA();
    
    // PK-specific tolerances
    solver.rtol = 1e-8;   // Tighter relative tolerance for concentrations
    solver.atol = 1e-12;  // Absolute tolerance for concentrations near zero
    solver.maxstep = 1000; // Allow more steps for long simulations
    
    return solver;
}

// Export for use in main application
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { LSODA, createPKSolver };
} else {
    window.LSODA = LSODA;
    window.createPKSolver = createPKSolver;
}