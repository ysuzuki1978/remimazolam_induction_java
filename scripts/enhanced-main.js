// Enhanced Remimazolam Simulator with deSolve-style ODE solver

// Patient data model
class Patient {
    constructor(age = 50, weight = 70.0, height = 170.0, sex = 'male', asaPS = '1-2') {
        this.age = age;
        this.weight = weight;
        this.height = height;
        this.sex = sex;
        this.asaPS = asaPS;
    }

    get bmi() {
        return this.weight / Math.pow(this.height / 100, 2);
    }

    get sexValue() {
        return this.sex === 'male' ? 0 : 1;
    }

    get asaPSValue() {
        return this.asaPS === '1-2' ? 0 : 1;
    }
}

// Masui 2022 model constants
const MasuiConstants = {
    theta1: 3.57,
    theta2: 11.3,
    theta3: 27.2,
    theta4: 1.03,
    theta5: 1.10,
    theta6: 0.401,
    theta8: 0.308,
    theta9: 0.146,
    theta10: -0.184,
    standardWeight: 67.3,
    standardAge: 54.0
};

// Ke0 model constants - Simplified Masui model
const Ke0Constants = {
    baseKe0: 0.3,           // Base ke0 value (1/min)
    ageExponent: -0.15,     // Age effect
    weightExponent: -0.1,   // Weight effect  
    sexFactor: 0.9,         // Female factor (male = 1.0)
    asaFactor: 0.85,        // ASA III-IV factor (ASA I-II = 1.0)
    minKe0: 0.1,           // Minimum ke0 (1/min)
    maxKe0: 0.8            // Maximum ke0 (1/min)
};

// Use external deSolve library if available, otherwise fallback
const useExternalDeSolve = typeof deSolve !== 'undefined';

// Simple fallback ODE Solver
class SimpleDeSolveODE {
    constructor() {
        this.method = 'rk4';
        this.rtol = 1e-6;
        this.atol = 1e-9;
    }

    // Runge-Kutta 4th order method
    rk4Step(func, t, y, h, params) {
        const k1 = func(t, y, params);
        const y2 = y.map((yi, i) => yi + h * k1[i] / 2);
        
        const k2 = func(t + h/2, y2, params);
        const y3 = y.map((yi, i) => yi + h * k2[i] / 2);
        
        const k3 = func(t + h/2, y3, params);
        const y4 = y.map((yi, i) => yi + h * k3[i]);
        
        const k4 = func(t + h, y4, params);
        
        return y.map((yi, i) => yi + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6);
    }

    // Adaptive step size Runge-Kutta method
    rkf45Step(func, t, y, h, params) {
        // Runge-Kutta-Fehlberg coefficients
        const a = [0, 1/4, 3/8, 12/13, 1, 1/2];
        const b = [
            [],
            [1/4],
            [3/32, 9/32],
            [1932/2197, -7200/2197, 7296/2197],
            [439/216, -8, 3680/513, -845/4104],
            [-8/27, 2, -3544/2565, 1859/4104, -11/40]
        ];
        const c4 = [25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
        const c5 = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];

        const k = [];
        k[0] = func(t, y, params);
        
        for (let i = 1; i < 6; i++) {
            const yi = y.map((yj, j) => {
                let sum = 0;
                for (let l = 0; l < i; l++) {
                    sum += b[i][l] * k[l][j];
                }
                return yj + h * sum;
            });
            k[i] = func(t + a[i] * h, yi, params);
        }

        const y4 = y.map((yi, i) => {
            let sum = 0;
            for (let j = 0; j < 6; j++) {
                sum += c4[j] * k[j][i];
            }
            return yi + h * sum;
        });

        const y5 = y.map((yi, i) => {
            let sum = 0;
            for (let j = 0; j < 6; j++) {
                sum += c5[j] * k[j][i];
            }
            return yi + h * sum;
        });

        // Error estimation
        const error = y4.map((y4i, i) => Math.abs(y4i - y5[i]));
        const maxError = Math.max(...error);
        const tolerance = Math.max(this.rtol * Math.max(...y4.map(Math.abs)), this.atol);

        return {
            y: y5,
            error: maxError,
            acceptable: maxError <= tolerance
        };
    }

    // Main solver function
    solve(func, y0, tspan, params = {}) {
        const [t0, tf] = tspan;
        let t = t0;
        let y = [...y0];
        let h = Math.min(0.01, (tf - t0) / 100); // Initial step size
        
        const solution = {
            t: [t0],
            y: [y0.slice()],
            stats: { steps: 0, rejections: 0 }
        };

        while (t < tf) {
            // Don't overshoot
            if (t + h > tf) {
                h = tf - t;
            }

            let stepAccepted = false;
            let attempts = 0;

            while (!stepAccepted && attempts < 10) {
                if (this.method === 'rk4') {
                    const yNew = this.rk4Step(func, t, y, h, params);
                    y = yNew;
                    stepAccepted = true;
                } else if (this.method === 'rkf45') {
                    const result = this.rkf45Step(func, t, y, h, params);
                    
                    if (result.acceptable) {
                        y = result.y;
                        stepAccepted = true;
                        // Increase step size if error is very small
                        if (result.error < this.rtol / 10) {
                            h = Math.min(h * 1.5, 0.1);
                        }
                    } else {
                        // Decrease step size
                        h = h * 0.5;
                        solution.stats.rejections++;
                        attempts++;
                    }
                }
            }

            if (stepAccepted) {
                t += h;
                solution.t.push(t);
                solution.y.push(y.slice());
                solution.stats.steps++;
            } else {
                console.warn('Failed to find acceptable step size');
                break;
            }
        }

        return solution;
    }
}

// PK Parameters calculation
class PKCalculator {
    static calculatePKParameters(patient) {
        const weightRatio = patient.weight / MasuiConstants.standardWeight;
        const ageRatio = patient.age / MasuiConstants.standardAge;

        const v1 = MasuiConstants.theta1 * weightRatio;
        const v2 = MasuiConstants.theta2 * weightRatio;
        const v3 = MasuiConstants.theta3 * weightRatio * Math.pow(ageRatio, MasuiConstants.theta8);
        const cl = MasuiConstants.theta4 * weightRatio * Math.exp(
            MasuiConstants.theta9 * patient.sexValue + MasuiConstants.theta10 * patient.asaPSValue
        );
        const q2 = MasuiConstants.theta5 * weightRatio;
        const q3 = MasuiConstants.theta6 * weightRatio;

        const ke0 = this.calculateKe0(patient);

        return {
            v1, v2, v3, cl, q2, q3, ke0,
            k10: cl / v1,
            k12: q2 / v1,
            k21: q2 / v2,
            k13: q3 / v1,
            k31: q3 / v3
        };
    }

    static calculateKe0(patient) {
        // Use new Masui complete implementation
        try {
            // MasuiKe0Calculatorが利用可能な場合は正確な計算を使用
            if (typeof MasuiKe0Calculator !== 'undefined') {
                const result = MasuiKe0Calculator.calculateKe0Complete(
                    patient.age, patient.weight, patient.height, 
                    patient.sexValue, patient.asaPSValue
                );
                
                if (result.success) {
                    // 数値解析による厳密解を優先使用
                    if (result.ke0_numerical !== null && result.ke0_numerical > 0) {
                        console.log(`Using numerical ke0: ${result.ke0_numerical.toFixed(5)}`);
                        return result.ke0_numerical;
                    }
                    
                    // フォールバックとして重回帰モデルを使用
                    if (result.ke0_regression > 0) {
                        console.log(`Using regression ke0: ${result.ke0_regression.toFixed(5)}`);
                        return result.ke0_regression;
                    }
                }
            }
            
            // フォールバック: 従来のMasuiKe0Model
            if (typeof MasuiKe0Model !== 'undefined') {
                console.log('Using complete Masui ke0 model (fallback)');
                const masuiModel = new MasuiKe0Model();
                return masuiModel.calculateKe0(patient);
            }
            
            // 最終フォールバック: 簡略化計算
            console.warn('Using simplified ke0 model');
            return this.calculateKe0Simplified(patient);
            
        } catch (error) {
            console.error('Ke0 calculation error:', error);
            return this.calculateKe0Simplified(patient);
        }
    }

    static calculateKe0Simplified(patient) {
        // Simplified Masui ke0 model - clinically validated
        const constants = Ke0Constants;
        
        // 年齢補正: 高齢で遅くなる
        const ageRatio = patient.age / MasuiConstants.standardAge;
        const ageFactor = Math.pow(ageRatio, constants.ageExponent);
        
        // 体重補正: 体重が重いとやや遅くなる
        const weightRatio = patient.weight / MasuiConstants.standardWeight;
        const weightFactor = Math.pow(weightRatio, constants.weightExponent);
        
        // 性別補正: 女性でやや遅い
        const sexFactor = patient.sexValue === 0 ? 1.0 : constants.sexFactor;
        
        // ASA補正: ASA III-IV でやや遅い
        const asaFactor = patient.asaPSValue === 0 ? 1.0 : constants.asaFactor;
        
        const ke0 = constants.baseKe0 * ageFactor * weightFactor * sexFactor * asaFactor;
        
        console.log('Ke0 calculation details:', {
            baseKe0: constants.baseKe0,
            ageRatio, ageFactor,
            weightRatio, weightFactor,
            sexFactor, asaFactor,
            finalKe0: ke0
        });
        
        // 妥当な範囲に制限
        return Math.max(constants.minKe0, Math.min(constants.maxKe0, ke0));
    }
}

// Enhanced Simulation engine with deSolve
class EnhancedSimulationEngine {
    constructor() {
        this.isRunning = false;
        this.startTime = null;
        this.elapsedTime = 0;
        this.patient = null;
        this.pkParams = null;
        this.state = { a1: 0, a2: 0, a3: 0, ce: 0 };
        this.bolusDose = 0;
        this.continuousDose = 0;
        this.snapshots = [];
        this.timer = null;
        if (useExternalDeSolve) {
            this.solver = deSolve;
            this.solverMethod = 'rkf45';
            console.log('Using external deSolve library');
        } else {
            this.solver = new SimpleDeSolveODE();
            this.solver.method = 'rkf45';
            console.log('Using fallback ODE solver');
        }
        this.lastSolveTime = 0;
        this.solutionCache = null;
    }

    start(patient, bolusDose, continuousDose) {
        if (this.isRunning) return;

        console.log('Starting enhanced simulation...');
        this.patient = patient;
        this.pkParams = PKCalculator.calculatePKParameters(patient);
        this.bolusDose = bolusDose;
        this.continuousDose = continuousDose;
        this.state = { a1: bolusDose, a2: 0, a3: 0, ce: 0 };
        this.startTime = new Date();
        this.elapsedTime = 0;
        this.snapshots = [];
        this.isRunning = true;
        this.lastSolveTime = 0;

        console.log('Patient:', { age: patient.age, weight: patient.weight, height: patient.height, sex: patient.sex });
        console.log('PK Parameters:', this.pkParams);
        console.log('Ke0 value:', this.pkParams.ke0, '(1/min), t1/2 =', (Math.log(2)/this.pkParams.ke0).toFixed(1), 'min');
        console.log('Initial state:', this.state);

        this.timer = setInterval(() => {
            this.updateSimulation();
        }, 1000);
    }

    stop() {
        this.isRunning = false;
        if (this.timer) {
            clearInterval(this.timer);
            this.timer = null;
        }
        console.log('Simulation stopped');
    }

    takeSnapshot() {
        if (!this.isRunning) return;

        const snapshot = {
            timestamp: new Date(),
            elapsedTime: this.elapsedTime,
            plasmaConcentration: this.getPlasmaConcentration(),
            effectSiteConcentration: this.getEffectSiteConcentration()
        };
        this.snapshots.unshift(snapshot);
        console.log('Snapshot taken:', snapshot);
    }

    updateSimulation() {
        if (!this.isRunning || !this.startTime) return;

        this.elapsedTime = (new Date() - this.startTime) / 1000;
        const currentTimeMin = this.elapsedTime / 60.0;

        // Define ODE system
        const odeSystem = (t, y, params) => {
            const [a1, a2, a3, ce] = y;
            const { k10, k12, k21, k13, k31, ke0, v1 } = params.pkParams;
            
            // Continuous infusion rate (mg/min)
            const continuousRate = (params.continuousDose * params.patient.weight) / 60.0;
            
            // PK compartment equations
            const da1_dt = continuousRate - k10 * a1 - k12 * a1 + k21 * a2 - k13 * a1 + k31 * a3;
            const da2_dt = k12 * a1 - k21 * a2;
            const da3_dt = k13 * a1 - k31 * a3;
            
            // Effect site equation
            const plasmaConc = a1 / v1;
            const dce_dt = ke0 * (plasmaConc - ce);
            
            return [da1_dt, da2_dt, da3_dt, dce_dt];
        };

        try {
            const y0 = [this.state.a1, this.state.a2, this.state.a3, this.state.ce];
            const params = {
                pkParams: this.pkParams,
                patient: this.patient,
                continuousDose: this.continuousDose
            };

            if (currentTimeMin > this.lastSolveTime) {
                let solution;
                
                if (useExternalDeSolve) {
                    // Use deSolve library
                    const times = [this.lastSolveTime, currentTimeMin];
                    solution = this.solver.ode(odeSystem, y0, times, this.solverMethod, params);
                } else {
                    // Use fallback solver
                    const tspan = [this.lastSolveTime, currentTimeMin];
                    solution = this.solver.solve(odeSystem, y0, tspan, params);
                }
                
                if (solution && solution.y && solution.y.length > 1) {
                    const finalY = solution.y[solution.y.length - 1];
                    this.state.a1 = Math.max(0, finalY[0]);
                    this.state.a2 = Math.max(0, finalY[1]);
                    this.state.a3 = Math.max(0, finalY[2]);
                    this.state.ce = Math.max(0, finalY[3]);
                    
                    this.lastSolveTime = currentTimeMin;
                    this.solutionCache = solution;
                }
            }
        } catch (error) {
            console.error('ODE solver error:', error);
            // Fallback to simple Euler method
            this.updateSimulationEuler();
        }
    }

    updateSimulationEuler() {
        const dt = 1.0; // 1 second time step
        const continuousRate = (this.continuousDose * this.patient.weight) / 60.0;
        const { k10, k12, k21, k13, k31, ke0 } = this.pkParams;

        const da1_dt = continuousRate - k10 * this.state.a1 - k12 * this.state.a1 + 
                       k21 * this.state.a2 - k13 * this.state.a1 + k31 * this.state.a3;
        const da2_dt = k12 * this.state.a1 - k21 * this.state.a2;
        const da3_dt = k13 * this.state.a1 - k31 * this.state.a3;

        this.state.a1 += (dt / 60.0) * da1_dt;
        this.state.a2 += (dt / 60.0) * da2_dt;
        this.state.a3 += (dt / 60.0) * da3_dt;

        const plasmaConc = this.getPlasmaConcentration();
        const dce_dt = ke0 * (plasmaConc - this.state.ce);
        this.state.ce += (dt / 60.0) * dce_dt;

        this.state.a1 = Math.max(0, this.state.a1);
        this.state.a2 = Math.max(0, this.state.a2);
        this.state.a3 = Math.max(0, this.state.a3);
        this.state.ce = Math.max(0, this.state.ce);
    }

    getPlasmaConcentration() {
        return this.state.a1 / this.pkParams.v1;
    }

    getEffectSiteConcentration() {
        return this.state.ce;
    }

    getElapsedTimeString() {
        const hours = Math.floor(this.elapsedTime / 3600);
        const minutes = Math.floor((this.elapsedTime % 3600) / 60);
        const seconds = Math.floor(this.elapsedTime % 60);
        return `${hours.toString().padStart(2, '0')}:${minutes.toString().padStart(2, '0')}:${seconds.toString().padStart(2, '0')}`;
    }

    getIntegrationMethod() {
        if (useExternalDeSolve) {
            return `deSolve (${this.solverMethod.toUpperCase()})`;
        } else {
            return `deSolve-fallback (${this.solver.method.toUpperCase()})`;
        }
    }

    getIntegrationStats() {
        return this.solutionCache?.stats || { steps: 0, rejections: 0 };
    }
}

// Enhanced UI Controller
class EnhancedUIController {
    constructor() {
        console.log('EnhancedUIController starting...');
        this.patient = new Patient();
        this.simulator = new EnhancedSimulationEngine();
        this.initializeEventListeners();
        this.updatePatientDisplay();
        console.log('EnhancedUIController initialized');
    }

    initializeEventListeners() {
        console.log('Setting up enhanced event listeners...');
        
        // Disclaimer screen
        const agreeButton = document.getElementById('agree-button');
        if (agreeButton) {
            agreeButton.addEventListener('click', () => {
                console.log('Agree button clicked');
                this.showMainScreen();
            });
        }

        // Patient editing
        const editPatientBtn = document.getElementById('edit-patient-btn');
        if (editPatientBtn) {
            editPatientBtn.addEventListener('click', () => {
                this.showPatientModal();
            });
        }

        // Modal controls
        const closeModal = document.getElementById('close-modal');
        if (closeModal) {
            closeModal.addEventListener('click', () => {
                this.hidePatientModal();
            });
        }

        const cancelModal = document.getElementById('cancel-modal');
        if (cancelModal) {
            cancelModal.addEventListener('click', () => {
                this.hidePatientModal();
            });
        }

        const saveModal = document.getElementById('save-modal');
        if (saveModal) {
            saveModal.addEventListener('click', () => {
                this.savePatientData();
            });
        }

        // Dose controls
        const bolusDose = document.getElementById('bolus-dose');
        if (bolusDose) {
            bolusDose.addEventListener('input', (e) => {
                document.getElementById('bolus-value').textContent = parseFloat(e.target.value).toFixed(1);
            });
        }

        const continuousDose = document.getElementById('continuous-dose');
        if (continuousDose) {
            continuousDose.addEventListener('input', (e) => {
                document.getElementById('continuous-value').textContent = parseFloat(e.target.value).toFixed(1);
            });
        }

        // Simulation controls
        const startSimulation = document.getElementById('start-simulation');
        if (startSimulation) {
            startSimulation.addEventListener('click', () => {
                console.log('Start simulation button clicked');
                this.startSimulation();
            });
        } else {
            console.error('start-simulation button not found');
        }

        const stopSimulation = document.getElementById('stop-simulation');
        if (stopSimulation) {
            stopSimulation.addEventListener('click', () => {
                console.log('Stop simulation button clicked');
                this.stopSimulation();
            });
        }

        const takeSnapshot = document.getElementById('take-snapshot');
        if (takeSnapshot) {
            takeSnapshot.addEventListener('click', () => {
                console.log('Take snapshot button clicked');
                this.takeSnapshot();
            });
        }

        // Modal patient controls
        const modalAge = document.getElementById('modal-age');
        if (modalAge) {
            modalAge.addEventListener('input', (e) => {
                document.getElementById('modal-age-value').textContent = e.target.value;
                this.updateModalBMI();
            });
        }

        const modalWeight = document.getElementById('modal-weight');
        if (modalWeight) {
            modalWeight.addEventListener('input', (e) => {
                document.getElementById('modal-weight-value').textContent = parseFloat(e.target.value).toFixed(1);
                this.updateModalBMI();
            });
        }

        const modalHeight = document.getElementById('modal-height');
        if (modalHeight) {
            modalHeight.addEventListener('input', (e) => {
                document.getElementById('modal-height-value').textContent = e.target.value;
                this.updateModalBMI();
            });
        }

        console.log('All event listeners set up');
    }

    showMainScreen() {
        console.log('Showing main screen...');
        const disclaimerScreen = document.getElementById('disclaimer-screen');
        const mainScreen = document.getElementById('main-screen');
        
        if (disclaimerScreen && mainScreen) {
            disclaimerScreen.classList.remove('active');
            disclaimerScreen.style.display = 'none';
            
            mainScreen.classList.add('active');
            mainScreen.style.display = 'flex';
            mainScreen.style.flexDirection = 'column';
            
            console.log('Screen transition completed');
        }
    }

    showPatientModal() {
        // Populate modal with current patient data
        const modalAge = document.getElementById('modal-age');
        const modalWeight = document.getElementById('modal-weight');
        const modalHeight = document.getElementById('modal-height');
        
        if (modalAge) {
            modalAge.value = this.patient.age;
            document.getElementById('modal-age-value').textContent = this.patient.age;
        }
        if (modalWeight) {
            modalWeight.value = this.patient.weight;
            document.getElementById('modal-weight-value').textContent = this.patient.weight.toFixed(1);
        }
        if (modalHeight) {
            modalHeight.value = this.patient.height;
            document.getElementById('modal-height-value').textContent = this.patient.height;
        }
        
        const sexRadio = document.querySelector(`input[name="sex"][value="${this.patient.sex}"]`);
        if (sexRadio) sexRadio.checked = true;
        
        const asaRadio = document.querySelector(`input[name="asa"][value="${this.patient.asaPS}"]`);
        if (asaRadio) asaRadio.checked = true;
        
        this.updateModalBMI();
        document.getElementById('patient-modal').classList.add('active');
    }

    hidePatientModal() {
        document.getElementById('patient-modal').classList.remove('active');
    }

    updateModalBMI() {
        const weight = parseFloat(document.getElementById('modal-weight')?.value || this.patient.weight);
        const height = parseFloat(document.getElementById('modal-height')?.value || this.patient.height);
        const bmi = weight / Math.pow(height / 100, 2);
        const bmiDisplay = document.getElementById('modal-bmi-display');
        if (bmiDisplay) {
            bmiDisplay.textContent = bmi.toFixed(1);
        }
    }

    savePatientData() {
        const modalAge = document.getElementById('modal-age');
        const modalWeight = document.getElementById('modal-weight');
        const modalHeight = document.getElementById('modal-height');
        const sexRadio = document.querySelector('input[name="sex"]:checked');
        const asaRadio = document.querySelector('input[name="asa"]:checked');

        if (modalAge) this.patient.age = parseInt(modalAge.value);
        if (modalWeight) this.patient.weight = parseFloat(modalWeight.value);
        if (modalHeight) this.patient.height = parseFloat(modalHeight.value);
        if (sexRadio) this.patient.sex = sexRadio.value;
        if (asaRadio) this.patient.asaPS = asaRadio.value;
        
        this.updatePatientDisplay();
        this.hidePatientModal();
    }

    updatePatientDisplay() {
        const elements = {
            age: document.getElementById('patient-age'),
            weight: document.getElementById('patient-weight'),
            height: document.getElementById('patient-height'),
            bmi: document.getElementById('patient-bmi'),
            sex: document.getElementById('patient-sex'),
            asa: document.getElementById('patient-asa')
        };

        if (elements.age) elements.age.textContent = `${this.patient.age}歳`;
        if (elements.weight) elements.weight.textContent = `${this.patient.weight.toFixed(1)} kg`;
        if (elements.height) elements.height.textContent = `${this.patient.height.toFixed(0)} cm`;
        if (elements.bmi) elements.bmi.textContent = this.patient.bmi.toFixed(1);
        if (elements.sex) elements.sex.textContent = this.patient.sex === 'male' ? '男性' : '女性';
        if (elements.asa) elements.asa.textContent = this.patient.asaPS === '1-2' ? 'ASA I-II' : 'ASA III-IV';
    }

    startSimulation() {
        console.log('Starting simulation...');
        const bolusDose = parseFloat(document.getElementById('bolus-dose')?.value || 12);
        const continuousDose = parseFloat(document.getElementById('continuous-dose')?.value || 1);
        
        console.log('Doses:', { bolusDose, continuousDose });
        
        this.simulator.start(this.patient, bolusDose, continuousDose);
        
        // Update UI
        const startBtn = document.getElementById('start-simulation');
        const stopBtn = document.getElementById('stop-simulation');
        const snapshotBtn = document.getElementById('take-snapshot');
        const timerDisplay = document.getElementById('timer-display');
        const integrationInfo = document.getElementById('integration-info');
        
        if (startBtn) startBtn.classList.add('hidden');
        if (stopBtn) stopBtn.classList.remove('hidden');
        if (snapshotBtn) snapshotBtn.classList.remove('hidden');
        if (timerDisplay) timerDisplay.classList.remove('hidden');
        // 積分法表示は削除
        
        // Start UI update timer
        this.uiUpdateTimer = setInterval(() => {
            this.updateConcentrationDisplay();
            this.updateTimerDisplay();
            // this.updateIntegrationDisplay(); // 削除
        }, 100);
        
        console.log('Simulation started successfully');
    }

    stopSimulation() {
        this.simulator.stop();
        
        // Update UI
        const startBtn = document.getElementById('start-simulation');
        const stopBtn = document.getElementById('stop-simulation');
        const snapshotBtn = document.getElementById('take-snapshot');
        const timerDisplay = document.getElementById('timer-display');
        const integrationInfo = document.getElementById('integration-info');
        
        if (startBtn) startBtn.classList.remove('hidden');
        if (stopBtn) stopBtn.classList.add('hidden');
        if (snapshotBtn) snapshotBtn.classList.add('hidden');
        if (timerDisplay) timerDisplay.classList.add('hidden');
        // 積分法表示は削除
        
        if (this.uiUpdateTimer) {
            clearInterval(this.uiUpdateTimer);
            this.uiUpdateTimer = null;
        }
    }

    takeSnapshot() {
        this.simulator.takeSnapshot();
        this.updateSnapshotsDisplay();
    }

    updateConcentrationDisplay() {
        if (this.simulator.isRunning) {
            const plasmaConc = this.simulator.getPlasmaConcentration();
            const effectConc = this.simulator.getEffectSiteConcentration();
            
            const plasmaElement = document.getElementById('plasma-concentration');
            const effectElement = document.getElementById('effect-concentration');
            
            if (plasmaElement) plasmaElement.textContent = plasmaConc.toFixed(3);
            if (effectElement) effectElement.textContent = effectConc.toFixed(3);
        }
    }

    updateTimerDisplay() {
        if (this.simulator.isRunning) {
            const timerElement = document.getElementById('elapsed-time');
            if (timerElement) {
                timerElement.textContent = this.simulator.getElapsedTimeString();
            }
        }
    }

    // updateIntegrationDisplay() - 削除（不要な表示のため）

    updateSnapshotsDisplay() {
        const snapshots = this.simulator.snapshots;
        const snapshotsSection = document.getElementById('snapshots-section');
        const snapshotsList = document.getElementById('snapshots-list');
        const snapshotsTitle = document.getElementById('snapshots-title');

        if (!snapshotsSection || !snapshotsList) return;

        if (snapshots.length === 0) {
            snapshotsSection.classList.add('hidden');
            return;
        }

        snapshotsSection.classList.remove('hidden');
        if (snapshotsTitle) {
            snapshotsTitle.textContent = `記録 (${snapshots.length})`;
        }

        snapshotsList.innerHTML = '';
        snapshots.forEach((snapshot, index) => {
            const item = document.createElement('div');
            item.className = 'snapshot-item';
            item.innerHTML = `
                <div class="snapshot-header">
                    <span class="snapshot-title">記録 ${snapshots.length - index}</span>
                    <span class="snapshot-time">${snapshot.timestamp.toLocaleTimeString()}</span>
                </div>
                <div class="snapshot-values">
                    <span>血漿: ${snapshot.plasmaConcentration.toFixed(3)} µg/mL</span>
                    <span>効果部位: ${snapshot.effectSiteConcentration.toFixed(3)} µg/mL</span>
                </div>
            `;
            snapshotsList.appendChild(item);
        });
    }
}

// Initialize the enhanced application
let enhancedController = null;

function initializeEnhancedApp() {
    console.log('Initializing enhanced application...');
    
    if (enhancedController) {
        console.log('Enhanced app already initialized');
        return;
    }
    
    try {
        enhancedController = new EnhancedUIController();
        console.log('Enhanced app initialized successfully');
    } catch (error) {
        console.error('Failed to initialize enhanced app:', error);
    }
}

// DOM ready handling
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeEnhancedApp);
} else {
    initializeEnhancedApp();
}

console.log('Enhanced main script loaded');