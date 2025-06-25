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

// Ke0 model constants
const Ke0Constants = {
    baseLogKe0: -9.06,
    ageAmplitude: 5.44,
    ageMean: 30.0,
    ageStandardDeviation: 15.0,
    weightAmplitude: 2.1,
    weightMean: 70.0,
    weightStandardDeviation: 20.0,
    heightAmplitude: 1.8,
    heightMean: 165.0,
    heightStandardDeviation: 15.0,
    bmiAmplitude: 1.2,
    bmiMean: 22.0,
    bmiStandardDeviation: 5.0,
    genderEffect: 0.15,
    asaEffect: 0.08,
    timeConversionFactor: 60.0
};

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
        const age = patient.age;
        const weight = patient.weight;
        const height = patient.height;
        const bmi = patient.bmi;
        const sex = patient.sexValue;
        const asaPS = patient.asaPSValue;

        const ageEffect = Ke0Constants.ageAmplitude * Math.exp(
            -Math.pow((age - Ke0Constants.ageMean) / Ke0Constants.ageStandardDeviation, 2) / 2
        );
        const weightEffect = Ke0Constants.weightAmplitude * Math.exp(
            -Math.pow((weight - Ke0Constants.weightMean) / Ke0Constants.weightStandardDeviation, 2) / 2
        );
        const heightEffect = Ke0Constants.heightAmplitude * Math.exp(
            -Math.pow((height - Ke0Constants.heightMean) / Ke0Constants.heightStandardDeviation, 2) / 2
        );
        const bmiEffect = Ke0Constants.bmiAmplitude * Math.exp(
            -Math.pow((bmi - Ke0Constants.bmiMean) / Ke0Constants.bmiStandardDeviation, 2) / 2
        );

        const logKe0 = Ke0Constants.baseLogKe0 + ageEffect + weightEffect + heightEffect + bmiEffect +
                      Ke0Constants.genderEffect * sex + Ke0Constants.asaEffect * asaPS;

        return Math.exp(logKe0) / Ke0Constants.timeConversionFactor;
    }
}

// Enhanced Simulation engine with LSODA integration
class SimulationEngine {
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
        this.useLSODA = false; // Disable LSODA for debugging
        this.lsodaSolver = null;
        this.integrationStats = null;
    }

    start(patient, bolusDose, continuousDose) {
        if (this.isRunning) return;

        this.patient = patient;
        this.pkParams = PKCalculator.calculatePKParameters(patient);
        this.bolusDose = bolusDose;
        this.continuousDose = continuousDose;
        this.state = { a1: bolusDose, a2: 0, a3: 0, ce: 0 };
        this.startTime = new Date();
        this.elapsedTime = 0;
        this.snapshots = [];
        this.isRunning = true;

        // Initialize LSODA solver if enabled
        if (this.useLSODA && typeof LSODA !== 'undefined') {
            this.lsodaSolver = createPKSolver();
            console.log('Using LSODA integration method');
        } else {
            console.log('Using Euler integration method');
        }

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
    }

    updateSimulation() {
        if (!this.isRunning || !this.startTime) return;

        this.elapsedTime = (new Date() - this.startTime) / 1000;
        
        if (this.useLSODA && this.lsodaSolver) {
            this.updateSimulationLSODA();
        } else {
            this.updateSimulationEuler();
        }
    }

    updateSimulationEuler() {
        const dt = 1.0; // 1 second time step

        // Continuous infusion rate (mg/min)
        const continuousRate = (this.continuousDose * this.patient.weight) / 60.0;

        // Rate constants
        const { k10, k12, k21, k13, k31, ke0 } = this.pkParams;

        // Differential equations
        const da1_dt = continuousRate - k10 * this.state.a1 - k12 * this.state.a1 + 
                       k21 * this.state.a2 - k13 * this.state.a1 + k31 * this.state.a3;
        const da2_dt = k12 * this.state.a1 - k21 * this.state.a2;
        const da3_dt = k13 * this.state.a1 - k31 * this.state.a3;

        // Update compartments (Euler integration)
        this.state.a1 += (dt / 60.0) * da1_dt; // Convert dt to minutes
        this.state.a2 += (dt / 60.0) * da2_dt;
        this.state.a3 += (dt / 60.0) * da3_dt;

        // Effect site concentration
        const plasmaConc = this.getPlasmaConcentration();
        const dce_dt = ke0 * (plasmaConc - this.state.ce);
        this.state.ce += (dt / 60.0) * dce_dt;

        // Ensure non-negative values
        this.state.a1 = Math.max(0, this.state.a1);
        this.state.a2 = Math.max(0, this.state.a2);
        this.state.a3 = Math.max(0, this.state.a3);
        this.state.ce = Math.max(0, this.state.ce);
    }

    updateSimulationLSODA() {
        const currentTimeMin = this.elapsedTime / 60.0;
        const nextTimeMin = currentTimeMin + 1.0/60.0; // 1 second step in minutes
        
        // Define the ODE system for LSODA
        const odeSystem = (t, y) => {
            const [a1, a2, a3, ce] = y;
            const { k10, k12, k21, k13, k31, ke0, v1 } = this.pkParams;
            
            // Continuous infusion rate (mg/min)
            const continuousRate = (this.continuousDose * this.patient.weight) / 60.0;
            
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
            // Solve ODE from current time to next time step
            const y0 = [this.state.a1, this.state.a2, this.state.a3, this.state.ce];
            const result = this.lsodaSolver.solve(odeSystem, y0, currentTimeMin, nextTimeMin, {
                rtol: 1e-8,
                atol: 1e-12,
                maxstep: 100
            });
            
            if (result.t.length > 1) {
                // Update state with the latest solution
                const finalY = result.y[result.y.length - 1];
                this.state.a1 = Math.max(0, finalY[0]);
                this.state.a2 = Math.max(0, finalY[1]);
                this.state.a3 = Math.max(0, finalY[2]);
                this.state.ce = Math.max(0, finalY[3]);
                
                // Store integration statistics
                this.integrationStats = result.stats;
            }
        } catch (error) {
            console.warn('LSODA integration failed, falling back to Euler method:', error);
            this.useLSODA = false;
            this.updateSimulationEuler();
        }
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
        if (this.useLSODA && this.integrationStats) {
            return `LSODA (${this.integrationStats.method})`;
        }
        return this.useLSODA ? 'LSODA (initializing)' : 'Euler';
    }

    getIntegrationStats() {
        return this.integrationStats || {
            nsteps: 0,
            nfev: 0,
            njev: 0,
            method: 'Euler'
        };
    }
}

// UI Controller
class UIController {
    constructor() {
        this.patient = new Patient();
        this.simulator = new SimulationEngine();
        this.initializeEventListeners();
        this.updatePatientDisplay();
    }

    initializeEventListeners() {
        console.log('Initializing event listeners...');
        
        // Disclaimer screen
        const agreeButton = document.getElementById('agree-button');
        if (agreeButton) {
            agreeButton.addEventListener('click', () => {
                console.log('Agree button clicked');
                this.showMainScreen();
            });
        } else {
            console.error('agree-button not found');
        }

        // Patient editing
        document.getElementById('edit-patient-btn').addEventListener('click', () => {
            this.showPatientModal();
        });

        document.getElementById('close-modal').addEventListener('click', () => {
            this.hidePatientModal();
        });

        document.getElementById('cancel-modal').addEventListener('click', () => {
            this.hidePatientModal();
        });

        document.getElementById('save-modal').addEventListener('click', () => {
            this.savePatientData();
        });

        // Dose controls
        document.getElementById('bolus-dose').addEventListener('input', (e) => {
            document.getElementById('bolus-value').textContent = parseFloat(e.target.value).toFixed(1);
        });

        document.getElementById('continuous-dose').addEventListener('input', (e) => {
            document.getElementById('continuous-value').textContent = parseFloat(e.target.value).toFixed(1);
        });

        // Modal patient controls
        document.getElementById('modal-age').addEventListener('input', (e) => {
            document.getElementById('modal-age-value').textContent = e.target.value;
            this.updateModalBMI();
        });

        document.getElementById('modal-weight').addEventListener('input', (e) => {
            document.getElementById('modal-weight-value').textContent = parseFloat(e.target.value).toFixed(1);
            this.updateModalBMI();
        });

        document.getElementById('modal-height').addEventListener('input', (e) => {
            document.getElementById('modal-height-value').textContent = e.target.value;
            this.updateModalBMI();
        });

        // Simulation controls
        document.getElementById('start-simulation').addEventListener('click', () => {
            this.startSimulation();
        });

        document.getElementById('stop-simulation').addEventListener('click', () => {
            this.stopSimulation();
        });

        document.getElementById('take-snapshot').addEventListener('click', () => {
            this.takeSnapshot();
        });

        // Modal click outside to close
        document.getElementById('patient-modal').addEventListener('click', (e) => {
            if (e.target.id === 'patient-modal') {
                this.hidePatientModal();
            }
        });
    }

    showMainScreen() {
        console.log('Showing main screen...');
        const disclaimerScreen = document.getElementById('disclaimer-screen');
        const mainScreen = document.getElementById('main-screen');
        
        if (disclaimerScreen && mainScreen) {
            disclaimerScreen.classList.remove('active');
            mainScreen.classList.add('active');
            console.log('Screen transition completed');
        } else {
            console.error('Screen elements not found:', {
                disclaimerScreen: !!disclaimerScreen,
                mainScreen: !!mainScreen
            });
        }
    }

    showPatientModal() {
        // Populate modal with current patient data
        document.getElementById('modal-age').value = this.patient.age;
        document.getElementById('modal-age-value').textContent = this.patient.age;
        document.getElementById('modal-weight').value = this.patient.weight;
        document.getElementById('modal-weight-value').textContent = this.patient.weight.toFixed(1);
        document.getElementById('modal-height').value = this.patient.height;
        document.getElementById('modal-height-value').textContent = this.patient.height;
        
        document.querySelector(`input[name="sex"][value="${this.patient.sex}"]`).checked = true;
        document.querySelector(`input[name="asa"][value="${this.patient.asaPS}"]`).checked = true;
        
        this.updateModalBMI();
        document.getElementById('patient-modal').classList.add('active');
    }

    hidePatientModal() {
        document.getElementById('patient-modal').classList.remove('active');
    }

    updateModalBMI() {
        const weight = parseFloat(document.getElementById('modal-weight').value);
        const height = parseFloat(document.getElementById('modal-height').value);
        const bmi = weight / Math.pow(height / 100, 2);
        document.getElementById('modal-bmi-display').textContent = bmi.toFixed(1);
    }

    savePatientData() {
        this.patient.age = parseInt(document.getElementById('modal-age').value);
        this.patient.weight = parseFloat(document.getElementById('modal-weight').value);
        this.patient.height = parseFloat(document.getElementById('modal-height').value);
        this.patient.sex = document.querySelector('input[name="sex"]:checked').value;
        this.patient.asaPS = document.querySelector('input[name="asa"]:checked').value;
        
        this.updatePatientDisplay();
        this.hidePatientModal();
    }

    updatePatientDisplay() {
        document.getElementById('patient-age').textContent = `${this.patient.age}歳`;
        document.getElementById('patient-weight').textContent = `${this.patient.weight.toFixed(1)} kg`;
        document.getElementById('patient-height').textContent = `${this.patient.height.toFixed(0)} cm`;
        document.getElementById('patient-bmi').textContent = this.patient.bmi.toFixed(1);
        document.getElementById('patient-sex').textContent = this.patient.sex === 'male' ? '男性' : '女性';
        document.getElementById('patient-asa').textContent = this.patient.asaPS === '1-2' ? 'ASA I-II' : 'ASA III-IV';
    }

    startSimulation() {
        const bolusDose = parseFloat(document.getElementById('bolus-dose').value);
        const continuousDose = parseFloat(document.getElementById('continuous-dose').value);
        
        this.simulator.start(this.patient, bolusDose, continuousDose);
        
        // Update UI
        document.getElementById('start-simulation').classList.add('hidden');
        document.getElementById('stop-simulation').classList.remove('hidden');
        document.getElementById('take-snapshot').classList.remove('hidden');
        document.getElementById('timer-display').classList.remove('hidden');
        document.getElementById('integration-info').classList.remove('hidden');
        
        // Start UI update timer
        this.uiUpdateTimer = setInterval(() => {
            this.updateConcentrationDisplay();
            this.updateTimerDisplay();
            this.updateIntegrationDisplay();
        }, 100);
    }

    stopSimulation() {
        this.simulator.stop();
        
        // Update UI
        document.getElementById('start-simulation').classList.remove('hidden');
        document.getElementById('stop-simulation').classList.add('hidden');
        document.getElementById('take-snapshot').classList.add('hidden');
        document.getElementById('timer-display').classList.add('hidden');
        document.getElementById('integration-info').classList.add('hidden');
        
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
            
            document.getElementById('plasma-concentration').textContent = plasmaConc.toFixed(3);
            document.getElementById('effect-concentration').textContent = effectConc.toFixed(3);
        }
    }

    updateTimerDisplay() {
        if (this.simulator.isRunning) {
            document.getElementById('elapsed-time').textContent = this.simulator.getElapsedTimeString();
        }
    }

    updateIntegrationDisplay() {
        if (this.simulator.isRunning) {
            document.getElementById('integration-method').textContent = this.simulator.getIntegrationMethod();
            const stats = this.simulator.getIntegrationStats();
            document.getElementById('integration-steps').textContent = stats.nsteps || 0;
            document.getElementById('integration-fev').textContent = stats.nfev || 0;
        }
    }

    updateSnapshotsDisplay() {
        const snapshots = this.simulator.snapshots;
        const snapshotsSection = document.getElementById('snapshots-section');
        const snapshotsList = document.getElementById('snapshots-list');
        const snapshotsTitle = document.getElementById('snapshots-title');

        if (snapshots.length === 0) {
            snapshotsSection.classList.add('hidden');
            return;
        }

        snapshotsSection.classList.remove('hidden');
        snapshotsTitle.textContent = `記録 (${snapshots.length})`;

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

// Prevent multiple initializations
let uiController = null;

// Initialize the application
document.addEventListener('DOMContentLoaded', () => {
    console.log('DOM loaded, initializing application...');
    if (!uiController) {
        try {
            uiController = new UIController();
            console.log('UIController initialized successfully');
        } catch (error) {
            console.error('Failed to initialize UIController:', error);
        }
    }
});

// Fallback initialization if DOMContentLoaded already fired
if (document.readyState === 'loading') {
    // DOM is still loading, wait for DOMContentLoaded
} else {
    // DOM is already loaded
    console.log('DOM already loaded, initializing application...');
    if (!uiController) {
        try {
            uiController = new UIController();
            console.log('UIController initialized successfully');
        } catch (error) {
            console.error('Failed to initialize UIController:', error);
        }
    }
}