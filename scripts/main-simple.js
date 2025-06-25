// 簡略版メインアプリケーション - デバッグ用

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
}

// Simple UI Controller
class SimpleUIController {
    constructor() {
        console.log('SimpleUIController starting...');
        this.patient = new Patient();
        this.initializeEventListeners();
        this.updatePatientDisplay();
        console.log('SimpleUIController initialized');
    }

    initializeEventListeners() {
        console.log('Setting up event listeners...');
        
        // Disclaimer screen
        const agreeButton = document.getElementById('agree-button');
        if (agreeButton) {
            agreeButton.addEventListener('click', () => {
                console.log('Agree button clicked');
                this.showMainScreen();
            });
            console.log('Agree button listener added');
        } else {
            console.error('agree-button not found');
        }
    }

    showMainScreen() {
        console.log('Attempting to show main screen...');
        const disclaimerScreen = document.getElementById('disclaimer-screen');
        const mainScreen = document.getElementById('main-screen');
        
        if (disclaimerScreen && mainScreen) {
            disclaimerScreen.classList.remove('active');
            disclaimerScreen.style.display = 'none';
            
            mainScreen.classList.add('active');
            mainScreen.style.display = 'flex';
            mainScreen.style.flexDirection = 'column';
            
            console.log('Screen transition completed');
        } else {
            console.error('Screen elements not found:', {
                disclaimerScreen: !!disclaimerScreen,
                mainScreen: !!mainScreen
            });
        }
    }

    updatePatientDisplay() {
        console.log('Updating patient display...');
        const elements = {
            age: document.getElementById('patient-age'),
            weight: document.getElementById('patient-weight'),
            height: document.getElementById('patient-height'),
            bmi: document.getElementById('patient-bmi'),
            sex: document.getElementById('patient-sex'),
            asa: document.getElementById('patient-asa')
        };

        // 安全に更新
        if (elements.age) elements.age.textContent = `${this.patient.age}歳`;
        if (elements.weight) elements.weight.textContent = `${this.patient.weight.toFixed(1)} kg`;
        if (elements.height) elements.height.textContent = `${this.patient.height.toFixed(0)} cm`;
        if (elements.bmi) elements.bmi.textContent = this.patient.bmi.toFixed(1);
        if (elements.sex) elements.sex.textContent = this.patient.sex === 'male' ? '男性' : '女性';
        if (elements.asa) elements.asa.textContent = this.patient.asaPS === '1-2' ? 'ASA I-II' : 'ASA III-IV';
        
        console.log('Patient display updated');
    }
}

// 単一の初期化ポイント
let appController = null;

function initializeApp() {
    console.log('Initializing application...');
    
    if (appController) {
        console.log('App already initialized');
        return;
    }
    
    try {
        appController = new SimpleUIController();
        console.log('App initialized successfully');
    } catch (error) {
        console.error('Failed to initialize app:', error);
        
        // フォールバック：最小限の機能
        const agreeButton = document.getElementById('agree-button');
        if (agreeButton) {
            agreeButton.addEventListener('click', function() {
                console.log('Fallback: Agree button clicked');
                document.getElementById('disclaimer-screen').style.display = 'none';
                document.getElementById('main-screen').style.display = 'flex';
            });
        }
    }
}

// DOM ready handling
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initializeApp);
} else {
    initializeApp();
}

console.log('Simple main script loaded');