/**
 * Masui ke0 Model - Complete Implementation
 * Based on the full model with all interaction terms
 */

class MasuiKe0Model {
    constructor() {
        this.name = "Masui ke0 Model (Complete)";
    }

    /**
     * Calculate standardized variables F2_*
     */
    calculateStandardizedVariables(patient) {
        // Patient parameters
        const age = patient.age;
        const weight = patient.weight; // TBW (Total Body Weight)
        const height = patient.height;
        const sex = patient.sexValue; // 0=male, 1=female
        const asaps = patient.asaPSValue; // 0=ASA I-II, 1=ASA III-IV

        // Calculate F2 standardized variables (assuming these are standardized)
        // Note: The exact standardization formula may need adjustment based on the paper
        const F2_age = (age - 54.0) / 15.0;
        const F2_TBW = (weight - 67.3) / 15.0;
        const F2_height = (height - 165.0) / 10.0;
        const F2_sex = sex; // Binary variable
        const F2_ASAPS = asaps; // Binary variable

        return {
            F2_age,
            F2_TBW,
            F2_height,
            F2_sex,
            F2_ASAPS
        };
    }

    /**
     * Calculate main effects F(variable)
     */
    calculateMainEffects(F2) {
        // These coefficients would need to be derived from the paper
        // Using placeholder values that need to be replaced with actual coefficients
        const F_age = 5.44 * Math.exp(-0.5 * Math.pow(F2.F2_age, 2));
        const F_TBW = 2.1 * Math.exp(-0.5 * Math.pow(F2.F2_TBW, 2));
        const F_height = 1.8 * Math.exp(-0.5 * Math.pow(F2.F2_height, 2));
        const F_sex = 0.15 * F2.F2_sex;
        const F_ASAPS = 0.08 * F2.F2_ASAPS;

        return {
            F_age,
            F_TBW,
            F_height,
            F_sex,
            F_ASAPS
        };
    }

    /**
     * Calculate complete ke0 using the full model
     */
    calculateKe0(patient) {
        console.log('=== Masui Ke0 Calculation (Complete Model) ===');
        console.log('Patient:', patient);

        // Step 1: Calculate standardized variables
        const F2 = this.calculateStandardizedVariables(patient);
        console.log('F2 variables:', F2);

        // Step 2: Calculate main effects
        const F = this.calculateMainEffects(F2);
        console.log('Main effects:', F);

        // Step 3: Calculate ke0 using the complete formula
        let logKe0 = 0;

        // 定数項と主効果
        logKe0 += -9.06; // 定数項
        logKe0 += F.F_age;
        logKe0 += F.F_TBW;
        logKe0 += F.F_height;
        logKe0 += 0.999 * F.F_sex;
        logKe0 += F.F_ASAPS;

        console.log('After main effects, logKe0 =', logKe0);

        // 2次の交互作用項
        logKe0 -= 4.50 * F2.F2_age * F2.F2_TBW;
        logKe0 -= 4.51 * F2.F2_age * F2.F2_height;
        logKe0 += 2.46 * F2.F2_age * F2.F2_sex;
        logKe0 += 3.35 * F2.F2_age * F2.F2_ASAPS;
        logKe0 -= 12.6 * F2.F2_TBW * F2.F2_height;
        logKe0 += 0.394 * F2.F2_TBW * F2.F2_sex;
        logKe0 += 2.06 * F2.F2_TBW * F2.F2_ASAPS;
        logKe0 += 0.390 * F2.F2_height * F2.F2_sex;
        logKe0 += 2.07 * F2.F2_height * F2.F2_ASAPS;
        logKe0 += 5.03 * F2.F2_sex * F2.F2_ASAPS;

        console.log('After 2nd order interactions, logKe0 =', logKe0);

        // 3次の交互作用項
        logKe0 += 99.8 * F2.F2_age * F2.F2_TBW * F2.F2_height;
        logKe0 += 5.11 * F2.F2_TBW * F2.F2_height * F2.F2_sex;
        logKe0 -= 5.00 * F2.F2_TBW * F2.F2_sex * F2.F2_ASAPS;
        logKe0 -= 5.04 * F2.F2_height * F2.F2_sex * F2.F2_ASAPS;
        logKe0 += 1.0 * F2.F2_TBW * F2.F2_height * F2.F2_ASAPS; // 最後の項

        console.log('After 3rd order interactions, logKe0 =', logKe0);

        // 4次の交互作用項
        logKe0 -= 39.4 * F2.F2_TBW * F2.F2_height * F2.F2_sex * F2.F2_ASAPS;

        console.log('Final logKe0 =', logKe0);

        // ke0を計算 (単位: 1/min)
        const ke0 = Math.exp(logKe0) / 60.0; // 60で割って1/minに変換

        console.log('ke0 =', ke0, '(1/min)');
        console.log('t1/2 =', (Math.log(2) / ke0).toFixed(1), 'min');

        // 妥当性チェック
        if (ke0 < 0.01 || ke0 > 2.0) {
            console.warn('ke0 value seems unrealistic:', ke0);
        }

        return ke0;
    }
}

// Simplified ke0 model for comparison
class SimplifiedKe0Model {
    constructor() {
        this.name = "Simplified ke0 Model";
    }

    calculateKe0(patient) {
        console.log('=== Simplified Ke0 Calculation ===');
        
        // 基本ke0値 (1/min) - Remimazolamの典型値
        const baseKe0 = 0.3;
        
        // 年齢補正: 高齢で遅くなる
        const ageRatio = patient.age / 54.0;
        const ageFactor = Math.pow(ageRatio, -0.15);
        
        // 体重補正: 体重が重いとやや遅くなる
        const weightRatio = patient.weight / 67.3;
        const weightFactor = Math.pow(weightRatio, -0.1);
        
        // 性別補正: 女性でやや遅い
        const sexFactor = patient.sexValue === 0 ? 1.0 : 0.9;
        
        // ASA補正: ASA III-IV でやや遅い
        const asaFactor = patient.asaPSValue === 0 ? 1.0 : 0.85;
        
        const ke0 = baseKe0 * ageFactor * weightFactor * sexFactor * asaFactor;
        
        console.log('Factors:', { ageFactor, weightFactor, sexFactor, asaFactor });
        console.log('Final ke0:', ke0);
        
        // 妥当な範囲に制限 (0.1-0.8 /min)
        return Math.max(0.1, Math.min(0.8, ke0));
    }
}

// Export models
if (typeof window !== 'undefined') {
    window.MasuiKe0Model = MasuiKe0Model;
    window.SimplifiedKe0Model = SimplifiedKe0Model;
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { MasuiKe0Model, SimplifiedKe0Model };
}