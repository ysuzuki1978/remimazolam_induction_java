# Remimazolam PK/PD Simulator (Induction)

レミマゾラム麻酔導入用の薬物動態・薬力学シミュレーター（導入期特化版）

## 概要

本アプリケーションは、レミマゾラムの麻酔導入期における薬物動態をリアルタイムでシミュレーションする教育・研究用ツールです。指示通りの**Masui PKモデル完全実装**により、数値解析による厳密なke₀計算と重回帰モデルによる近似計算の両方を提供します。

## 主要機能

- **リアルタイムシミュレーション**: 麻酔導入期の濃度変化を1秒間隔で更新
- **Masui完全実装**: 指示通りの2つのke₀計算方法を正確に実装
- **スナップショット機能**: 任意のタイミングで濃度値を記録・保存
- **患者情報編集**: リアルタイムでパラメータ変更可能
- **計算方法比較**: 数値解析・重回帰・従来手法の結果比較
- **レスポンシブUI**: あらゆるデバイスで最適表示

## 計算モデル

### Masui 2022 PK モデル

```
V1 = 3.57 × (ABW/67.3)
V2 = 11.3 × (ABW/67.3)  
V3 = (27.2 + 0.308×(AGE-54)) × (ABW/67.3)
CL = (1.03 + 0.146×SEX - 0.184×ASAPS) × (ABW/67.3)^0.75
Q2 = 1.10 × (ABW/67.3)^0.75
Q3 = 0.401 × (ABW/67.3)^0.75
```

### 体重計算

```
IBW = 45.4 + 0.89×(height-152.4) + 4.5×(1-sex)
ABW = IBW + 0.4×(TBW-IBW)
```

### Ke0 計算モデル - 完全実装

**方法A - 数値解析（厳密解）:**
1. 3次方程式解: x³ + a₂x² + a₁x + a₀ = 0 をCardanoの公式で解く
2. 係数計算: 血漿濃度式の指数関数係数A, B, Cを求める
3. 数値的根探索: f(ke₀) = 0をブレント法で解く（区間[0.15, 0.26]）
4. 最大効果到達時間: t_peak = 2.6分で最適化

**方法B - 重回帰モデル（近似解）:**
- 補助関数F(x): 年齢・体重・身長・性別・ASA-PSの多項式関数
- 補助変数F2(x): 中心化された変数による交互作用項
- 重回帰式: 25項からなる複雑な多変量回帰モデル
- 個別化計算: 患者特性に応じた高精度ke₀推定

**計算優先順位:** 数値解析 → 重回帰 → フォールバック（従来式）

### VHAC効果部位濃度計算

効果部位濃度の計算には、血漿濃度の変化パターンに応じて最適な手法を自動選択する**Variable-step Hybrid Algorithm for Ce (VHAC)**を採用：

**1. 定常状態（血漿濃度一定）:**
```
Ce(t) = Cp + (Ce₀ - Cp) × e^(-ke₀×Δt)
```
- 解析解による高精度計算
- 血漿濃度変化 < 1×10⁻⁶ で適用

**2. 線形変化（一次関数的変化）:**
```
一般解析解: Ce(t) = Cp(t) + (Ce₀ - Cp₀ + slope/ke₀) × e^(-ke₀×Δt) - slope/ke₀
```
- 線形補間 + 解析解の組み合わせ
- 通常の投与パターンで最も使用される手法

**3. 微小時間ステップ（ke₀×Δt < 0.001）:**
```
テイラー展開: Ce(t) = Ce₀ + Δt×ke₀×(Cp₀ - Ce₀) + (Δt²×ke₀×slope)/2
```
- 極小時間ステップ用の数値安定性確保
- 高頻度サンプリング時の精度向上

**手法選択ロジック:**
```javascript
if (|Cp_current - Cp_prev| < 1e-6) {
    // 定常状態 → 解析解
    Ce = Cp + (Ce_prev - Cp) * exp(-ke0 * dt);
} else if (|ke0 * dt| < 0.001) {
    // 微小ステップ → テイラー展開
    Ce = Ce_prev + dt * ke0 * (Cp_prev - Ce_prev) + dt² * ke0 * slope / 2;
} else {
    // 線形変化 → 一般解析解
    Ce = Cp + (Ce_prev - Cp_prev + slope/ke0) * exp(-ke0 * dt) - slope/ke0;
}
```

**VHAC法の利点:**
- **高精度**: 各状況に最適化された数学的解法
- **数値安定性**: 条件に応じた手法切り替えによる安定計算
- **計算効率**: 解析解優先による高速処理
- **適応性**: あらゆる投与パターンに対応

## ファイル構造

```
remimazolam_java_induction/
├── index.html                    # メインアプリケーション
├── debug.html                    # 開発・デバッグ用UI
├── test.html                     # 基本機能テスト用
├── ke0_comparison.html           # Ke0計算方法比較ツール
├── ke0_debug.html               # Ke0計算詳細デバッグ
├── masui-ke0-test.html          # Masui完全実装テストページ
├── package.json                  # プロジェクト設定・依存関係
├── README.md                     # このファイル
├── LSODA_IMPLEMENTATION.md       # LSODA積分法実装詳細
├── scripts/                      # JavaScriptファイル
│   ├── main.js                  # メインアプリケーションロジック
│   ├── main-simple.js           # 簡略版メインロジック（開発用）
│   ├── enhanced-main.js         # 拡張版メインロジック（開発用）
│   ├── masui-ke0-exact.js       # Masui完全実装計算エンジン
│   └── utils/                   # ユーティリティライブラリ
│       ├── lsoda.js             # LSODA積分法実装
│       ├── desolve.js           # 汎用ODE解法ライブラリ
│       └── ke0-masui.js         # Ke0計算比較用モデル
└── styles/                       # CSSファイル
    ├── main.css                 # メインスタイルシート
    └── components.css           # UIコンポーネントスタイル
```

### ファイル詳細説明

#### HTMLファイル

- **index.html**: 本番用メインアプリケーション。リアルタイムシミュレーション機能搭載
- **debug.html**: 開発段階で使用したデバッグ用UI。現在は使用停止、履歴保存目的で保持
- **test.html**: 基本機能のテスト用。開発初期段階で使用、現在は参考資料として保持
- **ke0_comparison.html**: 異なるKe0計算方法の比較ツール。開発時の検証用、現在も参考可能
- **ke0_debug.html**: Ke0計算の詳細デバッグページ。開発時の問題解析用、現在は使用停止
- **masui-ke0-test.html**: Masui完全実装の検証用テストページ。現在も動作確認に使用可能

#### JavaScriptファイル

**メインロジック:**
- **main.js**: 現在の本番メインロジック。リアルタイムシミュレーション実装
- **main-simple.js**: 開発初期の簡略版。基本機能のみ実装、現在は使用停止
- **enhanced-main.js**: 拡張機能開発版。追加機能テスト用、現在は使用停止

**計算エンジン:**
- **masui-ke0-exact.js**: 指示通りのMasui完全実装。数値解析と重回帰の両方を実装

**ユーティリティライブラリ:**
- **lsoda.js**: LSODA積分法の完全実装。血漿濃度計算に使用
- **desolve.js**: 汎用ODE解法ライブラリ。複数の積分手法を提供
- **ke0-masui.js**: Ke0計算方法比較用。開発時の検証目的、現在は比較ツールで使用

#### CSSファイル

- **main.css**: メインスタイルシート。アプリケーション全体のデザイン定義
- **components.css**: UIコンポーネント専用スタイル。ボタン、フォーム等の詳細デザイン

#### 設定・ドキュメント

- **package.json**: Node.js環境での依存関係管理（開発環境用）
- **LSODA_IMPLEMENTATION.md**: LSODA積分法の詳細実装解説
- **README.md**: プロジェクト全体の説明（このファイル）

## 技術仕様

- **言語**: Vanilla JavaScript (ES6+)
- **スタイル**: CSS3 (Grid, Flexbox)
- **計算手法**: **Masui完全実装エンジン**
  - **Ke0計算**: 数値解析（厳密解）+ 重回帰モデル（近似解）
  - **3次方程式**: Cardanoの公式による解析解
  - **数値的根探索**: ブレント法による高精度計算
- **血漿濃度計算**: LSODA積分法（Adams法・BDF法の自動切り替え）
- **効果部位濃度**: VHAC法 (Variable-step Hybrid Algorithm for Ce)
  - 定常状態: 解析解 Ce(t) = Cp + (Ce₀-Cp)×e^(-ke₀×Δt)
  - 線形変化: 一般解析解（線形補間+指数関数解）
  - 微小ステップ: テイラー展開による数値安定化
  - 判定条件: 血漿濃度変化量・時間ステップサイズで自動選択
- **更新頻度**: 1秒間隔
- **ブラウザ対応**: Chrome, Firefox, Safari, Edge (最新版)

## 開発者向け情報

### カスタマイズ

患者パラメータや計算モデルの変更は以下のファイルを編集してください：

- **PKパラメータ**: `scripts/masui-ke0-exact.js` の MASUI_THETA 定数
- **Ke0計算**: `scripts/masui-ke0-exact.js` の MasuiKe0Calculator クラス
- **UI設定**: `scripts/main.js` の定数部分

### 新機能追加

1. `scripts/main.js` にロジックを追加
2. `index.html` にUI要素を追加  
3. `styles/` にスタイルを追加

### 開発環境セットアップ

```bash
# 依存関係インストール（開発時のみ）
npm install

# ローカルサーバー起動
npx http-server -p 8080

# ブラウザでアクセス
open http://localhost:8080
```

## 免責事項

本アプリケーションは教育・研究用ツールです。実際の臨床判断の根拠として使用しないでください。すべての臨床判断は、資格を持つ医療専門家の責任において行われるべきです。

## ライセンス

MIT License

## 参考文献

1. Masui, K., et al. (2022). A population pharmacokinetic model of remimazolam for general anesthesia and consideration of remimazolam dose in clinical practice. Journal of Anesthesia, 36(4), 493-505.

2. Masui, K., & Hagihira, S. (2022). Equilibration rate constant, ke0, to determine effect-site concentration for the Masui remimazolam population pharmacokinetic model in general anesthesia patients. Journal of Anesthesia, 36(6), 733-742.

## 開発履歴

- **v1.0**: 基本機能実装（main-simple.js）
- **v1.1**: 拡張機能追加（enhanced-main.js）
- **v2.0**: Masui完全実装（masui-ke0-exact.js）
- **v2.1**: LSODA積分法統合（lsoda.js）
- **v3.0**: VHAC法統合（Variable-step Hybrid Algorithm for Ce）
- **v3.1**: 本番リリース（main.js + 完全実装）

---

**開発者**: YASUYUKI SUZUKI  
**所属**: 愛媛大学大学院医学系研究科薬理学  
**開発環境**: JavaScript/HTML/CSS  
**開発支援**: Claude Code (Anthropic)