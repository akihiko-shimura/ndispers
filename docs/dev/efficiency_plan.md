# SFG 変換効率（半解析的）— 実装計画

> 状況（2026-08-30）: 計画のみ。未着手。
>
> `dfg_plan.md` の「スコープ外」節は変換効率を別パッケージ扱いとしていたが、
> 本計画はそれを改訂する: **半解析的な簡易見積もりに限り** ndispers に入れる。
> スプリットステップ数値伝搬（SSFT）は企業内 IP として外に留める — 入れるのは
> 教科書公式（Armstrong 1962; Guha & Gonzalez 2014, eq. 4.118）の実装だけで、
> 原型は社内の自作スクリプト（2021, 自著）にある。
>
> **理論ノート**: `efficiency_theory.tex/pdf`（実装前レビュー用）。導出の過程で
> 原型スクリプトの a₁ における d↔d₁ の対応交差の疑義を発見（ノート §3.A）。
> 縮退運用では不可視。実装はノートの式を正とし、Guha 本 eq. (4.118) との
> 照合と波交換対称性テストで確定させる。
> SHG の正規化 I_a^SHG は Boyd 4th ed. §2.7（原典確認済み、2026-08-30）の
> Δk=0 解で確定 — SFG の I_a の縮退代入と同形、追加因子なし。
> d↔d₁ の対応は Armstrong (1962) 原論文 §VI（原典確認済み、同日）の
> eqs. (6.6)/(6.10) で決着: ノートの式が正、原型スクリプトは転置。
> **理論のレビュー事項はゼロ**。著者レビュー後、段階 1 の実装に入れる。

## 0. 物理と近似 — 何を計算し、何を無視するか

パルス SFG のエネルギー変換効率 η_e。Jacobi 楕円関数 sn による
Armstrong 解を、ガウス型の空間・時間プロファイルにわたり数値積分する
（Guha & Gonzalez eq. 4.118）。

含む: **ポンプ枯渇**（低変換近似ではない）、位相不整合 σ = Δk·L、
2 入力波の任意のビーム径比 ρ・パルス幅比 τ・強度比 p₁。

無視する（docstring に明記し、§3 の自己診断警告で守る）:

| 無視する効果 | 妥当条件 | ndispers 自身で検査可能か |
|---|---|---|
| 空間ウォークオフ | ρ_wo·L ≪ r₀ | `woa_theta` × L で可 |
| 時間ウォークオフ (GVM) | GVM·L ≪ t₀ | `GVM` × L で可 |
| 回折 | L ≪ z_R | n と r₀ から可 |
| GVD によるパルス変形 | GVD·L ≪ t₀² | `GVD` × L で可 |

この「モデルが自分の適用限界を検査できる」ことが SNLO 的な電卓に対する
差別化点で、ValidityWarning / TemperatureWarning の文化とも一貫する。

添字の非対称性に注意: 無次元化は波 1 を基準に行う（d = λ₁/λ₃, K₁ = I₁₀/I_a）。
波 1 と 2 の入れ替えで結果は不変のはずなので、これはテストで確認する（§4）。

## 1. API 案（2026-08-30 改訂 2 — Beam dataclass 方式で確定）

ビーム・パルスの諸元は frozen dataclass に束ね、メソッドは 2 本にする。
pulse か CW かは**ビームの型**で決まる（設定をオブジェクトに束ねる思想は
SSFT パッケージと同じだが、オブジェクト互換までは意図しない）。

```python
from dataclasses import dataclass

@dataclass(frozen=True)
class PulsedBeam:
    wl_um: float    # 波長 (µm)
    E_uJ: float     # パルスエネルギー (µJ)
    w_um: float     # 1/e² 強度半径 = ガウシアン w0 (µm), I ∝ exp(−2r²/w²)
    t_fs: float     # FWHM パルス幅 (fs)
    pol: str        # 'o' | 'e'

@dataclass(frozen=True)
class CWBeam:
    wl_um: float
    P_W: float      # パワー (W)
    w_um: float
    pol: str
```

置き場所は `ndispers/_efficiency.py`、`ndispers.PulsedBeam` /
`ndispers.CWBeam` としてトップレベルにエクスポート。`__post_init__` で
正値・pol の検証（教育的 ValueError）。

```python
eta_sfg(beam1, beam2, theta_rad, phi_rad, T_degC, pol3, L_mm,
        deff_pmV=None, dk_offset=0.0, n_grid=50, details=False)

eta_shg(beam, theta_rad, phi_rad, T_degC, pol3, L_mm,
        deff_pmV=None, dk_offset=0.0, n_grid=50, details=False)
```

- 結晶メソッド側に残るのは**結晶の性質だけ**: 角度対 (theta, phi)、温度、
  和周波の偏光 pol3、結晶長。ビームの物理は全部 Beam に。
- `eta_sfg` は 2 ビームの型一致を検査。PulsedBeam と CWBeam の混在は
  教育的 TypeError（「両方 PulsedBeam か両方 CWBeam に」）。
- 型で分岐: PulsedBeam 対 → (r, t) 2 次元積分、CWBeam 対 → r の 1 次元
  積分。解の骨格は共有し、時間ガウシアンの有無だけを変える。

### 引数の設計判断（改訂前から維持）

- **角度は (theta_rad, phi_rad) の対** — `deff_sfg` と同形。n と Δk は
  位相整合角だけで決まるが、d_eff は方位角にも依存する（BBO の d22 の
  φ 依存。原型スクリプトが φ = π/4 を別途持っていた理由）。主平面クラス
  では面が固定する方に None。`deff_pmV` 明示時は φ は None でよい。
- **単位は実験室仕様を直接受ける**: 換算（r₀ = w/√2,
  t₀ = t_FWHM/(2√ln2)）は内部で行い docstring に明記。
- **f_rep は取らない**。E = P_avg/f_rep はユーザーの 1 行。
- **`dk_offset`**（rad/µm）: 残差不整合の指定用。
- **`n_grid`**: 積分メッシュ。`tuning_dfg` と同名規約。
- pol1/pol2 は Beam の `pol` に移った。`pol3` だけがメソッド引数。

### SHG の扱い（物理の要点）

`eta_shg` は `eta_sfg` に同じビームを 2 度渡す糖衣**ではない**。

- **Type I（beam.pol == 両基本波の偏光）**: 単一場の Armstrong SHG 解
  （Δk=0 で tanh²、一般 Δk で sn）を単一ガウシアンで積分。縮退因子
  （d_eff の定義と結合方程式の因子 2）はここで一度だけ正しく処理し、
  規約を docstring に明記。
- **Type II**: 1 本のビームの o/e 分割なので、内部で `eta_sfg` に委譲
  （E/2 ずつの PulsedBeam を 2 つ合成）。呼び方: `eta_shg` に
  `pol_fund=('o','e')` のような指定を足すか、Type II はユーザーが直接
  `eta_sfg` を呼ぶ規約にするか — **後者**（Type II SHG は縮退 SFG として
  `eta_sfg` で書けるので、`eta_shg` は Type I 専用とし docstring で誘導）。
  API を増やさない。
- テスト: 波長 1 ppm ずらした Type II 縮退 `eta_sfg` と Type I `eta_shg`
  の連続性ではなく、Type I は tanh² 平面波極限と、Type II 経路は
  `eta_sfg` の対称性テストで別々に縛る。

### 返り値

- `details=False`（既定）: float。**η ≡ E₃/(E₁+E₂)**（CW は P₃/(P₁+P₂)、
  SHG は E_SH/E_fund）。
- `details=True`: dict —

```python
{'eta': float,
 'eta_photon_1': float,     # 波 1 の光子変換率（枯渇率）
 'eta_photon_2': float,     # 波 2 の同上（SHG では省略）
 'E3_uJ' | 'P3_W': float,
 'I1_peak_MWcm2': float, 'I2_peak_MWcm2': float,
 'deff_pmV': float, 'dk_radum': float, 'K1': float,
 'model_ratios': {'walkoff': x, 'gvm': x, 'diffraction': x, 'gvd': x}}
```

  強ポンプ＋弱シグナルのアップコンバージョンでは η より
  `eta_photon_2`（シグナル量子効率）が本命なので両方返す。
  `model_ratios` は §3 の自己診断の生値。

### スコープ外の確認

DFG/OPA/OPO の利得・効率（別解）、集光ビーム（Boyd–Kleinman）、
繰り返し・熱効果。Boyd–Kleinman は将来の自然な拡張だが含めない。

## 2. 依存関係 — scipy をどうするか

必要なのは `scipy.special.ellipj`（sn）と 2 次元 Simpson 積分。

- Simpson は numpy で 10 行 — 問題ない。
- **ellipj は自前実装しない**。特殊関数の転写ミスはこのリポジトリが最も
  恐れる種類のバグで、AGM/Landen 実装の検証コストが利得に見合わない。
- 結論: **optional extra `ndispers[eff]` で scipy** を入れる。
  `pip install ndispers` の「依存は numpy のみ」という Feature は守られる。
  import は遅延させ、無ければ教育的 ImportError
  （"pip install ndispers[eff]"）を出す — シグネチャエラーと同じ流儀。

## 3. 自己診断警告

§0 の表の 4 条件を計算し、比が閾値（暫定 0.3）を超えたら
`ModelValidityWarning`（新設、UserWarning 系）で「無視している効果 X が
このパラメータでは小さくない（比 Y）。η は過大評価になる」と 1 回警告。
閾値と比は details dict にも入れる。数値パッケージ（SSFT）への誘導は
書かない（IP 境界）。

## 4. 数値と検証

数値の注意:

- 積分メッシュ（原型: [0,5]×[0,5], 50×50）はデフォルトを収束テストで決め、
  `n_grid` 引数として公開。
- 楕円パラメータ m = b₁/b₂ → 1（完全整合・強変換）で ellipj は tanh に
  退化する。丸めで m > 1 に出たら clip。K₁ → 0 の縮退（b₂ → σ²/4 など）も
  踏む。エッジはテストで踏み抜く。

検証（`docs/validation.md` に行を足す）:

1. **原型スクリプトの再現**: 縮退ケースの実行値と 1e-3 相対で一致
   （実施済み。パラメータと数値は社内情報のためリポジトリには残さない）。
2. **無枯渇極限**: K₁ ≪ 1 で解析式 η ∝ K₁·sinc²(σ/2) 系に収束。
3. **Manley-Rowe 上限**: 任意パラメータで η_e ≤ d（波 1 の光子が全て
   変換されたときのエネルギー比）を超えない。ランダムパラメータで確認。
4. **1↔2 対称性**: 波 1 と 2 の入れ替え（と p₁, ρ, τ の逆数化）で η 不変。
5. **SHG 退化**: λ₁ = λ₂, 同一ビームで Boyd の tanh² 公式（平面波・
   パルス平均なし極限）と整合。
6. 文献値: Guha & Gonzalez の作例の数値が本文にあれば 1 点照合
   （書籍 PDF が手元にあるか要確認）。

## 5. ドキュメント

- conventions.md「Second-order nonlinearity」の下に効率の節
  （定義・単位・無視効果の表）。
- llms.txt に 2 行（シグネチャと「無視効果は警告される」）。llms-full 再生成。
- ノートブック §7 として CLBO の作例（原型スクリプトのパラメータ）。
- 位相整合アプリへの搭載は当面見送り（入力欄が 10 個増える。需要が出たら）。
- agent-evals.md に課題 4「CLBO 532→266 の効率見積もり」を追加し、
  コールドスタート試行を 1 回実施。

## 6. 段階

| 段階 | 内容 | 成果物 | 検証 |
|---|---|---|---|
| 1 | `_efficiency.py` 移植 + `NonlinearGroup.eta_sfg_pulse` + extra `[eff]` | 動く η | 検証 1, 2 |
| 2 | エッジケース（m→1, K₁→0, σ 大）と `n_grid` 収束 | 既定メッシュ確定 | 検証 3, 4, 5 |
| 3 | 自己診断警告 `ModelValidityWarning` | §3 | 閾値の妥当性を GVM の実例で確認 |
| 4 | ドキュメント一式 + validation.md + agent-eval 1 回 | §5 | ドリフトテスト通過 |
| 5 | リリース v0.20.0（新機能なので minor） | — | RELEASING.md の手順 |

見積もり: 段階 1–2 で計算コア、3–4 が同量。原型が動いているので
物理の実装リスクは低く、主リスクは単位系の取り違え（→ §1 の換算を
docstring とテスト両方に書く）と検証 5 の極限対応の導出。
