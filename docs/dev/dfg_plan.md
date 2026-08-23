# DFG / OPA / OPO の位相整合 — 実装計画

> 状況（2026-08-23）: 段階 1–4 実装済み（v0.17.0）。`wl_idler`, `dk_dfg`, `pmFactor_dfg`,
> `qpm_period_dfg`, `pmAngles_dfg`, `tuning_dfg`, `d_dfg`, `deff_dfg`; アプリに Process 切替、
> チューニング曲線（Δk のゼロ等高線; 根探索は曲線描画には不要だった）、QPM 周期曲線、OPA 帯域、
> GVM 3 種。段階 5（NOPA）は未着手 — 理論ノートから。
> 文献値による検証は KTP x カット NCPM OPO の 1 点のみ（`docs/validation.md`）。Dmitriev / Myers 1995 /
> Cheng 1988 の PDF が手元に無く、他のチューニング曲線は自己無撞着テスト（SFG との厳密一致、
> 角度⇄波長の往復）で押さえている。PDF が入手できたら表に追加する。

## 0. 物理の整理 — 何が同じで何が違うか

三波混合 ω₃ = ω₁ + ω₂ の位相整合条件 **Δk = k₃ − k₁ − k₂ = 0 は SFG・DFG・OPA・OPO で同一**。
違うのは「何を与えて何を求めるか」だけ:

| 過程 | 与えるもの | 求めるもの | 帯域の固定条件 |
|---|---|---|---|
| SFG/SHG | λ₁, λ₂ | λ₃, θ_pm | λ₂ 固定で λ₁ を振る（λ₃ も動く） |
| DFG | λ_p(=λ₃), λ_s(=λ₁) | λ_i(=λ₂), θ_pm | 同上の読み替え |
| OPA | λ_p, θ（または T, Λ） | **λ_s, λ_i** | **λ_p 固定**で λ_s を振る → λ_i は逆向きに動く |
| OPO | λ_p, θ/T/Λ | λ_s(θ), λ_s(T), λ_s(Λ) **チューニング曲線** | 同上 |

したがって:

- **共線 DFG の角度計算は `pmAngles_sfg(λ_s, λ_i)` の薄いラッパ**で済む（段階 1）。
- 新しい計算は 2 つ: **(a) 波長についての根探索**（チューニング曲線、段階 2）と
  **(b) ポンプ固定・エネルギー保存付きの帯域**（段階 3）。SFG の「λ₁ のみ」帯域は
  λ₃ を動かすので OPA の帯域とは別物。OPA の信号帯域は Δk(λ_s; λ_i = f(λ_p, λ_s)) の幅で、
  展開の一次項が **信号と遊休の群速度差 1/v_s − 1/v_i** に比例する。退化点と NOPA で
  それが消えて広帯域になる — 超短パルス用途の本質。
- d_eff は Kleinman 対称のもとで三波の置換に不変、Miller 則の χ(ω₃)χ(ω₁)χ(ω₂) も対称。
  **`deff_dfg` は `deff_sfg` と同値**（段階 4 は API だけ）。QPM 周期も同じ Δk から出る。
- 添字の順番: ndispers は (1, 2, 3) = (低, 低, 高) 周波数。SNLO の QMIX も (red1, red2, blue)。
  OPA では **(1, 2, 3) = (signal, idler, pump)** と読む。'ooe' = Type I OPA（s, i 同偏光、
  pump が e）、'oee'/'eoe' = Type II。ここは docstring とアプリで明示する。

### スコープ外（明記する）

利得・閾値・変換効率・パルス伝搬（SNLO の 2D-mix / PW-OPO 相当）は扱わない。
ndispers は分散と位相整合のパッケージで、強度を持たない。将来やるとしても別パッケージ。

---

## 1. 共線 DFG API（ラッパ）

`_baseclass.Medium` に追加。すべて `wl_p, wl_s` を取り、内部で `wl_i = 1/(1/wl_p − 1/wl_s)` を
作って `_sfg` 版を `(wl_s, wl_i, …)` で呼ぶ。

```python
def wl_idler(self, wl_p, wl_s): ...              # 1/λi = 1/λp − 1/λs、λs ≤ λp なら ValueError
def dk_dfg(self, wl_p, wl_s, angle, T, pol_s, pol_i, pol_p)
def pmAngles_dfg(self, wl_p, wl_s, T, tol_deg=0.001, deg=False)   # 返り値 dict に 'wl_i' を追加
def pmFactor_dfg(self, wl_p, wl_s, angle, T, pol_s, pol_i, pol_p, L_mm)
def qpm_period_dfg(self, wl_p, wl_s, angle, T, pol_s, pol_i, pol_p, order=1)
```

`NonlinearGroup` に `deff_dfg(wl_p, wl_s, theta, phi, T, pol_s, pol_i, pol_p)` と
`d_dfg(il, wl_p, wl_s, T)`（同値の別名。Miller の基準が (wl1, wl2) で書かれている
ことを docstring に書く）。

偏光引数の並びは **(pol_s, pol_i, pol_p)** で SFG の (pol1, pol2, pol3) と同じ位置。
辞書のキー 'ooe' 等も同じ意味になるので `TYPE_OF` が再利用できる。

テスト: `pmAngles_dfg(λp, λs) == pmAngles_sfg(λs, λi)` を全結晶で厳密一致（ラッパの証明）。
`wl_idler` の境界: λs = λp で ZeroDivision ではなく ValueError。

---

## 2. チューニング曲線（新規の根探索）

```python
def tuning_dfg(self, wl_p, angle_rad, T_degC, pol_s, pol_i, pol_p,
               wl_s_range=None, tol=1e-6) -> list[(wl_s, wl_i)]
```

固定 (θ, T) で Δk(λ_s) = 0 の λ_s を返す。`pmAngles_sfg` と同じ方針:
**粗格子で符号反転を探して `helper.brentq` で詰める**。根は複数あり得る
（Type II の 2 分枝、角度チューニングの折り返し）ので全部返す。

- 探索区間: λ_s ∈ (2λ_p, λ_s,max)。2λ_p は退化点（λ_s = λ_i）。信号と遊休は入れ替え対称なので
  **λ_s ≥ λ_i の側だけ**探して重複を避ける。上限は docstring の validity range か
  引数で上書き。範囲外へ外挿していることは `pmAngles` と同じ扱い（返すが警告しない；
  アプリ側で validity を表示）。
- 退化点そのものは端点解。`pmAngles_sfg` に入れた端点・傾き判定と同じ論理を使う。
- 格子の密度: 波長では Δk が急に変わる領域（退化近傍、吸収端近く）があるので、λ_s ではなく
  **ν_s = 1/λ_s（等間隔）**で格子を切る。ν_i = ν_p − ν_s が線形になるので素直。

これを使って 3 種類の曲線をアプリ側で描く: λ_s(θ)（角度チューニング）、λ_s(T)（温度チューニング、
NCPM 結晶と PPLN）、λ_s(Λ)（QPM: `qpm_period_dfg` の逆問題。Λ は `tuning` と同じ根探索で
Δk − 2π/Λ = 0 を解く。`qpm_period_dfg` を使うなら順方向で λ_s 格子 → Λ を出して描くだけ。
**順方向で十分**: 描くのは曲線で、逆問題の根は要らない）。

テスト（往復不変条件）: 任意の結晶・(λp, λs) で θ_pm = `pmAngles_dfg` → `tuning_dfg(λp, θ_pm)` が
λs（または λi）を 1e-6 µm 以内で返す。

---

## 3. OPA 帯域・群速度ミスマッチ

アプリ側（`acceptance()` は既にある）で固定条件を変えるだけ:

| 行 | 動かすもの | 固定 |
|---|---|---|
| 角度 / 温度 | θ / T | λ_p, λ_s（λ_i も） |
| **信号（λ_i 自由）** | λ_s | λ_p（λ_i = f(λ_p, λ_s) が追従） |
| **ポンプ（λ_s 固定）** | λ_p | λ_s（λ_i 追従） |

「信号帯域」が OPA の本命で、SNLO の "mix acceptance" とは違う量になることを画面で注記。
レポートには **GVM 3 種** 1/v_s − 1/v_i、1/v_p − 1/v_s、1/v_p − 1/v_i と、
L での walk-through 時間を出す（超短パルス OPA の設計量）。

退化点 (λ_s = λ_i = 2λ_p, Type I) の自動検出: 退化では 1/v_s − 1/v_i = 0 で帯域は二次項で決まる。
格子の「符号反転」法は退化で Δk が接するだけで横切らないケースがあるので、`acceptance()` の
囲い込みは sinc² の最初の零点を探す現行方式のままで良い（Δk ではなく sinc² を見ている）。

---

## 4. アプリ

`phasematching.py` に **Process** セレクタ: `SFG / SHG` | `DFG / OPA / OPO`。

DFG モードの入力: λ_p, λ_s（または λ_i — ラジオで切替）, T, L。表示:
- 位相整合角の表（6 通り、`TYPE_OF` 流用、見出しは "(s, i, p)"）
- 選択解での n / n_g / v_g / GVD / walk-off の表（SFG と同じセル、ラベルだけ s/i/p）
- **チューニング曲線**: λ_s, λ_i vs θ（または φ）を 0–90° で、現在の解に縦線。
  T をスライダで動かせば温度チューニングが見える。QPM 結晶（SLN/MgOLN/SLT/KTP）では
  Λ vs λ_s も。
- 帯域: §3 の 4 行
- d_eff: 既存セルをそのまま（`deff_dfg` = `deff_sfg`）
- レポート: 既存 REPORT に PROCESS 行と GVM 3 種を追加

Explorer はそのまま。

---

## 5. 非共線（NOPA）— 独立した段階

ベクトル位相整合 **k_p = k_s + k_i**。ポンプ方向 θ_p、信号–ポンプ内部角 α を与えると
遊休方向 β は閉じる条件から決まり、未知数 (θ_p, β) に対し Δk_∥ = Δk_⊥ = 0 の 2 式。
各波の n は各波の伝播角で評価する — `n(wl, theta, T, pol)` が波ごとに角度を取るので
**既存の n で書ける**（主平面内の非共線に限る。面外は φ も動くので後回し）。

```python
def dk_ncpm(self, wl_p, wl_s, theta_p, alpha, T, pols)  -> (dk_par, dk_perp)
def pmAngles_nopa(self, wl_p, wl_s, alpha, T) -> θ_p の解
def alpha_broadband(self, wl_p, wl_s, T, pols)  # v_s = v_i cos(β−α) 条件（群速度整合）
```

広帯域条件 v_g,s = v_g,i·cos(Ω)（Ω = 信号–遊休角）は解析式で出せるので、α_opt を直接返し、
アプリではその α で θ_p と信号帯域を表示する。検証値: **BBO, 400 nm ポンプ, Type I の
NOPA: θ ≈ 31.5°, α ≈ 3.7°（内部角; Cerullo & De Silvestri 2003 / Wilhelm, Piel, Riedle 1997）**
— 文献 PDF で数値と「内部/外部」を確定してから実装（記憶から転写しない）。

ここは理論ノート（`deff_theory.tex` とは別に `ncpm_theory.tex`）を先に書く価値がある:
ベクトル Δk の定義、歩行角の扱い（非共線ではポインティングではなく k で閉じる）、
広帯域条件の導出、Type II での分枝。d_eff も各波の偏光ベクトルが別方向を向くので、
**既存の `deff_sfg` の縮約は各波の k̂ を別々に受ければそのまま使える**（`_base.py` の
偏光ベクトル生成を wave ごとの角度に一般化）。

---

## 6. 段階と順序

| # | 内容 | 新規コード | 検証 |
|---|---|---|---|
| 0 | `docs/dev/dfg_plan.md`（本書）、`deff_theory.tex` に「三波の置換対称性」の 1 節追記 | 文書 | — |
| 1 | §1 ラッパ + `wl_idler` + テスト | `_baseclass.py` ~60 行、`groups/_base.py` ~15 行 | SFG との厳密一致 |
| 2 | §2 `tuning_dfg` + 往復テスト | `_baseclass.py` ~60 行 | 往復不変条件、文献チューニング曲線 |
| 3 | §3–4 アプリ DFG モード、帯域、GVM、レポート | `apps/phasematching.py` | SNLO / 文献値との照合を `docs/validation.md` に追記 |
| 4 | リリース v0.17.0 | — | — |
| 5 | §5 NOPA: 理論ノート → 実装 → アプリの α スライダ | 別計画 | BBO NOPA 文献値 |

段階 1–3 は一つのリリース。段階 5 は理論ノートのレビューを挟むので別リリース。

### 検証に使う文献値（手元 PDF で確認してから数値を固定する）

- BBO Type I OPO/OPA, 355 nm ポンプの角度チューニング曲線（Dmitriev ハンドブック、または Cheng/Bosenberg/Tang 1988）
- KTP Type II OPO, 1064 nm ポンプ, xz 面のチューニング（Dmitriev; SNLO）
- LBO Type I, 355 nm ポンプ, xy 面 / 温度チューニング NCPM（Dmitriev; Castech）
- PPLN OPO, 1064 nm ポンプ, Λ = 28–31 µm の温度/周期チューニング（Myers et al. JOSA B 12, 2102 (1995)）— SLN/MgOLN の QPM の検証。`SLN` は温度依存を持つのでこれが効く
- BBO NOPA（段階 5）: Cerullo & De Silvestri, Rev. Sci. Instrum. 74, 1 (2003); Wilhelm, Piel, Riedle, Opt. Lett. 22, 1494 (1997)

各値は `docs/validation.md` の表に「計算値 / 文献値 / 出典」で 1 行ずつ追加する。

### 決めておくこと

- 命名: `_dfg` で統一し、OPA/OPO は docstring とアプリの見出しに出す（物理過程名として DFG が
  SFG の対；API が 2 系統に割れない）。`_opa` の別名は作らない。
- `pmAngles_dfg` の返り値: `pmAngles_sfg` と同形 + `'wl_i'`。既存の `'wl3'` に相当。
- λ_s > λ_p や λ_s = λ_p は ValueError（遊休が定義できない）。アプリは入力欄で弾く。
