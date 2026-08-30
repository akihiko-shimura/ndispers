<img src="https://raw.githubusercontent.com/akihiko-shimura/ndispers/main/docs/assets/logo.svg" alt="" width="132" align="right">

# ndispers

[English](https://github.com/akihiko-shimura/ndispers/blob/main/README.md) | **日本語**

[![PyPI](https://img.shields.io/pypi/v/ndispers)](https://pypi.org/project/ndispers/)
[![test](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml/badge.svg)](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml)
[![apps](https://github.com/akihiko-shimura/ndispers/actions/workflows/pages.yml/badge.svg)](https://akihiko-shimura.github.io/ndispers/)

*ndispers* は、非線形光学・超高速光学で用いられる結晶とガラスの屈折率分散を
計算する Python パッケージです。文献で報告された Sellmeier 方程式と
熱光学係数を実装しています。

計算できる量:

- 屈折率
- 群遅延、群速度、群屈折率
- 群速度分散、3 次分散
- 2 波間の群速度不整合 (GVM)
- ウォークオフ角
- dn/dT、d²n/dT²

変数:

1. 波長
2. 誘電主軸に対する波数ベクトルの極角 (θ) または方位角 (φ)
3. 温度
4. 偏光（常光線・異常光線）

非中心対称結晶ではさらに、位相不整合 Δk、位相整合角、位相整合因子
sinc²(Δk·L/2)、許容幅（スペクトル・角度・温度）、擬似位相整合周期、
OPO チューニング曲線、および Miller 則で使用波長にスケーリングした
実効非線形係数 d_eff を計算します。

## なぜ ndispers か

この種の計算の既存ツールはアプリケーションです: SNLO
([as-photonics.com](https://as-photonics.com/products/snlo/)、Windows)、
[refractiveindex.info](https://refractiveindex.info/) と
[toolbox.lightcon.com](http://toolbox.lightcon.com/)（Web）、
[iPhasematch](https://apps.apple.com/jp/app/iphasematch/id492370060)（iOS）。
アプリケーションは入力された問いに答えますが、プログラムから呼び出すことは
できません。

ndispers は Python ライブラリです。SNLO の Qmix 機能が報告する量 —
屈折率、群速度と分散、ウォークオフ、位相整合条件、許容幅、実効非線形係数 —
を、numpy 配列を受けて返す関数として計算します。したがって波長・角度・温度の
グリッド上で評価し、より大きな計算の内部で使えます。すべての係数は出典まで
遡れます: 各クラスの docstring が、論文、有効範囲、および Sellmeier 方程式・
熱光学係数・d テンソルそれぞれの出典を記載しています。クラス名とパッケージの
バージョンを指定すれば、計算は完全に特定されます。

ライブラリが内省可能でデータの出典が明示されているため、AI アシスタントが
これを操作できます: 媒質は `dir()` または `catalog()` で列挙でき、docstring
が必要なメタデータを持ち、[llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt)
が API を 1 ページに記述しています。

対象範囲は SNLO より広く、線形の複屈折結晶（α-BBO、方解石、水晶、MgF₂、
YVO₄、サファイア）と等方媒質（ガラス、CaF₂、Si、Ge など）を同一の
インターフェースで含みます。

**特長**

- 全メソッドが numpy 配列を受けてブロードキャストします。評価される関数は
  事前生成されており、実行時依存は numpy のみです。medium オブジェクトは
  picklable で `multiprocessing` と併用できます。
- 温度微分 dn/dT・d²n/dT² は屈折率式の記号微分から得ており、数値差分では
  ありません。
- 位相整合: Δk、位相整合角、sinc² 因子、許容幅、QPM 周期、OPO チューニング
  曲線。全非中心対称結晶で Miller 波長スケーリングつきの d_eff。
- 同一結晶の複数のパラメータ化は別クラスです（LBO は 5 つ）。出典間の
  不一致は、1 つのデフォルトの背後に隠さず、見えるようにしています。
- 係数は引用論文から転写し、独立の再抽出で照合しています。計算値は回帰
  テストで固定されています（[検証記録](https://ndispers.readthedocs.io/en/latest/validation/)）。
- 結晶の追加は 1 ファイルです: 点群の基底クラス＋ Sellmeier と d 係数
  （[AGENTS.md](AGENTS.md)）。要望・貢献は
  [GitHub](https://github.com/akihiko-shimura/ndispers/issues) へ。
- インストールは任意のプラットフォームで `pip install ndispers`、
  Python ≥ 3.10 です。

## インストール

```
pip install ndispers
```

Python 3.10 以降。依存は numpy のみです。分散式は sympy の式として書かれて
いますが、評価される関数は事前に生成して同梱しています。式そのもの
（`n_expr` など）を扱う場合や、媒質を自分でサブクラス化して評価する場合は
`pip install ndispers[sym]` で sympy を追加してください。

## クイックスタート

[チュートリアルノートブック](examples/basic_usage.ipynb)が、Sellmeier 方程式の
確認から位相整合曲線までの手順をプロットつきで示しています。要点:

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
>>> bbo.n(0.532, 0, 25, pol='o')
1.674884049110459
>>> bbo.n(0.532, 3.1416/2, 25, pol='e')
1.5554658787539917
```

4 つの引数は

1. 波長 (µm)
2. θ 角 (rad)
3. 温度 (°C)
4. 偏光（`pol='o'` または `'e'`。デフォルト `'o'`）

θ = 0（光学軸に沿った伝搬）では常光線と異常光線の屈折率は一致します。
スカラー入力は float を返し、numpy 配列入力はブロードキャストして配列を
返します:

```python
>>> import numpy as np
>>> bbo.n(np.arange(0.2, 1.5, 0.2), 0, 25, pol='o')
array([1.89001202, 1.69328828, 1.66985875, 1.6612891 , 1.65633946,
       1.65252664, 1.64903624])
```

他の分散メソッド（`GD`, `GV`, `ng`, `GVD`, `TOD`, `GVM`, `woa_theta`,
`woa_phi`, `dndT`, `dndT2`）も同じ引数を取ります:

```python
>>> bbo.GVD(0.8, 0, 25, pol='o')   # fs^2/mm
71.86403019943364
```

光学的に等方な媒質 — ガラスと立方晶系の結晶 — は波長と温度のみを取ります:

```python
>>> silica = nd.media.glasses.FusedSilica()
>>> silica.n(1.064, 20)
1.4495857898590634
```

`help(bbo)`（IPython/Jupyter では `bbo?`）が材料データ、Sellmeier 方程式、
参考文献を表示します。`bbo.constants` は係数値を返します。

## 位相整合

`pmAngles_sfg` が和周波発生の位相整合角を解きます。β-BBO の 1064 nm SHG、
25 °C の場合:

```python
>>> bbo.pmAngles_sfg(1.064, 1.064, 25, deg=True)
{'wl3': 0.532,
 'ooe': {'theta': [22.884169498625802], 'phi': None},   # タイプ I
 'eeo': {'theta': [], 'phi': None},
 'oee': {'theta': [32.56045545648089], 'phi': None},    # タイプ II
 'eoe': {'theta': [32.56045545648089], 'phi': None},
 'eoo': {'theta': [], 'phi': None},
 'oeo': {'theta': [], 'phi': None}}
```

`dk_sfg` と `pmFactor_sfg` は任意角度での位相不整合と sinc² 因子を、
`acceptance_sfg` は sinc² 曲線の波長・角度・温度方向の FWHM（許容幅）を、
`qpm_period_sfg` は擬似位相整合周期を与えます。5% MgO:LiNbO₃ の x 軸伝搬・
全波異常光線 (d₃₃) での 1064 nm SHG では、
`MgOLN_Zelmon1997().qpm_period_sfg(1.064, 1.064, np.pi/2, 20, 'e', 'e', 'e')`
が 7.0 µm を返します。

差周波発生と光パラメトリック過程は、同じ相互作用にポンプ側から入ります:
`pmAngles_dfg(wl_p, wl_s, T)` が位相整合角とアイドラー波長を、
`tuning_dfg(wl_p, angle, T, pol_s, pol_i, pol_p)` が固定角度で位相整合する
シグナル/アイドラー対 — OPO チューニング曲線の 1 点 — を与え、
`qpm_period=` を渡せば周期分極反転結晶の温度チューニングになります。
`dk_dfg`, `pmFactor_dfg`, `qpm_period_dfg`, `deff_dfg` は SFG 版と対応し、
波の順序は（シグナル, アイドラー, ポンプ）です。x カット KTP の 1064 nm
ポンプでは `KTP_zx().tuning_dfg(1.064, np.pi/2, 25, 'o', 'e', 'o')` が
ノンクリティカル対 1571/3298 nm を返します。

位相整合角における実効非線形係数（φ = 90° カット）:

```python
>>> bbo.deff_sfg(1.064, 1.064, np.radians(22.88), np.radians(90), 25, 'o', 'o', 'e')
-1.9937...  # pm/V。全体の符号は規約による。d22 = 2.2 pm/V（1.064 µm SHG,
            # Shoji et al. 1999）、ウォークオフ込み
>>> bbo.d_sfg("d22", 0.8, 0.8, 25)
2.3300...   # Miller 則による 800 nm SHG での d22
```

符号と角度の規約は
[Conventions](https://ndispers.readthedocs.io/en/latest/conventions/)
ページに定義されています。各結晶の `_d_note` が d 係数の出典と、その結晶で
Miller スケーリングがどこまで検証されているかを記載しています。

## 収録媒質

各媒質は `nd.media.crystals` または `nd.media.glasses` のクラスです。
[媒質カタログ](https://ndispers.readthedocs.io/en/latest/api/crystals/)が
各クラスの Sellmeier 方程式・有効範囲・参考文献を一覧し、`nd.catalog()` が
同じ情報をデータとして返します。複数のパラメータ化がある結晶では、出典ごとに
別クラスになっています。

### 非線形光学結晶

非中心対称。位相整合、許容幅、d_eff が利用できます。

| 材料 | 略称 | 組成式 | 点群 | 光学的分類 |
|---|---|---|---|---|
| β-ホウ酸バリウム | β-BBO | β-BaB₂O₄ | 3m | 負の一軸性 |
| トリホウ酸リチウム | LBO | LiB₃O₅ | mm2 | 二軸性、3 主平面 |
| リン酸チタニルカリウム | KTP | KTiOPO₄ | mm2 | 二軸性、3 主平面 |
| トリホウ酸ビスマス | BiBO | BiB₃O₆ | 2 | 二軸性（単斜晶）、3 主平面 |
| ホウ酸セシウムリチウム | CLBO | CsLiB₆O₁₀ | 4̄2m | 負の一軸性 |
| リン酸二水素カリウム | KDP | KH₂PO₄ | 4̄2m | 負の一軸性 |
| 重水素化リン酸二水素カリウム | DKDP, KD*P | KD₂PO₄ | 4̄2m | 負の一軸性 |
| フルオロホウ酸ベリリウムカリウム | KBBF | KBe₂BO₃F₂ | 32 | 負の一軸性 |
| フルオロホウ酸ベリリウムルビジウム | RBBF | RbBe₂BO₃F₂ | 32 | 負の一軸性 |
| 四ホウ酸リチウム | LB4 (LBT) | Li₂B₄O₇ | 4mm | 負の一軸性 |
| ヨウ素酸リチウム | — | LiIO₃ | 6 | 負の一軸性 |
| リン化ゲルマニウム亜鉛 | ZGP | ZnGeP₂ | 4̄2m | 正の一軸性、中赤外 |
| チオガリウム酸銀 | AGS | AgGaS₂ | 4̄2m | 負の一軸性、中赤外 |
| セレン化ガリウム銀 | AGSe | AgGaSe₂ | 4̄2m | 負の一軸性、中赤外 |
| α-水晶 | — | SiO₂ | 32 | 正の一軸性 |
| ニオブ酸リチウム（5% MgO 添加コングルエント） | MgO:LN | MgO:LiNbO₃ | 3m | 負の一軸性、両光線 |
| ニオブ酸リチウム（1% MgO 添加ストイキオメトリック） | MgO:SLN | MgO:LiNbO₃ | 3m | 負の一軸性、異常光線のみ |
| タンタル酸リチウム（1% MgO 添加ストイキオメトリック） | MgO:SLT | MgO:LiTaO₃ | 3m | 負の一軸性 |

### 複屈折光学結晶

中心対称のため 2 次非線形性はありません。窓材・偏光子・補償板向けの分散、
ウォークオフ、熱光学を計算できます。

| 材料 | 略称 | 組成式 | 点群 | 光学的分類 |
|---|---|---|---|---|
| α-ホウ酸バリウム | α-BBO | α-BaB₂O₄ | 3̄m | 負の一軸性 |
| 方解石 | — | CaCO₃ | 3̄m | 負の一軸性 |
| サファイア | — | α-Al₂O₃ | 3̄m | 負の一軸性 |
| フッ化マグネシウム | — | MgF₂ | 4/mmm | 正の一軸性 |
| バナジン酸イットリウム | YVO₄ | YVO₄ | 4/m | 正の一軸性 |

### 光学的等方媒質

屈折率は 1 つで、メソッドは `(wl_um, T_degC)` を取ります。

| 材料 | 組成式 |
|---|---|
| 溶融石英 | SiO₂（非晶質） |
| フッ化カルシウム | CaF₂（立方晶、m3̄m） |
| フッ化リチウム | LiF（立方晶、m3̄m） |
| フッ化バリウム | BaF₂（立方晶、m3̄m） |
| イットリウム・アルミニウム・ガーネット (YAG) | Y₃Al₅O₁₂（立方晶、m3̄m） |
| N-BK7, SF10, SF11, SF57（SCHOTT） | ホウケイ酸クラウン・高密度フリントガラス |
| セレン化亜鉛、硫化亜鉛 | ZnSe, ZnS（立方晶、4̄3m。CVD グレード） |
| シリコン、ゲルマニウム | Si, Ge（立方晶、m3̄m） |
| ダイヤモンド | C（立方晶、m3̄m） |

同一材料の複数のパラメータ化は交換可能ではありません。各クラスはそれぞれの
出典に忠実であり、出典ごとに範囲と精度が異なります。比較の記録は
[検証ページ](https://ndispers.readthedocs.io/en/latest/validation/)に
あります。

二軸結晶は誘電主平面ごとに 1 クラスです。角度引数はその面内で変化する方
（xy 面では φ、yz・zx 面では θ）であり、インスタンスの `theta_rad` /
`phi_rad` 属性がどちらかを示します。Sellmeier 方程式に温度項を持たない媒質
も、シグネチャ統一のため温度引数を受け取り、無視します。その `dndT` は 0 を
返し、初回使用時に `TemperatureWarning` を発します。Sellmeier の有効範囲外の
波長では `ValidityWarning` を発し、外挿値をそのまま返します。

## 並列処理

medium オブジェクトは使用後も picklable で、`multiprocessing`、
`concurrent.futures`、`joblib` のワーカーに渡せ、`pickle`/`shelve` で保存
できます。評価関数はクラスごとにキャッシュされるため、同じ結晶のインスタンス
を多数作るコストは小さいです。

## ドキュメント

- [ブラウザアプリ](https://akihiko-shimura.github.io/ndispers/) — 屈折率エクスプローラーと位相整合計算機。クライアントサイドで動作します。
- [チュートリアルノートブック](examples/basic_usage.ipynb) — 基本の作業手順。GitHub 上で閲覧できます。
- [検証](https://ndispers.readthedocs.io/en/latest/validation/) — 文献との比較。数値、出典、注意点。
- [ndispers.readthedocs.io](https://ndispers.readthedocs.io/en/latest/) — 規約（単位・角度・符号）と[媒質カタログ](https://ndispers.readthedocs.io/en/latest/api/crystals/)。
- [llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt) / [llms-full.txt](https://ndispers.readthedocs.io/en/latest/llms-full.txt) — AI アシスタント向けの API リファレンス。

## 開発

```
git clone https://github.com/akihiko-shimura/ndispers.git
cd ndispers
uv sync
uv run pytest
```

チュートリアルノートブックを実行するには `notebook` グループ（matplotlib,
IPython, JupyterLab, marimo）を追加します:

```
uv sync --group notebook
uv run jupyter lab examples/basic_usage.ipynb
```

リリースはバージョンタグの push で PyPI に公開されます。
[docs/RELEASING.md](docs/RELEASING.md) を参照してください。結晶の追加手順は
[AGENTS.md](AGENTS.md) にあります。

## ライセンス

MIT — [LICENSE](LICENSE) を参照してください。
