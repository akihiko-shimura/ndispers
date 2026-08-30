<img src="https://raw.githubusercontent.com/akihiko-shimura/ndispers/main/docs/assets/logo.svg" alt="" width="132" align="right">

# ndispers

[English](https://github.com/akihiko-shimura/ndispers/blob/main/README.md) | **日本語**

[![PyPI](https://img.shields.io/pypi/v/ndispers)](https://pypi.org/project/ndispers/)
[![test](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml/badge.svg)](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml)
[![apps](https://github.com/akihiko-shimura/ndispers/actions/workflows/pages.yml/badge.svg)](https://akihiko-shimura.github.io/ndispers/)

*ndispers* は、非線形光学・超高速光学で使われる結晶やガラスの屈折率分散を計算する Python パッケージです。文献で報告された Sellmeier 方程式と熱光学係数 (dn/dT) に基づいています。

計算できる量:

- 屈折率
- 群遅延
- 群速度
- 群屈折率
- 群速度分散 (GVD)
- 3 次分散 (TOD)
- ウォークオフ角
- dn/dT
- d²n/dT²

変数:

1. 光の波長
2. 異方性結晶の誘電主軸に対する波数ベクトルの極角 (theta) または方位角 (phi)
3. 媒質の温度
4. 光の偏光（常光線・異常光線）

結晶には非線形光学のメソッドもあります:

- 位相不整合 Δk
- 位相整合角
- 位相整合因子 sinc²(Δk·L/2)
- 実効非線形係数 d_eff（テンソル成分は Miller 則で使用波長にスケーリング。パッケージ内の非中心対称結晶すべてで利用可）

## なぜ ndispers か

この種の計算のツールは既にいくつかあります。最も有名で網羅的なのは
Arlee V. Smith による *SNLO*
（[as-photonics.com](https://as-photonics.com/products/snlo/)、Windows 用 GUI）です。
ほかにも屈折率の Web アプリ（[refractiveindex.info](https://refractiveindex.info/)）、
位相整合計算の Web アプリ（[toolbox.lightcon.com](http://toolbox.lightcon.com/)）、
iOS アプリ（[iPhasematch](https://apps.apple.com/jp/app/iphasematch/id492370060)）が
あります。デスクトップ・Web・モバイルと「電卓」は揃っていますが、いずれも
**あなたのプログラムの部品にはなりません**。

*ndispers* は、SNLO 型計算の中核 — 屈折率・分散・位相整合・d_eff — を
**Python ライブラリとして**提供します。すべてのメソッドは numpy 配列を受けて
配列を返す普通の関数なので、波長・角度・温度を一度に掃引して、自作の
数値シミュレーション（パルス伝搬、OPO 設計、熱解析）や Jupyter ノートブックに
そのまま組み込めます。そして全係数が出典つきで公開されています。`bbo?` と
打てば「どの論文のどの表の値か」まで確認でき、論文の Methods にクラス名を
1 行書けば計算が再現可能になります。

普通の、内省可能な Python ライブラリであることには、もう 1 つの帰結があります。
**ndispers はエージェント対応です**。AI アシスタントが必要とするものはすべて
自己記述的になっています — 媒質は `dir()` で列挙でき、各 docstring は出典と
有効範囲を持ち、[llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt)
が API 全体・単位・落とし穴を 1 ページに凝縮しています。このリンクを
アシスタントに渡して、日本語で頼むだけ —「BBO の 1064 nm タイプ I SHG の
位相整合角を求めて、温度依存性をプロットして」— あなたは 1 行もコードを
書かずに、インストールから計算・プロットまで代行させられます。GUI や Web の
電卓をこの方法で駆動することはできません。

対象は非線形結晶だけではありません。SNLO が扱わない線形の複屈折結晶
（α-BBO、方解石、水晶、YVO₄ など）や光学ガラス（溶融石英、CaF₂、SF10 など）も
同じインターフェースで揃えており、波長板・偏光子・分散管理まで含めた
超高速光学の計算を 1 つのパッケージで完結させることを目指しています。

**特長の一覧**

- **組み込み可能** — GUI ではなくライブラリ。全メソッドが numpy 配列を受け、
  事前コンパイル済みの numpy 関数として動くため、実行時依存は numpy のみ
  （sympy 不要、初回のコード生成待ちなし）。medium オブジェクトは picklable で、
  `multiprocessing` / `joblib` にそのまま渡せます。
- **エージェント対応** — GUI が隠すものすべてが内省可能な Python として
  露出しており、[llms.txt](https://ndispers.readthedocs.io/en/latest/llms.txt)
  が AI アシスタントに API 全体を 1 ページで渡します。リンクを渡して自然言語で
  頼むだけ。電卓が対話型になります。
- **出典の透明性** — 各媒質の docstring に Sellmeier の出典と波長・温度の
  有効範囲を明記。`constants` プロパティで係数値そのものも見られます。
  値の出所が追えない、ということがありません。
- **同一結晶の複数 Sellmeier** — 例えば LBO は Kato 1994 / Kato–Kuroda 2018 /
  Ghosh 1995 / メーカー値（Castech, Newlight）を別クラスとして併載。
  文献間の食い違いを自分で比較・評価できます。
- **SNLO にない材料** — 線形複屈折結晶（α-BBO、方解石、水晶、サファイア、
  MgF₂、YVO₄）と等方媒質（溶融石英、CaF₂、BaF₂、YAG、N-BK7、SF ガラス、
  Si、Ge など）を同一 API で提供。
- **温度微分は解析的** — dn/dT・d²n/dT² を数値差分ではなく式の微分から計算。
  温度チューニングや熱レンズ解析に。
- **非線形光学の道具立て** — 位相不整合 Δk、位相整合角の直接解、
  sinc² 位相整合因子、Miller 則で使用波長にスケールした d_eff
  （全非中心対称結晶）。
- **拡張が容易** — 点群ごとの基底クラスを継承し、Sellmeier と d 係数を
  1 ファイル書けば新しい結晶を追加できます。追加リクエストや貢献は
  [GitHub](https://github.com/akihiko-shimura/ndispers/issues) へ。
- **クロスプラットフォーム** — `pip install ndispers` だけで
  Linux / macOS / Windows、クラスタでも Colab でも動きます。
- **転写の検証つき** — 係数は文献からの独立再抽出で照合し、回帰テストで固定。
  写し間違いは運ではなく仕組みで防いでいます。

## インストール

Python 3.10 以降が必要です。

```
pip install ndispers
```

依存パッケージは numpy だけです。分散式は sympy の式として書かれていますが、実際に評価される関数は事前に生成して同梱しています。式そのもの（`n_expr` など）を扱う場合や、媒質を自分でサブクラス化して評価する場合は `pip install ndispers[sym]` で sympy を追加してください。

## クイックスタート

結晶の Sellmeier 方程式の確認から位相整合曲線まで、図付きの案内は [チュートリアルノートブック](examples/basic_usage.ipynb)（GitHub 上でそのまま読めます）を参照してください。要点:

β-BBO 結晶のオブジェクトを作る:

```python
>>> import ndispers as nd
>>> bbo = nd.media.crystals.BetaBBO_Eimerl1987()
```

屈折率を計算する:

```python
>>> bbo.n(0.532, 0, 25, pol='o')
1.674884049110459
>>> bbo.n(0.532, 3.1416/2, 25, pol='e')
1.5554658787539917
```

4 つの引数は順に

1. 波長 (µm)
2. theta 角 (rad)
3. 温度 (°C)
4. 偏光（`pol='o'` または `'e'`、常光線/異常光線。既定は `'o'`）

theta = 0（光軸に沿った伝搬）では常光線と異常光線の屈折率は一致します。スカラー入力には float を返し、各引数は numpy 配列も受け付けて同じ形の配列を返します:

```python
>>> import numpy as np
>>> bbo.n(np.arange(0.2, 1.5, 0.2), 0, 25, pol='o')
array([1.89001202, 1.69328828, 1.66985875, 1.6612891 , 1.65633946,
       1.65252664, 1.64903624])
```

他の分散メソッド（`GD`, `GV`, `ng`, `GVD`, `TOD`, `woa_theta`, `woa_phi`, `dndT`, `dndT2`）も同じ引数を取ります:

```python
>>> bbo.GVD(0.8, 0, 25, pol='o')   # fs^2/mm
71.86403019943364
```

光学的に等方な媒質（ガラスと立方晶の CaF₂）は角度も偏光も無いので、波長と温度だけを取ります:

```python
>>> silica = nd.media.glasses.FusedSilica()
>>> silica.n(1.064, 20)
1.4495857898590634
```

材料情報、Sellmeier 方程式、係数の出典文献を見るには `help(bbo)`、IPython/Jupyter なら `bbo?` とします。`bbo.constants` は係数の値そのものを返します。

## 位相整合

和周波発生の位相整合角を直接解きます。β-BBO、25 °C、1064 nm の Type-I SHG:

```python
>>> bbo.pmAngles_sfg(1.064, 1.064, 25, deg=True)
{'wl3': 0.532,
 'ooe': {'theta': [22.884169498625802], 'phi': None},
 'eeo': {'theta': [], 'phi': None},
 'oee': {'theta': [32.56045545648089], 'phi': None},
 'eoe': {'theta': [32.56045545648089], 'phi': None},
 'eoo': {'theta': [], 'phi': None},
 'oeo': {'theta': [], 'phi': None}}
```

差周波発生・光パラメトリック増幅/発振では同じ相互作用をポンプ側から入ります: `pmAngles_dfg(wl_p, wl_s, T)` が位相整合角（と遊休波長）、`tuning_dfg(wl_p, angle, T, pol_s, pol_i, pol_p)` が固定角度で位相整合する信号/遊休の組（OPO チューニング曲線の 1 点）を返し、`dk_dfg`・`pmFactor_dfg`・`qpm_period_dfg`・`deff_dfg` は SFG 版と対になります（波の順序は signal, idler, pump）。1064 nm 励起の x カット KTP なら `KTP_zx().tuning_dfg(1.064, np.pi/2, 25, 'o', 'e', 'o')` がおなじみの非臨界の組 1571 / 3298 nm を返します。`qpm_period=` を渡せば同じ関数で周期分極反転結晶の温度チューニングが得られます。

`dk_sfg` と `pmFactor_sfg` は任意の角度での位相不整合と sinc² 位相整合因子を、`qpm_period_sfg` は擬似位相整合周期を返します。5% MgO:LiNbO₃ の x 軸伝搬・全波異常光線（d33）での 1064 nm SHG なら `MgOLN_Zelmon1997().qpm_period_sfg(1.064, 1.064, np.pi/2, 20, 'e', 'e', 'e')` でおなじみの約 7.0 µm です。

その角度、φ = 90° カットでの実効非線形係数:

```python
>>> bbo.deff_sfg(1.064, 1.064, np.radians(22.88), np.radians(90), 25, 'o', 'o', 'e')
-1.9937...  # pm/V（全体の符号は規約）; d22 = 2.2 pm/V at 1.064 µm SHG
            # (Shoji et al. 1999)、ウォークオフ込み
>>> bbo.d_sfg("d22", 0.8, 0.8, 25)
2.3300...   # Miller 則による 800 nm SHG の d22
```

非中心対称の結晶はすべてこれを持ちます（KDP, CLBO, LBO, KTP, ...）。符号と角度の規約は `help(bbo.deff_sfg)` と *Conventions* ページに、各結晶の係数の出典と Miller スケーリングの検証状況は `_d_note` にあります。

## 収録媒質

すべての媒質は `nd.media.crystals` または `nd.media.glasses` のクラスです。クラス名と、各媒質の Sellmeier 方程式・有効範囲・文献は [媒質カタログ](https://ndispers.readthedocs.io/en/latest/api/crystals/) にあります。いくつかの結晶は複数のパラメータ化を持ち、係数の出典（文献またはベンダー）の名前が付いています。信頼したい出典を選んでください。

### 非線形光学結晶

非中心対称。位相整合、許容幅、d_eff が利用できます。

| 材料 | 略称 | 組成 | 点群 | 光学分類 |
|---|---|---|---|---|
| β-ホウ酸バリウム | β-BBO | β-BaB₂O₄ | 3m | 負の一軸性 |
| 三ホウ酸リチウム | LBO | LiB₃O₅ | mm2 | 二軸性、主平面 3 面 |
| チタン酸リン酸カリウム | KTP | KTiOPO₄ | mm2 | 二軸性、主平面 3 面 |
| 三ホウ酸ビスマス | BiBO | BiB₃O₆ | 2 | 二軸性（単斜晶）、主平面 3 面 |
| ホウ酸セシウムリチウム | CLBO | CsLiB₆O₁₀ | 4̄2m | 負の一軸性 |
| リン酸二水素カリウム | KDP | KH₂PO₄ | 4̄2m | 負の一軸性 |
| 重水素化リン酸二水素カリウム | DKDP, KD*P | KD₂PO₄ | 4̄2m | 負の一軸性 |
| フッ化ホウ酸ベリリウムカリウム | KBBF | KBe₂BO₃F₂ | 32 | 負の一軸性 |
| フッ化ホウ酸ベリリウムルビジウム | RBBF | RbBe₂BO₃F₂ | 32 | 負の一軸性 |
| 四ホウ酸リチウム | LB4（LBT とも） | Li₂B₄O₇ | 4mm | 負の一軸性 |
| ヨウ素酸リチウム | — | LiIO₃ | 6 | 負の一軸性 |
| リン化ゲルマニウム亜鉛 | ZGP | ZnGeP₂ | 4̄2m | 正の一軸性、中赤外 |
| 硫化ガリウム銀 | AGS | AgGaS₂ | 4̄2m | 負の一軸性、中赤外 |
| セレン化ガリウム銀 | AGSe | AgGaSe₂ | 4̄2m | 負の一軸性、中赤外 |
| α-水晶 | — | SiO₂ | 32 | 正の一軸性 |
| ニオブ酸リチウム、5% MgO 添加コングルエント | MgO:LN | MgO:LiNbO₃ | 3m | 負の一軸性、両光線 |
| ニオブ酸リチウム、1% MgO 添加定比 | MgO:SLN | MgO:LiNbO₃ | 3m | 負の一軸性、異常光線のみ |
| タンタル酸リチウム、1% MgO 添加定比 | MgO:SLT | MgO:LiTaO₃ | 3m | 負の一軸性 |

### 複屈折光学結晶

中心対称なので 2 次非線形性はありません。窓・偏光子・補償板向けの分散、ウォークオフ、熱光学。

| 材料 | 略称 | 組成 | 点群 | 光学分類 |
|---|---|---|---|---|
| α-ホウ酸バリウム | α-BBO | α-BaB₂O₄ | 3̄m | 負の一軸性 |
| 方解石 | — | CaCO₃ | 3̄m | 負の一軸性 |
| サファイア | — | α-Al₂O₃ | 3̄m | 負の一軸性 |
| フッ化マグネシウム | — | MgF₂ | 4/mmm | 正の一軸性 |
| バナジン酸イットリウム | YVO₄ | YVO₄ | 4/m | 正の一軸性 |

### 光学的等方媒質

屈折率は 1 つで、角度・偏光の引数はありません。メソッドは `(wl_um, T_degC)` を取ります。

| 材料 | 組成 |
|---|---|
| 溶融石英 | SiO₂（非晶質） |
| フッ化カルシウム | CaF₂（立方晶, m3̄m） |
| フッ化リチウム | LiF（立方晶, m3̄m） |
| フッ化バリウム | BaF₂（立方晶, m3̄m） |
| イットリウム・アルミニウム・ガーネット (YAG) | Y₃Al₅O₁₂（立方晶, m3̄m） |
| N-BK7, SF10, SF11, SF57（SCHOTT） | ホウケイ酸クラウン・高分散フリントガラス |
| セレン化亜鉛、硫化亜鉛 | ZnSe, ZnS（立方晶, 4̄3m; CVD 材） |
| シリコン、ゲルマニウム | Si, Ge（立方晶, m3̄m） |
| ダイヤモンド | C（立方晶, m3̄m） |

複数のパラメータ化がある材料では、それらは**等価ではありません**。各クラスはそれぞれの出典に忠実で、出典の質には差があります。どれを使うべきかは [検証記録](https://ndispers.readthedocs.io/en/latest/validation/) を参照してください。

二軸性結晶は誘電主平面ごとに 1 クラスを提供します。角度引数はその面で変化する方（xy では φ、yz と zx では θ）で、インスタンスの `theta_rad` / `phi_rad` 属性がどちらかを示します。Sellmeier 方程式に温度項のない媒質（α-BBO、方解石、サファイア、水晶、MgF₂、5% MgO:LiNbO₃、BiBO、中赤外結晶、YAG、N-BK7、LiF、BaF₂ ほか等方材料一式）もシグネチャ統一のため温度引数を取りますが無視し、`dndT` は 0 を返します。

## 並列処理

媒質オブジェクトは使用後も含めて pickle 可能なので、`multiprocessing`、`concurrent.futures`、`joblib` のワーカーにそのまま渡せ、`pickle`/`shelve` で保存できます。lambdify された分散関数はクラスごとにキャッシュされるため、同じ結晶のインスタンスを多数作っても安価です。

## ドキュメント

- **[ブラウザアプリ](https://akihiko-shimura.github.io/ndispers/)** — 屈折率エクスプローラと位相整合計算機。インストール不要、クライアントサイドで動作。
- [チュートリアルノートブック](examples/basic_usage.ipynb) — 結晶の Sellmeier 方程式の確認から位相整合曲線までの基本の流れ。GitHub 上でそのまま閲覧可。
- [検証記録 (Validation)](https://ndispers.readthedocs.io/en/latest/validation/) — 文献値との照合結果を数値と出典つきで、および監査で判明した注意点。
- [ndispers.readthedocs.io](https://ndispers.readthedocs.io/en/latest/) — 規約（単位、角度、符号の規約）と、全結晶・ガラスの式・有効範囲・文献を載せた [媒質カタログ](https://ndispers.readthedocs.io/en/latest/api/crystals/)。

## 開発

```
git clone https://github.com/akihiko-shimura/ndispers.git
cd ndispers
uv sync
uv run pytest
```

チュートリアルノートブックを動かすには `notebook` グループ（matplotlib, IPython, JupyterLab, marimo）を追加します:

```
uv sync --group notebook
uv run jupyter lab examples/basic_usage.ipynb
```

リリースはバージョンタグの push で PyPI に公開されます。[docs/RELEASING.md](docs/RELEASING.md) を参照。

## ライセンス

MIT — [LICENSE](LICENSE) を参照。
