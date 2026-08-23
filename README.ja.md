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

## インストール

Python 3.9 以降が必要です。

```
pip install ndispers
```

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

`dk_sfg` と `pmFactor_sfg` は任意の角度での位相不整合と sinc² 位相整合因子を返します。

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

### 光学的等方媒質

屈折率は 1 つで、角度・偏光の引数はありません。メソッドは `(wl_um, T_degC)` を取ります。

| 材料 | 組成 |
|---|---|
| 溶融石英 | SiO₂（非晶質） |
| フッ化カルシウム | CaF₂（立方晶, m3̄m） |
| イットリウム・アルミニウム・ガーネット (YAG) | Y₃Al₅O₁₂（立方晶, m3̄m） |
| N-BK7（SCHOTT） | ホウケイ酸クラウンガラス |

二軸性結晶は誘電主平面ごとに 1 クラスを提供します。角度引数はその面で変化する方（xy では φ、yz と zx では θ）で、インスタンスの `theta_rad` / `phi_rad` 属性がどちらかを示します。Sellmeier 方程式に温度項のない媒質（α-BBO、方解石、サファイア、水晶、5% MgO:LiNbO₃、BiBO、中赤外結晶、YAG、N-BK7）もシグネチャ統一のため温度引数を取りますが無視し、`dndT` は 0 を返します。

## 並列処理

媒質オブジェクトは使用後も含めて pickle 可能なので、`multiprocessing`、`concurrent.futures`、`joblib` のワーカーにそのまま渡せ、`pickle`/`shelve` で保存できます。lambdify された分散関数はクラスごとにキャッシュされるため、同じ結晶のインスタンスを多数作っても安価です。

## ドキュメント

- **[ブラウザアプリ](https://akihiko-shimura.github.io/ndispers/)** — 屈折率エクスプローラと位相整合計算機。インストール不要、クライアントサイドで動作。
- [チュートリアルノートブック](examples/basic_usage.ipynb) — 結晶の Sellmeier 方程式の確認から位相整合曲線までの基本の流れ。GitHub 上でそのまま閲覧可。
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
