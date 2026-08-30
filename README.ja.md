<img src="https://raw.githubusercontent.com/akihiko-shimura/ndispers/main/docs/assets/logo.svg" alt="" width="132" align="right">

# ndispers

[English](https://github.com/akihiko-shimura/ndispers/blob/main/README.md) | **日本語**

[![PyPI](https://img.shields.io/pypi/v/ndispers)](https://pypi.org/project/ndispers/)
[![test](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml/badge.svg)](https://github.com/akihiko-shimura/ndispers/actions/workflows/test.yml)
[![apps](https://github.com/akihiko-shimura/ndispers/actions/workflows/pages.yml/badge.svg)](https://akihiko-shimura.github.io/ndispers/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/akihiko-shimura/ndispers/blob/main/examples/basic_usage.ipynb)

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
実効非線形係数 d<sub>eff</sub> を計算します。

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
  曲線。全非中心対称結晶で Miller 波長スケーリングつきの d<sub>eff</sub>。
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

[チュートリアルノートブック](https://colab.research.google.com/github/akihiko-shimura/ndispers/blob/main/examples/basic_usage.ipynb)が、
Sellmeier 方程式の確認から位相整合曲線までの手順をプロットつきで示して
います。Colab で開くのでインストール不要、ブラウザからそのまま実行できます。
要点:

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

スカラー入力は float を返し、numpy 配列入力はブロードキャストして配列を
返します。他の分散メソッド（`GD`, `GV`, `ng`, `GVD`, `TOD`, `GVM`,
ウォークオフ, `dndT`）と位相整合メソッド（`pmAngles_sfg`, `acceptance_sfg`,
`qpm_period_sfg`, `tuning_dfg`, `deff_sfg` など）も同じ形式の引数を取ります。
`help(bbo)` が材料データと参考文献を表示します。

分散の掃引、位相整合角、OPO チューニング、d<sub>eff</sub> などの実例は
[ドキュメント](https://ndispers.readthedocs.io/en/latest/)と
[チュートリアルノートブック](https://colab.research.google.com/github/akihiko-shimura/ndispers/blob/main/examples/basic_usage.ipynb)にあります。

## 収録媒質

非線形結晶: β-BBO, LBO, KTP, BiBO, CLBO, KDP, DKDP, KBBF, RBBF, LB4,
LiIO₃, ZGP, AgGaS₂, AgGaSe₂, 水晶, MgO 添加 LiNbO₃（コングルエント・
ストイキオメトリック）, MgO:LiTaO₃。線形複屈折結晶: α-BBO, 方解石,
サファイア, MgF₂, YVO₄。等方媒質: 溶融石英, CaF₂, LiF, BaF₂, YAG, N-BK7,
SF10/11/57, ZnSe, ZnS, Si, Ge, ダイヤモンド。

各媒質は `nd.media.crystals` または `nd.media.glasses` のクラスです。複数の
パラメータ化がある結晶は出典ごとに別クラスになっています。
[媒質カタログ](https://ndispers.readthedocs.io/en/latest/api/crystals/)が
全クラスの Sellmeier 方程式・有効範囲・参考文献を一覧し、`nd.catalog()` が
同じ情報をデータとして返します。

## ドキュメント

- [ブラウザアプリ](https://akihiko-shimura.github.io/ndispers/) — 屈折率エクスプローラーと位相整合計算機。クライアントサイドで動作します。
- [チュートリアルノートブック](https://colab.research.google.com/github/akihiko-shimura/ndispers/blob/main/examples/basic_usage.ipynb) — 基本の作業手順。Colab で開きます（[ソース](examples/basic_usage.ipynb)）。
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
