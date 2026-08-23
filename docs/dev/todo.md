# やり残し・次にやること（2026-08-23 時点）

> 各項目の背景は括弧内の計画書に。ここは「何が残っているか」の一覧だけ。

## 機能

- **NOPA（非共線 OPA）** — `dfg_plan.md` §5。理論ノート `ncpm_theory.tex` を先に書く
  （ベクトル Δk、群速度整合条件 v_s = v_i cos Ω、Type II の分枝、d_eff の波ごとの k̂）。
  検証値は BBO 400 nm 励起 θ ≈ 31.5°, α ≈ 3.7° を Cerullo & De Silvestri 2003 /
  Wilhelm, Piel, Riedle 1997 の PDF（手元にあり）で内部角・外部角を確定してから。
- **QPM の温度チューニング曲線** λ_s(T) at fixed Λ — アプリは Λ(λ_s) しか描いていない。
  PPLN ユーザーが一番見る図。`tuning_dfg` に `qpm_period=` 引数（Δk − 2π/Λ = 0 を解く）を
  足せば一関数で済む。
- **二軸結晶の主平面外**（θ, φ とも自由）— `Medium` は面ごとのクラスで主平面に限定。
  NOPA の面外版と合わせて要る。
- 位相整合アプリの **Δk/sinc² 曲線データの CSV 書き出し**（レポート検討時に見送り）。
- `pmAngles_*` は 0–90° のみ。単斜晶 BiBO の 90–180° 側は d_eff セルが別扱いしているが、
  `pmAngles` の返り値には含まれない。

## 検証（docs/validation.md に行を足す）

- OPO チューニング曲線の文献照合: BBO 355 nm Type I（Cheng/Bosenberg/Tang 1988 または Dmitriev）、
  KTP xz Type II 1064 nm、LBO xy 355 nm と NCPM 温度チューニング、PPLN（Myers et al. 1995, Λ = 28–31 µm）。
  いずれも PDF が手元に無い。現状は KTP x カット NCPM の 1 点（1570.8/3297.6 nm）と自己無撞着テストのみ。
- `deff_theory.tex` に「三波の置換対称性（Kleinman + Miller）」の節を 1 つ追記（`dfg_plan.md` 段階 0、未実施）。

## 材料（`materials_plan.md` の保留リスト）

KTA（Fenimore 1995）, RTP（Mikami 2009）, CBO, SBO, KNbO₃（軸規約の衝突あり）, GaSe/LGS,
Yb:KGW/KYW, OP-GaAs/GaP, **LBGO**（Oxide 社ラインナップで唯一追加価値あり — データシート待ち）。

## 保守

- `AGS_Takaoka1999` は deprecated、1.0 で削除。
- 数値が変わるリリースでは `docs/validation.md` の表とフッタの版数を同じコミットで更新する
  （フッタの版数 ≠ `__version__` はテストが止める。表の中身の再確認は人手）。
- 媒質・`*_expr` を触ったら `uv run python tools/compile_media.py`（テストが止める）。
- ~~apps の `marimo check` 警告~~ → `__generated_with` を付けて解消。
- ~~Pages の CDN キャッシュ~~ → `apps/README.md` に確認手順を記載。
