# 広報の手順書（要: 本人アカウント）

2026-08-30 時点。準備済みの成果物: PyPI メタデータ・Colab バッジ（反映済み）、
JOSS 草稿 `paper/`、conda レシピ `docs/dev/conda-forge/meta.yaml`。
以下は上から順に。1 と 2 だけ先にやれば残りは急がない。

## 1. Zenodo 連携（5 分）

1. <https://zenodo.org> に GitHub でログイン
2. GitHub タブ → リポジトリ一覧で `ndispers` を ON
3. 連携済み（2026-08-30）。次のリリースタグ push で自動的に DOI が発行される
   （既存の v0.19.0 に付けたい場合は GitHub で Release を「編集→保存」し直すと
   Zenodo が拾うことがある。だめなら次リリースからで良い）
4. 発行された Concept DOI（全バージョン共通の方）を `CITATION.cff` に
   `doi:` として追記し、README のバッジに追加してもよい

## 1.5 RTD リダイレクト（5 分）

llms.txt 規約に従うエージェントはドメイン直下
`https://ndispers.readthedocs.io/llms.txt` を最初に試す。RTD ダッシュボード →
プロジェクト → Redirects で `/llms.txt` → `/en/latest/llms.txt` の
Exact Redirect を 1 本追加（`/llms-full.txt` も同様に）。

## 2. JOSS 投稿（30 分 + 査読数週間）

1. `paper/paper.md` の ORCID（今は 0000-0000-0000-0000）と所属を実物に差し替え
   → commit & push
2. <https://joss.theoj.org/papers/new> に GitHub でログインし、
   リポジトリ URL と `paper/` のブランチ (main) を指定して submit
3. 査読は GitHub の issue 上で進む。機能追加要求ではなく
   ドキュメント・テスト・インストール確認が中心
4. 採択されたら JOSS の DOI を `CITATION.cff` の `preferred-citation` に追加

## 3. conda-forge（30 分 + レビュー数日）

1. <https://github.com/conda-forge/staged-recipes> を自分のアカウントに fork
2. `recipes/ndispers/meta.yaml` として `docs/dev/conda-forge/meta.yaml` を
   コピー（sha256 は 0.19.0 用に記入済み。別バージョンを出すなら
   `curl -sL <sdist URL> | shasum -a 256` で更新）
3. PR を出す。テンプレの checklist に沿って自分を maintainer として明記
4. マージ後は feedstock リポジトリが作られ、以後の PyPI リリースに自動追従
5. **マージ後にやること**: README（英日）の Installation に
   `conda install -c conda-forge ndispers` を追記（Claude に頼めば即時）

## 4. 掲載依頼（各 10 分、効果は検索流入）

- **RP Photonics**: <https://www.rp-photonics.com> の Encyclopedia
  “nonlinear crystal materials” 等の記事末尾にあるソフトウェア掲載の
  問い合わせ先へ、1 段落の紹介文と URL をメール
- **awesome リスト**: GitHub で `awesome-photonics` / `awesome-optics` を検索し、
  該当セクションに 1 行追加する PR（説明は README の 1 文目を流用）

## 5. 実例コンテンツ（継続、1 本 30 分）

- Photonics/Physics StackExchange で phase matching / Sellmeier の質問に
  ndispers での 10 行解答を書く（宣伝ではなく解答が主、リンクは末尾に）
- Qiita/Zenn に日本語記事 1 本。切り口の案:
  「SNLO を Python から使いたかった人へ」「AI エージェントに位相整合計算を
  やらせる（llms.txt を渡すだけ）」
