# 活動報告の更新

記事は `activity-data.js` の `window.activityRecords` に追加・修正します。
トップと活動報告ページの HTML を個別に編集する必要はありません。

```js
{
  "date": "2026-11-10",
  "category": "Research",
  "text": {
    "ja": "活動報告の本文。",
    "en": "English description."
  },
  "links": [
    {
      "href": "https://example.com/report",
      "label": { "ja": "詳細", "en": "Details" }
    }
  ]
}
```

- 日付は `YYYY-MM-DD`。表示時に日付の新しい順に並びます。
- `en` が空文字の場合は日本語を表示します。HTML タグは使わず、本文は文字列で記述します。
- リンクがなければ `"links": []` にします。
- サイト内リンクはルートの `index.html` からの相対パスで指定します（例: `Data/F300/F340/index.html`）。
- トップの件数は `index.html` の `data-limit="8"` で変更できます。検索対象は表示した N 件です。
- 活動報告ページは全件を表示し、年別リンクと総件数も自動更新します。

`activity-data.js` → `activities.js` → `portfolio.js` の順で読み込みます。
データを通常の script として読み込むため、公開サイトでも `file://` で直接開く場合でも使えます。
JavaScript を無効にすると記事一覧は表示されません。
