/* Render both activity views before portfolio.js initializes translation and search. */
(() => {
  'use strict';

  const script = document.currentScript;
  const siteRoot = new URL('../', script.src);
  const records = [...(window.activityRecords || [])]
    .sort((a, b) => b.date.localeCompare(a.date));

  function element(tag, className, text) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  }

  function translated(value) {
    const span = element('span', '', value.ja);
    span.dataset.ja = value.ja;
    span.dataset.en = value.en || value.ja;
    return span;
  }

  function category(record, className) {
    const label = window.activityCategories?.[record.category]
      || { ja: record.category, en: record.category };
    const node = translated(label);
    if (className) node.className = className;
    return node;
  }

  function content(record) {
    const fragment = document.createDocumentFragment();
    fragment.append(translated(record.text));
    for (const link of record.links || []) {
      const url = new URL(link.href, siteRoot);
      if (!['http:', 'https:', 'file:'].includes(url.protocol)) continue;
      const anchor = element('a');
      anchor.href = url.href;
      anchor.append(translated(link.label));
      fragment.append(document.createTextNode('　'), anchor);
    }
    return fragment;
  }

  function date(record) {
    const time = element('time', '', record.date.replaceAll('-', '/'));
    time.dateTime = record.date;
    return time;
  }

  for (const view of document.querySelectorAll('[data-activity-view]')) {
    if (view.dataset.activityView === 'recent') {
      const limit = Number(view.dataset.limit);
      const count = Number.isInteger(limit) && limit >= 0 ? limit : 8;
      for (const record of records.slice(0, count)) {
        const row = element('article', 'news-row');
        const meta = element('div', 'news-date');
        meta.append(date(record), category(record));
        const body = element('div');
        body.append(content(record));
        row.append(meta, body);
        view.append(row);
      }
      continue;
    }

    if (view.dataset.activityView !== 'archive') continue;
    const years = document.querySelector('[data-activity-years]');
    if (years) years.append(element('span', '', 'YEAR'));
    let currentYear;
    let list;
    for (const record of records) {
      const year = record.date.slice(0, 4);
      if (year !== currentYear) {
        currentYear = year;
        const section = element('section', 'report-year');
        section.id = `year-${year}`;
        section.setAttribute('aria-labelledby', `heading-${year}`);
        const heading = element('h2', '', year);
        heading.id = `heading-${year}`;
        heading.append(element('span', '', 'ACTIVITY REPORTS'));
        list = element('dl', 'report-list');
        section.append(heading, list);
        view.append(section);
        if (years) {
          const link = element('a', '', year);
          link.href = `#year-${year}`;
          years.append(link);
        }
      }
      const row = element('div', 'report-entry');
      const meta = element('dt');
      meta.append(date(record), category(record, 'report-category'));
      const body = element('dd');
      body.append(content(record));
      row.append(meta, body);
      list.append(row);
    }
    if (years) years.append(element('span', 'report-count', `${records.length} REPORTS`));
  }
})();
