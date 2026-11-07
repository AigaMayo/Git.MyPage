/* Portfolio interactions. Loaded with defer after the document is parsed. */
(() => {
  'use strict';

  const languageStorageKey = 'kai-ogura-x2-language';
  const menu = document.querySelector('.menu');
  const navigation = document.querySelector('#navigation');
  const navigationLinks = [...navigation.querySelectorAll('a')];
  const languageButtons = [...document.querySelectorAll('[data-language]')];
  const activityRows = [...document.querySelectorAll('.news-row')];
  const search = document.querySelector('#search');
  const resultCount = document.querySelector('#result-count');
  const emptyState = document.querySelector('#empty');
  const dialog = document.querySelector('dialog');
  const largeImage = dialog.querySelector('img');
  const imageCaption = dialog.querySelector('p');
  const photos = [...document.querySelectorAll('.photo')];
  const progress = document.querySelector('.progress');
  const sections = [
    document.querySelector('.hero'),
    ...document.querySelectorAll('section.content'),
  ];
  const navigationLabels = {
    en: ['HOME', 'PROFILE', 'ACTIVITY', 'GALLERY ↗'],
    ja: ['ホーム', 'プロフィール', '活動履歴', 'ギャラリー ↗'],
  };
  const descriptions = {
    ja: '小倉 可維のポートフォリオ。MEMSジャイロスコープ、宇宙開発、ものづくりの活動記録。',
    en: 'Kai Ogura portfolio: MEMS gyroscope research, space engineering, and other activities.',
  };
  let scrollUpdatePending = false;

  // Mobile navigation
  function setMenuOpen(open) {
    const japanese = document.documentElement.lang === 'ja';
    navigation.classList.toggle('open', open);
    menu.setAttribute('aria-expanded', String(open));
    menu.textContent = japanese
      ? (open ? '閉じる −' : 'メニュー ＋')
      : (open ? 'Close −' : 'Menu ＋');
  }

  menu.addEventListener('click', () => {
    setMenuOpen(!navigation.classList.contains('open'));
  });
  navigationLinks.forEach(link => {
    link.addEventListener('click', () => setMenuOpen(false));
  });
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape') setMenuOpen(false);
  });

  // Reveal content when it enters the viewport.
  if ('IntersectionObserver' in window) {
    document.documentElement.classList.add('js');
    const observer = new IntersectionObserver(entries => {
      entries.forEach(entry => {
        if (!entry.isIntersecting) return;
        entry.target.classList.add('visible');
        observer.unobserve(entry.target);
      });
    }, { threshold: 0.06 });
    document.querySelectorAll('.reveal').forEach(element => observer.observe(element));
  }

  // Scroll progress and current navigation item
  function updateScrollState() {
    const scrollableHeight = document.documentElement.scrollHeight - innerHeight;
    progress.style.width = `${scrollableHeight > 0 ? scrollY / scrollableHeight * 100 : 0}%`;
    let activeIndex = 0;
    sections.forEach((section, index) => {
      if (section.getBoundingClientRect().top < innerHeight * 0.4) activeIndex = index;
    });
    navigationLinks.forEach((link, index) => {
      const active = index === activeIndex;
      link.classList.toggle('active', active);
      if (active) link.setAttribute('aria-current', 'location');
      else link.removeAttribute('aria-current');
    });
    scrollUpdatePending = false;
  }

  window.addEventListener('scroll', () => {
    if (scrollUpdatePending) return;
    scrollUpdatePending = true;
    requestAnimationFrame(updateScrollState);
  }, { passive: true });
  window.addEventListener('resize', updateScrollState);

  // Search the activity records in the currently displayed language.
  function filterActivities() {
    const query = search.value.trim().toLocaleLowerCase();
    let visibleCount = 0;
    activityRows.forEach(row => {
      row.hidden = !row.textContent.toLocaleLowerCase().includes(query);
      if (!row.hidden) visibleCount++;
    });
    const suffix = document.documentElement.lang === 'ja' ? ' 件' : ' RECORDS';
    resultCount.textContent = `${visibleCount} / ${activityRows.length}${suffix}`;
    emptyState.hidden = visibleCount !== 0;
  }
  search.addEventListener('input', filterActivities);

  // Gallery dialog
  function updateImageCaption(photo) {
    largeImage.alt = photo.dataset.caption;
    imageCaption.textContent = photo.dataset.caption;
  }

  photos.forEach(photo => {
    photo.addEventListener('click', () => {
      largeImage.src = photo.dataset.full;
      updateImageCaption(photo);
      dialog.showModal();
    });
  });
  dialog.querySelector('.close').addEventListener('click', () => dialog.close());
  dialog.addEventListener('click', event => {
    if (event.target !== dialog) return;
    const bounds = dialog.getBoundingClientRect();
    const outside = event.clientX < bounds.left || event.clientX > bounds.right
      || event.clientY < bounds.top || event.clientY > bounds.bottom;
    if (outside) dialog.close();
  });
  dialog.addEventListener('close', () => largeImage.removeAttribute('src'));

  // Translation data lives alongside the corresponding HTML content.
  function setLanguage(language) {
    const lang = language === 'ja' ? 'ja' : 'en';
    document.documentElement.lang = lang;
    document.querySelectorAll('[data-en][data-ja]').forEach(element => {
      element.textContent = element.dataset[lang];
    });
    for (const attribute of ['alt', 'aria-label', 'placeholder', 'data-caption']) {
      const translationAttribute = `data-i18n-${attribute}-${lang}`;
      document.querySelectorAll(`[${translationAttribute}]`).forEach(element => {
        element.setAttribute(attribute, element.getAttribute(translationAttribute));
      });
    }
    navigationLinks.forEach((link, index) => {
      link.textContent = navigationLabels[lang][index];
    });
    languageButtons.forEach(button => {
      button.setAttribute('aria-pressed', String(button.dataset.language === lang));
    });
    document.title = lang === 'ja' ? 'Kai Ogura — ポートフォリオ' : 'Kai Ogura — Portfolio';
    document.querySelector('meta[name="description"]').content = descriptions[lang];
    setMenuOpen(false);
    filterActivities();
    if (dialog.open) {
      const photo = photos.find(item => (
        new URL(item.dataset.full, document.baseURI).href === largeImage.src
      ));
      if (photo) updateImageCaption(photo);
    }
    try {
      localStorage.setItem(languageStorageKey, lang);
    } catch {
      // Language switching also works when browser storage is unavailable.
    }
  }

  languageButtons.forEach(button => {
    button.addEventListener('click', () => setLanguage(button.dataset.language));
  });
  let savedLanguage = 'ja';
  try {
    savedLanguage = localStorage.getItem(languageStorageKey) || 'ja';
  } catch {
    // Use Japanese when browser storage is unavailable.
  }
  // A fixed year scale; event descriptions open only on request.
  function initializeTimeline() {
    const track = document.querySelector('#experience-track');
    if (!track) return;
    const history = track.closest('.history');
    const controls = history.querySelector('.timeline-controls');
    const previous = controls.querySelector('.timeline-prev');
    const next = controls.querySelector('.timeline-next');
    const position = controls.querySelector('.timeline-position');
    const events = [...track.querySelectorAll('.timeline-event')];
    const panel = history.querySelector('.timeline-detail-panel');
    const descriptions = [...panel.querySelectorAll('.experience-detail')];
    const close = panel.querySelector('.timeline-detail-close');
    const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)');
    let selected = null;
    let pending = false;

    function update() {
      const maximum = Math.max(0, track.scrollWidth - track.clientWidth);
      previous.disabled = track.scrollLeft <= 2;
      next.disabled = track.scrollLeft >= maximum - 2;
      const first = Math.max(2017, Math.min(2026, 2017 + Math.floor((track.scrollLeft - 50) / 220)));
      const last = Math.max(first, Math.min(2026, 2017 + Math.floor((track.scrollLeft + track.clientWidth - 50) / 220)));
      position.textContent = `${first} – ${last}`;
      pending = false;
    }

    function move(direction) {
      track.scrollBy({ left: direction * 440, behavior: reducedMotion.matches ? 'instant' : 'smooth' });
    }
    function closeDetails(restoreFocus = false) {
      events.forEach(event => event.setAttribute('aria-expanded', 'false'));
      descriptions.forEach(description => { description.hidden = true; });
      panel.hidden = true;
      if (restoreFocus && selected) selected.focus({ preventScroll: true });
      selected = null;
    }
    events.forEach(event => {
      event.addEventListener('click', () => {
        const wasOpen = event.getAttribute('aria-expanded') === 'true';
        closeDetails();
        if (wasOpen) return;
        selected = event;
        event.setAttribute('aria-expanded', 'true');
        document.getElementById(event.getAttribute('aria-controls')).hidden = false;
        panel.hidden = false;
        panel.scrollIntoView({ block: 'nearest', behavior: reducedMotion.matches ? 'instant' : 'smooth' });
      });
    });
    close.addEventListener('click', () => closeDetails(true));
    history.addEventListener('keydown', event => {
      if (event.key === 'Escape' && !panel.hidden) closeDetails(true);
    });
    previous.addEventListener('click', () => move(-1));
    next.addEventListener('click', () => move(1));
    track.addEventListener('keydown', event => {
      if (event.target !== track) return;
      if (event.key === 'ArrowLeft' || event.key === 'ArrowRight') {
        event.preventDefault();
        move(event.key === 'ArrowLeft' ? -1 : 1);
      }
    });
    track.addEventListener('scroll', () => {
      if (pending) return;
      pending = true;
      requestAnimationFrame(update);
    }, { passive: true });
    window.addEventListener('resize', update);
    controls.hidden = false;
    requestAnimationFrame(() => {
      track.scrollLeft = track.scrollWidth - track.clientWidth;
      update();
    });
  }
  initializeTimeline();
  setLanguage(savedLanguage);
  updateScrollState();
})();
