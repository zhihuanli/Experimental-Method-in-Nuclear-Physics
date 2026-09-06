document.querySelectorAll('[data-code-tabs]').forEach((group) => {
  const buttons = group.querySelectorAll('[data-lang]');
  const panels = group.querySelectorAll('[data-code-panel]');

  buttons.forEach((button) => {
    button.addEventListener('click', () => {
      const language = button.dataset.lang;

      buttons.forEach((item) => {
        item.classList.toggle('active', item === button);
      });
      panels.forEach((panel) => {
        panel.hidden = panel.dataset.codePanel !== language;
      });
    });
  });
});
