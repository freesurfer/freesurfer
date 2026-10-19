#ifndef IMAGEGRIDVIEW_H
#define IMAGEGRIDVIEW_H

#include <QWidget>
#include <QStringList>
#include <QPixmap>
#include <QColor>
#include <QVector>

class QStackedWidget;
class QLabel;
class QPushButton;
class QThumbnailWidget;
class QListWidgetItem;
class QListWidget;

// Start page preview in widgetMain: a thumbnail grid (loaded asynchronously by
// QThumbnailWidget/MakeThumbnailThread). Double-click opens a lightbox with a
// carousel strip and a Back button.
//
// MainWindow uses SetImages()/Clear(), the photo-picking API (SetSelectable(),
// SelectedFiles(), selectionChanged) and shows or hides the widget.
class ImageGridView : public QWidget
{
  Q_OBJECT
public:
  explicit ImageGridView(QWidget* parent = nullptr);

  // Loads the files into the grid and carousel and returns to the grid page.
  void SetImages(const QStringList& files);

  // Clears all loaded thumbnails/images and resets to an empty grid page,
  // e.g. after the user clears the input folder on the Start page.
  void Clear();

  // Applies the light/dark panel tones used by MainWindow::SetDarkMode().
  void SetDarkMode(bool bDark);

  // Photo picking. While selectable, each cell has a tick badge: click it (or
  // press Space) to leave a photo out or bring it back; left-out cells are
  // dimmed. A bar counts the selection and offers Select all / Select none.
  // SetImages() selects every photo.
  void SetSelectable(bool bSelectable);
  bool IsSelectable() const { return m_bSelectable; }
  bool IsSelected(int n) const;
  void SetSelected(int n, bool bSelected);
  void SetAllSelected(bool bSelected);
  // The files passed to SetImages(), in the same order, minus the ones left
  // out. All of them while not selectable.
  QStringList SelectedFiles() const;
  int SelectedCount() const;
  const QStringList& Files() const { return m_listFiles; }

signals:
  // Emitted whenever a photo is picked or left out (also by SetImages()).
  void selectionChanged(int nSelected, int nTotal);

protected:
  void resizeEvent(QResizeEvent* e) override;
  // Watches the grid viewport's resizes (see UpdateGridMetrics()).
  bool eventFilter(QObject* obj, QEvent* e) override;
  // Paints the panel background itself: a plain QWidget's stylesheet background is
  // not guaranteed to paint.
  void paintEvent(QPaintEvent* e) override;

protected slots:
  void OnGridItemDoubleClicked(QListWidgetItem* item);
  void OnCarouselIndexChanged(int n);
  void OnBackClicked();

private:
  void ShowLightbox(int n);
  void UpdateLightboxPixmap();
  // Sizes cells so leftover width is shared evenly between columns.
  void UpdateGridMetrics();
  // Where the tick badge sits inside a cell whose item rect is 'rcItem'.
  // bSmall: the carousel's smaller badge.
  static QRect BadgeRect(const QRect& rcItem, bool bSmall = false);
  // Refreshes everything that shows the selection: the grid bar, both
  // views' badges and the lightbox's Selected button and count.
  void UpdateSelectionBar();
  // Badge clicks on a thumbnail view (grid or carousel). True if handled.
  bool HandleBadgeMouse(QListWidget* pView, QEvent* e, bool bSmall);

  friend class GridSelectionDelegate;

  QStringList       m_listFiles;
  int               m_nCurrentIndex;
  bool              m_bUpdatingGridMetrics;
  QPixmap           m_fullPixmap;
  bool              m_bSelectable;
  QVector<bool>     m_vSelected;
  bool              m_bDark;

  QStackedWidget*   m_pStack;
  QThumbnailWidget* m_pGrid;
  QWidget*          m_pGridPage;
  QWidget*          m_pSelectionBar;
  QLabel*           m_pSelectionLabel;
  QPushButton*      m_pSelectAll;
  QPushButton*      m_pSelectNone;
  QWidget*          m_pLightboxPage;
  // Lightbox bands: top bar, image area and carousel band, each a widget so it
  // can paint its own background.
  QWidget*          m_pTopBar;
  QWidget*          m_pImageArea;
  QWidget*          m_pCarouselBand;
  QLabel*           m_pImageLabel;
  // Lightbox top bar, right side: whether the photo on screen is selected
  // (click to change) and how many are.
  QLabel*           m_pLightboxCount;
  QPushButton*      m_pLightboxSelect;
  QPushButton*      m_pBackButton;
  QThumbnailWidget* m_pCarousel;
  QColor            m_colorPanelBg;
};

#endif // IMAGEGRIDVIEW_H
