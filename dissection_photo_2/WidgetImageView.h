#ifndef WIDGETIMAGEVIEW_H
#define WIDGETIMAGEVIEW_H

#include <QWidget>
#include <QImage>

class QPropertyAnimation;

#ifndef RECT_REGION
#define RECT_REGION QPair<QPoint,QPoint>
#endif

class WidgetImageView : public QWidget
{
  Q_OBJECT
  // Animated by m_pPreviewFadeAnim; a property so tests can set it directly.
  Q_PROPERTY(qreal previewFadeOpacity READ GetPreviewFadeOpacity WRITE SetPreviewFadeOpacity)
public:
  // True when DISSECTION_PHOTO_REDUCED_MOTION=1, or on the offscreen platform
  // unless it is "0" (same rule as AppMotion::Reduced(), which the other apps
  // don't build). Read once.
  static bool MotionReduced();

  explicit WidgetImageView(QWidget *parent = nullptr);

  // Arms a one-shot crossfade for the next LoadImage(): the current frame fades
  // out over the newly loaded one. Called only by
  // MainWindow::OnButtonTogglePreview(), since navigation also loads with
  // bPreview. Any other load cancels the fade and drops the retained frame.
  void BeginPreviewCrossfade();

  qreal GetPreviewFadeOpacity() const
  {
    return m_dPreviewFadeOpacity;
  }

  void SetPreviewFadeOpacity(qreal d)
  {
    m_dPreviewFadeOpacity = d;
    update();
  }

  bool LoadImage(const QString& filename, const QString& mask = "",
                 const QList<QPoint>& points = QList<QPoint>(),
                 const QList<RECT_REGION>& regions = QList<RECT_REGION>());

  void SetEditedPoints(const QList<QPoint>& pts);

  void paintEvent(QPaintEvent* e);
  void mousePressEvent(QMouseEvent* e);
  void mouseMoveEvent(QMouseEvent* e);
  void mouseReleaseEvent(QMouseEvent* e);
  void wheelEvent(QWheelEvent* e);
  void resizeEvent(QResizeEvent* e);
  void keyPressEvent(QKeyEvent* e);
  void keyReleaseEvent(QKeyEvent* e);

  bool eventFilter(QObject *watched, QEvent *event);

  QList<QPoint> GetEditedPoints()
  {
    return m_listPoints;
  }

  QList<RECT_REGION> GetEditedRegions()
  {
    return m_listRegions;
  }

  QString GetFilename()
  {
    return m_sFilename;
  }

  QImage GetImageInView()
  {
    return m_image;
  }

  // EM_SELECT_MASK is no longer used; kept so the values stay put.
  enum EditMode { EM_POINT = 0, EM_REGION, EM_CALIBRATION, EM_SELECT_MASK, EM_EDIT_MASK };

  // The tool selected in the toolbar. Changed only by SetMaskTool(), never by
  // modifiers.
  enum MaskTool { MT_BRUSH = 0, MT_ERASER, MT_BUCKET };

  // What a click would do now. Separate from MaskTool so flood overrides never
  // change the toolbar selection.
  //
  // Mask-editing bindings (EM_EDIT_MASK only). Left drag pans and right drag
  // zooms in every mode, so mask operations need a modifier:
  //
  //   Ctrl  + Left         apply the selected tool (Q / W / E pick it)
  //   Left, Shift + Left   pan (see mousePressEvent)
  //   Ctrl  + Right        flood fill
  //   Ctrl + Shift + Right flood erase
  enum MaskAction { MA_NONE = -1, MA_DRAW = 0, MA_ERASE, MA_FLOOD_FILL, MA_FLOOD_ERASE };

signals:
  void LastRegionEdited(int n);
  void PointEdited();
  // A point click on a corrected preview, where points are locked.
  void PointRejected();
  // The mask on screen went from saved to edited.
  void MaskEdited();

public slots:
  // Background behind the image, visible in the fit margin. MainWindow sets it to
  // the panel tone for the theme.
  void SetBackgroundColor(const QColor& c)
  {
    m_colorBackground = c;
    update();
  }

  // No points on a corrected preview: they would be in the wrong image's coordinates.
  void SetPointsLocked(bool b) { m_bPointsLocked = b; }
  // Drawn on since loaded or saved.
  bool IsMaskEdited() const { return m_bMaskEdited; }

  void SetNumberOfExpectedPoints(int n)
  {
    m_nNumberOfExpectedPoints = n;
    // Repaint so the connecting line appears when switching point count.
    update();
  }

  void SetEditMode(int n);

  // See MaskTool above. The only thing that ever changes the selection.
  void SetMaskTool(int n)
  {
    m_nMaskTool = n;
    m_nActiveAction = MaskActionFor(Qt::NoButton, Qt::NoModifier);
    update();
  }

  void ShowMessage(const QString& msg)
  {
    m_sMessage = msg;
    update();
  }

  void HideMessage()
  {
    m_sMessage.clear();
    update();
  }

  void SetOverlay(const QImage& overlay_image);
  void ClearEdits();

  void SetMaskOpacity(double val);

  void SetBrushSize(int n)
  {
    m_nBrushSize = n;
  }

  // Eraser size, separate from the brush size (see CurrentToolSize()).
  void SetEraserSize(int n)
  {
    m_nEraserSize = n;
  }

  // The size for the current action: brush size while drawing, eraser size while
  // erasing. Keyed to m_nActiveAction, not m_bErasing (only true mid-stroke), so
  // the hover preview is sized right.
  int CurrentToolSize() const
  {
    return (m_nActiveAction == MA_ERASE) ? m_nEraserSize : m_nBrushSize;
  }

  bool SaveMaskIfEdited();

  void UndoLastMaskEdit();

  void RedoLastMaskEdit();

  void Reset();

  void Clear();

  void UndoLastPoint();

  void RedoLastPoint();

  // Removes/restores the last CC region rectangle (EM_REGION). MainWindow undoes
  // the matching MaskProcessor fill.
  void UndoLastRegion();

  void RedoLastRegion();

  void SetCombinedMaskOpacity(double val);

private:
  void PrepareImage();
  void UpdateScaledImage(bool bSmooth = false);
  QPoint ScreenToImage(const QPoint& pt);
  void SetAlphaByMask(QImage& image);
  void FreeHandOnMaskImage(const QPoint& scr_pt1, const QPoint& scr_pt2);
  void UpdatePointOnMaskImage(const QPoint& pt_in);
  // Marks the mask edited; MaskEdited() on the first edit since load or save.
  void SetMaskEdited();
  void FloodFillMaskImage(const QPoint& scr_pt, bool bErase);
  void UpdateAll();

  // Stops any fade and releases the retained frame. Called by LoadImage() and
  // Clear().
  void CancelPreviewCrossfade();

  // Runs the fade armed by BeginPreviewCrossfade(), once the incoming image
  // is in place. No-op if nothing was retained.
  void StartPreviewCrossfade();

  // Resolves the bindings documented on MaskAction. Qt::NoButton asks what a left
  // click would do with the current modifiers (for the hover preview).
  int MaskActionFor(Qt::MouseButton button, Qt::KeyboardModifiers mods) const;

  // Applies one mask binding. False if button + modifiers are not a mask binding,
  // so the caller pans or zooms.
  bool BeginMaskAction(QMouseEvent* e);

  // Recomputes m_nActiveAction from the held modifiers. No-op outside
  // EM_EDIT_MASK.
  void RefreshMaskAction(Qt::KeyboardModifiers mods);

  // Draws the tool glyph beside the hover circle (nothing if the icon is missing).
  // Separate from paintEvent() so tests can call it at chosen positions.
  // dCircleRadius is the drawn circle's radius in logical px; the badge is offset
  // off it.
  void DrawToolCursorIcon(QPainter& p, const QPointF& pt, int nAction, double dCircleRadius);

  // Cached white glyph for a MaskAction, or a null pixmap. Loaded on first use,
  // not in paintEvent().
  const QPixmap& ToolCursorIcon(int nAction);

  QString   m_sFilename;
  QString   m_sMaskFilename;
  QImage    m_image;
  QImage    m_imageMask;
  QImage    m_imageScaled;
  QImage    m_imageOverlay;
  QImage    m_imageOriginal;
  double    m_dScale;
  double    m_dOldScale;
  QPoint    m_ptOffset;
  QPoint    m_ptOldOffset;
  QPoint    m_ptPress;
  QPoint    m_ptPrev;
  bool      m_bPanning;
  bool      m_bZooming;
  bool      m_bDrawing;
  bool      m_bErasing;

  int       m_nMaskTool = MT_BRUSH;      // toolbar selection (see MaskTool)
  int       m_nActiveAction = MA_DRAW;   // what a click would do now (see MaskAction)
  int       m_nNumberOfExpectedPoints;
  int       m_nEditMode;
  bool      m_bPointsLocked = false;
  QList<QPoint> m_listPoints;
  QList<QPoint> m_listRedoPoints;
  QList<RECT_REGION> m_listRegions;
  QList<RECT_REGION> m_listRedoRegions;
  QColor    m_colorPen;
  QColor    m_colorBackground = Qt::black;

  QString   m_sMessage;
  double    m_dMaskOpacity;
  double    m_dCombinedMaskOpacity;
  int      m_nBrushSize;
  int      m_nEraserSize = 10; // matches sliderEraserSize's initial value in MainWindow's constructor

  bool    m_bMaskEdited;
  bool    m_bModifierDown;
  QList<QImage>  m_listMaskUndoBuffer;
  QList<QImage>  m_listMaskRedoBuffer;

  // Crossfade state. m_imageFadeFrom is the retained outgoing frame (null when no
  // fade runs); m_dPreviewFadeOpacity is the incoming image's share (1.0 = done).
  QImage    m_imageFadeFrom;
  qreal     m_dPreviewFadeOpacity = 1.0;

  // Glyph cache indexed by MaskAction (0..3). m_bToolIconsLoaded limits a missing
  // resource to one failed load. Members, not function statics: a static QPixmap
  // would outlive QGuiApplication.
  QPixmap   m_pixToolCursorIcons[4];
  bool      m_bToolIconsLoaded = false;
  // Null pixmap returned by reference for an out-of-range action.
  QPixmap   m_pixNullToolIcon;
  bool      m_bPreviewFadePending = false;
  QPropertyAnimation* m_pPreviewFadeAnim = nullptr;
};

#endif // WIDGETIMAGEVIEW_H
