#include "ReconPage.h"
#include "CommonDef.h"
#include "SegmentedControl.h"
#include "ProjectManager.h"

#include <QApplication>
#include <QBoxLayout>
#include <QButtonGroup>
#include <QCheckBox>
#include <QClipboard>
#include <QCryptographicHash>
#include <QDesktopServices>
#include <QDir>
#include <QDoubleValidator>
#include <QEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontMetrics>
#include <QFrame>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QLocale>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProcessEnvironment>
#include <QProgressBar>
#include <QPushButton>
#include <QRadioButton>
#include <QRegularExpression>
#include <QResizeEvent>
#include <QScrollArea>
#include <QSettings>
#include <QShowEvent>
#include <QStyle>
#include <QTextBlock>
#include <QTextCursor>
#include <QTextLayout>
#include <QTimer>
#include <QUrl>

namespace
{
// The output folder the project gets by default. The other two are what the
// earlier steps in this pipeline actually write.
const char* kOutputSubdir = "recon";
const char* kPhotoSubdir = "deformed";
const char* kMaskSubdir = "connected_components";

// The resolution (mm per pixel) the correction step writes its photos at:
// always passed, never asked for.
const char* kDefaultPhotoResolution = "0.1";

const int kLabelWidth = 128;      // the side label column when the column is wide
const int kWideMin = 380;         // below this form width the labels go above their controls
const int kControlHeight = 28;
const int kMaxLogBlocks = 5000;
const int kLogMinHeight = 120;
const int kPreviewMinLines = 7;

// Reconstruction's step colour (Correction is red, Segmentation olive). Fills
// only; dark-theme ink uses the lighter kPurpleInkDark for contrast on #2a2a2a.
const char* kPurple = "#89448e";
const char* kPurpleHover = "#783b7c";
const char* kPurplePressed = "#68336b";
const char* kPurpleInkDark = "#c58fca";

// The command's "expert" loss weights: one float flag each, emitted only
// when filled in (utils.py get_arguments, cp_spacing_2d .. k_regularizer_sz).
const char* kExpertKeys[] = {"cp_spacing_2d", "cp_spacing_3d", "k_lncc_mri", "k_dice_mri", "k_dif_slice_loss",
                             "k_mesh_loss", "k_regularizer", "k_regularizer_nonlin", "k_regularizer_nonlin3d",
                             "k_regularizer_sz", nullptr};
const char* kExpertLabels[] = {"Control point spacing (2D)", "Control point spacing (3D)", "LNCC weight (MRI)",
                               "Dice weight (MRI)", "Slice difference weight", "Mesh loss weight",
                               "Regulariser weight", "Nonlinear regulariser", "Nonlinear regulariser (3D)",
                               "Size regulariser", nullptr};

// Fields holding paths. The FreeSurfer launcher forwards an unquoted $@, so a
// path with a space reaches the parser as two arguments.
const char* kPathKeys[] = {"ref_mesh", "ref_mri",
                           "ref_mri_synthseg", "ref_mri_synthsr", "deform_recon_dir", "input_roi_dir", "weights",
                           nullptr};

QString ShellQuote(const QString& s)
{
  static const QRegularExpression safe("^[A-Za-z0-9_./=:,+-]+$");
  if (safe.match(s).hasMatch())
    return s;
  QString q = s;
  q.replace("'", "'\\''");
  return "'" + q + "'";
}

bool ParseNumber(const QString& text, double* out)
{
  bool ok = false;
  *out = QLocale::c().toDouble(text.trimmed(), &ok);
  return ok;
}

bool ParseInt(const QString& text, int* out)
{
  bool ok = false;
  *out = QLocale::c().toInt(text.trimmed(), &ok);
  return ok;
}

QString FreesurferHome()
{
  return QProcessEnvironment::systemEnvironment().value("FREESURFER_HOME");
}

bool Overridden()
{
  return !QProcessEnvironment::systemEnvironment().value("DISSECTION_PHOTO_RECON_CMD").isEmpty();
}

QFont MonoFont(const QFont& base)
{
  QFont f = base;
  FsSetFontFamilies(f, QStringList() << "Menlo" << "DejaVu Sans Mono" << "Noto Sans Mono" << "Liberation Mono");
  f.setStyleHint(QFont::Monospace);
  f.setPixelSize(12);
  return f;
}

// A property the page's sheet selects on, re-polished so the change shows.
void SetKind(QWidget* w, const char* name, const QString& value)
{
  if (w->property(name).toString() == value)
    return;
  w->setProperty(name, value);
  w->style()->unpolish(w);
  w->style()->polish(w);
  w->update();
}

}   // namespace

// ==== Local widgets ====

// One form row: a label and a cell (the control over its caption). Side by
// side while the column is wide, the label above the control when it is not.
class ReconFormRow : public QWidget
{
public:
  ReconFormRow(const QString& label, QWidget* control, QLabel* caption, bool bIndent, QWidget* parent = nullptr)
    : QWidget(parent), m_pLabel(nullptr), m_bIndent(bIndent)
  {
    m_pBox = new QBoxLayout(QBoxLayout::LeftToRight, this);
    m_pBox->setContentsMargins(0, 0, 0, 0);
    m_pBox->setSpacing(8);
    if (!label.isEmpty() || bIndent)
    {
      m_pLabel = new QLabel(label, this);
      m_pLabel->setObjectName("reconLabel");
      m_pLabel->setWordWrap(true);
      m_pLabel->setMinimumHeight(kControlHeight);
      m_pLabel->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      m_pBox->addWidget(m_pLabel, 0, Qt::AlignTop);
    }
    QWidget* cell = new QWidget(this);
    QVBoxLayout* v = new QVBoxLayout(cell);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(4);
    v->addWidget(control);
    if (caption)
      v->addWidget(caption);
    m_pBox->addWidget(cell, 1);
    SetNarrow(false);
  }

  void SetNarrow(bool bNarrow)
  {
    if (!m_pLabel)
      return;
    m_pBox->setDirection(bNarrow ? QBoxLayout::TopToBottom : QBoxLayout::LeftToRight);
    m_pBox->setSpacing(bNarrow ? 2 : 8);
    if (bNarrow)
    {
      m_pLabel->setMinimumWidth(0);
      m_pLabel->setMaximumWidth(QWIDGETSIZE_MAX);
      m_pLabel->setMinimumHeight(0);
      // An indent-only label (the inventory line) has nothing to say on a
      // line of its own.
      m_pLabel->setVisible(!m_bIndent);
    }
    else
    {
      m_pLabel->setFixedWidth(kLabelWidth);
      m_pLabel->setMinimumHeight(kControlHeight);
      m_pLabel->setVisible(true);
    }
  }

private:
  QBoxLayout* m_pBox;
  QLabel* m_pLabel;
  bool m_bIndent;
};

// ==== Construction ====

ReconPage::ReconPage(QWidget* parent)
  : QWidget(parent)
{
  const bool bDark = palette().color(QPalette::Window).lightness() < 128;

  m_pValidate = new QTimer(this);
  m_pValidate->setSingleShot(true);
  m_pValidate->setInterval(150);
  connect(m_pValidate, &QTimer::timeout, this, &ReconPage::OnValidate);
  m_pSave = new QTimer(this);
  m_pSave->setSingleShot(true);
  m_pSave->setInterval(500);
  connect(m_pSave, &QTimer::timeout, this, &ReconPage::SaveSettings);
  m_pTick = new QTimer(this);
  m_pTick->setInterval(1000);
  connect(m_pTick, &QTimer::timeout, this, &ReconPage::OnTick);
  m_pCopied = new QTimer(this);
  m_pCopied->setSingleShot(true);
  m_pCopied->setInterval(1200);
  QVBoxLayout* v = new QVBoxLayout(this);
  v->setContentsMargins(0, 0, 0, 0);
  v->setSpacing(0);

  m_pScroll = new QScrollArea(this);
  m_pScroll->setObjectName("reconScroll");
  m_pScroll->setWidgetResizable(true);
  m_pScroll->setFrameShape(QFrame::NoFrame);
  m_pScroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  m_pScroll->setWidget(BuildForm());
  m_pScroll->viewport()->installEventFilter(this);
  v->addWidget(m_pScroll, 1);
  v->addWidget(BuildFooter());

  m_pPane = BuildPane();
  connect(m_pCopied, &QTimer::timeout, this, [this]() { m_pCopy->setText("Copy command"); });

  SetDarkMode(bDark);
  ResetForm();
  UpdateVisibility();
  UpdatePreview();
  SetState(Idle);
  CheckAvailability();
  OnValidate();
}

ReconPage::~ReconPage()
{
  if (m_pSave->isActive())
    SaveSettings();
  // A running reconstruction outlives the page (Switch Project); ReconRunner
  // stops it at quit.
  if (m_pRunner)
    disconnect(m_pRunner, nullptr, this, nullptr);
  // The pane was handed to MainWindow's layout; if it never was (a harness
  // building the page on its own), it is still ours to free.
  if (m_pPane && !m_pPane->parent())
    delete m_pPane;
}

void ReconPage::AddSection(QVBoxLayout* into, const QString& strText, bool bFirst)
{
  QLabel* l = new QLabel(strText);
  l->setObjectName("reconSection");
  if (!bFirst)
    into->addSpacing(12);
  into->addWidget(l);
}

void ReconPage::AddRow(QVBoxLayout* into, const QString& strKey, const QString& strLabel,
                       QWidget* pControl, const QString& strHelp)
{
  QLabel* cap = nullptr;
  if (!strKey.isEmpty())
  {
    cap = new QLabel(strHelp);
    cap->setObjectName("reconCaption");
    cap->setWordWrap(true);
    cap->setMinimumHeight(16);
    // No help text: the caption only appears for an error (SetCaptionText()).
    cap->setVisible(!strHelp.isEmpty());
    m_captions[strKey] = cap;
    m_help[strKey] = strHelp;
  }
  const bool bIndent = strKey.isEmpty();
  ReconFormRow* row = new ReconFormRow(strLabel, pControl, cap, bIndent);
  into->addWidget(row);
  if (!strKey.isEmpty())
    m_rows[strKey] = row;
  m_listRows << row;
}

QLineEdit* ReconPage::MakeNumberEdit(const QString& strKey, const QString& strPlaceholder, bool bInteger)
{
  QLineEdit* e = new QLineEdit;
  e->setObjectName("reconMonoField");
  e->setFont(MonoFont(font()));
  e->setPlaceholderText(strPlaceholder);
  e->setFixedHeight(kControlHeight);
  if (bInteger)
    e->setValidator(new QIntValidator(e));
  else
  {
    QDoubleValidator* dv = new QDoubleValidator(e);
    dv->setLocale(QLocale::c());
    dv->setNotation(QDoubleValidator::StandardNotation);
    e->setValidator(dv);
  }
  m_edits[strKey] = e;
  connect(e, &QLineEdit::textChanged, this, &ReconPage::OnFormEdited);
  return e;
}

// The loss weights are often written like 1e-3; the parser's float() takes
// that as typed, and so does ParseNumber(), so the text is passed on as is.
QLineEdit* ReconPage::MakeExpertEdit(const QString& strKey)
{
  QLineEdit* e = MakeNumberEdit(strKey, QString(), false);   // UpdateVisibility() fills in the default
  QDoubleValidator* dv = new QDoubleValidator(e);
  dv->setLocale(QLocale::c());
  dv->setNotation(QDoubleValidator::ScientificNotation);
  e->setValidator(dv);
  return e;
}

QWidget* ReconPage::MakeUnitRow(QLineEdit* pEdit, const QString& strUnit, int nWidth)
{
  QWidget* w = new QWidget;
  QHBoxLayout* h = new QHBoxLayout(w);
  h->setContentsMargins(0, 0, 0, 0);
  h->setSpacing(8);
  pEdit->setMinimumWidth(0);
  pEdit->setMaximumWidth(nWidth);
  h->addWidget(pEdit, 1);
  if (!strUnit.isEmpty())
  {
    QLabel* u = new QLabel(strUnit);
    u->setObjectName("reconUnit");
    h->addWidget(u);
  }
  h->addStretch(0);
  return w;
}

QWidget* ReconPage::MakePathRow(const QString& strKey, const QString& strPlaceholder, bool bDirectory)
{
  QWidget* w = new QWidget;
  QHBoxLayout* h = new QHBoxLayout(w);
  h->setContentsMargins(0, 0, 0, 0);
  h->setSpacing(8);
  QLineEdit* e = new QLineEdit;
  e->setObjectName("reconField");
  e->setPlaceholderText(strPlaceholder);
  e->setFixedHeight(kControlHeight);
  e->setMinimumWidth(60);
  m_edits[strKey] = e;
  connect(e, &QLineEdit::textChanged, this, &ReconPage::OnFormEdited);
  connect(e, &QLineEdit::textChanged, e, [e](const QString& t) { e->setToolTip(t); });
  h->addWidget(e, 1);
  QPushButton* b = new QPushButton(QString::fromUtf8("Choose\xe2\x80\xa6"));
  b->setProperty("reconKind", "ghost");
  b->setFixedHeight(kControlHeight);
  b->setCursor(Qt::PointingHandCursor);
  h->addWidget(b);
  connect(b, &QPushButton::clicked, this, [this, e, bDirectory, strKey]()
  {
    const QString start = e->text().isEmpty() ? m_strRoot : e->text();
    QString path;
    if (bDirectory)
      path = QFileDialog::getExistingDirectory(window(), "Choose a folder", start);
    else
    {
      QString filter;
      if (strKey == "ref_mesh") filter = "Surface meshes (*.stl);;All files (*)";
      else if (strKey.startsWith("ref_mri")) filter = "Volumes (*.nii *.nii.gz *.mgz *.mgh);;All files (*)";
      else if (strKey == "weights") filter = "CSV files (*.csv *.txt);;All files (*)";
      path = QFileDialog::getOpenFileName(window(), "Choose a file", start, filter);
    }
    if (!path.isEmpty())
      e->setText(path);
  });
  return w;
}

QCheckBox* ReconPage::MakeCheck(const QString& strKey, const QString& strText)
{
  QCheckBox* c = new QCheckBox(strText);
  c->setObjectName("reconCheck");
  c->setCursor(Qt::PointingHandCursor);
  m_checks[strKey] = c;
  connect(c, &QCheckBox::toggled, this, &ReconPage::OnFormEdited);
  return c;
}

SegmentedControl* ReconPage::MakePill(const QStringList& items, int nDefault)
{
  SegmentedControl* s = new SegmentedControl;
  // Items and the default first, the connect after: SetItems() resets to
  // index 0, so setting the default would otherwise count as an edit.
  s->SetItems(items);
  s->SetCurrentIndex(nDefault, false);
  s->SetAccent(QColor(kPurple));
  // Chosen against the left column (#f0f0f0 light, #2a2a2a dark) the way the
  // defaults were chosen against the top bar: a visibly recessed track.
  s->SetTrackColors(QColor("#dadada"), QColor("#1b1b1b"));
  s->setFixedHeight(kControlHeight + 2);
  s->setCursor(Qt::PointingHandCursor);
  connect(s, &SegmentedControl::currentIndexChanged, this, &ReconPage::OnFormEdited);
  return s;
}

QWidget* ReconPage::BuildForm()
{
  QWidget* content = new QWidget;
  content->setObjectName("reconContent");
  QVBoxLayout* sl = new QVBoxLayout(content);
  // Right margin leaves room for the scroll bar when the form overflows.
  sl->setContentsMargins(0, 4, 12, 16);
  sl->setSpacing(8);

  QLabel* subtitle = new QLabel("Compute a 3D photo reconstruction of the preprocessed photographs. "
                                "Choose a reference to register to, set parameters, and run the reconstruction.");
  subtitle->setObjectName("reconSubtitle");
  subtitle->setWordWrap(true);
  sl->addWidget(subtitle);

  // Back sits under the description, as on the other steps' panels.
  m_pBack = new QPushButton(QString::fromUtf8("\xe2\x86\x90 Back"));
  m_pBack->setObjectName("reconBack");
  m_pBack->setProperty("reconKind", "ghost");
  m_pBack->setCursor(Qt::PointingHandCursor);
  m_pBack->setFixedHeight(kControlHeight);
  connect(m_pBack, &QPushButton::clicked, this, &ReconPage::BackRequested);
  sl->addWidget(m_pBack, 0, Qt::AlignLeft);
  sl->addSpacing(4);

  // The form proper: disabled as a whole while a run is in progress.
  m_pForm = new QWidget;
  QVBoxLayout* fl = new QVBoxLayout(m_pForm);
  fl->setContentsMargins(0, 0, 0, 0);
  fl->setSpacing(8);
  sl->addWidget(m_pForm);

  // Always this project's corrected photos and saved slabs; only the output
  // folder's name can change.
  AddSection(fl, "INPUTS");
  m_pInventory = new QLabel(QString::fromUtf8("\xe2\x80\x94"));
  m_pInventory->setObjectName("reconInventory");
  m_pInventory->setFont(MonoFont(font()));
  AddRow(fl, "inputs", "This project", m_pInventory,
         "The corrected photos and the slabs saved in Connected components.");
  {
    QLineEdit* e = new QLineEdit;
    e->setObjectName("reconMonoField");
    e->setFont(MonoFont(font()));
    e->setPlaceholderText(kOutputSubdir);
    e->setFixedHeight(kControlHeight);
    m_edits["output_name"] = e;
    connect(e, &QLineEdit::textChanged, this, &ReconPage::OnFormEdited);
    AddRow(fl, "output_name", "Output folder", MakeUnitRow(e, QString(), 192),
           "A folder inside this project, created if needed.");
  }

  AddSection(fl, "SLABS");
  m_pHemisphere = MakePill(QStringList() << "Left" << "Right" << "Both", 2);
  AddRow(fl, "hemisphere", "Hemisphere", m_pHemisphere, "Which side of the brain the slabs come from.");
  AddRow(fl, "slice_thickness", "Slab thickness",
         MakeUnitRow(MakeNumberEdit("slice_thickness", "for example 10", false), "mm", 96),
         "How thick each slab was cut. Be as exact as you can; it is refined during the run.");
  m_pFace = MakePill(QStringList() << "Anterior" << "Posterior", 0);
  AddRow(fl, "photos_of_posterior_side", "Photographed face", m_pFace,
         "Which cut face is in the photos: the front (anterior) or the back (posterior).");
  m_pOrder = MakePill(QStringList() << "Anterior first" << "Posterior first", 0);
  AddRow(fl, "order_posterior_to_anterior", "Photo order", m_pOrder,
         "The order the photos were numbered in: front to back, or back to front.");
  m_pTissue = MakePill(QStringList() << "Fixed" << "Fresh", 0);
  AddRow(fl, "fresh_tissue", "Tissue", m_pTissue, "Fresh tissue is allowed to deform more.");

  AddSection(fl, "REFERENCE");
  QLabel* refCap = new QLabel("A reference shape guides the stacking. Without one, the MNI atlas is used.");
  refCap->setObjectName("reconCaption");
  refCap->setWordWrap(true);
  fl->addWidget(refCap);
  {
    QWidget* radios = new QWidget;
    QGridLayout* rg = new QGridLayout(radios);
    rg->setContentsMargins(0, 0, 0, 0);
    rg->setHorizontalSpacing(16);
    rg->setVerticalSpacing(6);
    m_pRefNone = new QRadioButton(QString::fromUtf8("None \xe2\x80\x94 MNI atlas"));
    m_pRefMesh = new QRadioButton("Surface mesh");
    m_pRefMri = new QRadioButton("MRI scan");
    m_pRefBoth = new QRadioButton("Mesh and MRI");
    QButtonGroup* group = new QButtonGroup(this);
    int id = 0;
    foreach (QRadioButton* r, QList<QRadioButton*>() << m_pRefNone << m_pRefMesh << m_pRefMri << m_pRefBoth)
    {
      r->setObjectName("reconRadio");
      r->setCursor(Qt::PointingHandCursor);
      group->addButton(r, id++);
    }
    m_pRefNone->setChecked(true);
    rg->addWidget(m_pRefNone, 0, 0);
    rg->addWidget(m_pRefMesh, 0, 1);
    rg->addWidget(m_pRefMri, 1, 0);
    rg->addWidget(m_pRefBoth, 1, 1);
    rg->setColumnStretch(2, 1);
    fl->addWidget(radios);
    connect(group, static_cast<void (QButtonGroup::*)(QAbstractButton*)>(&QButtonGroup::buttonClicked),
            this, [this](QAbstractButton*) { OnFormEdited(); });
  }

  // A mode's fields keep their text while another mode is chosen; they are
  // simply not shown, not validated and not emitted.
  m_pMeshBlock = new QWidget;
  {
    QVBoxLayout* ml = new QVBoxLayout(m_pMeshBlock);
    ml->setContentsMargins(0, 4, 0, 0);
    ml->setSpacing(8);
    AddRow(ml, "ref_mesh", "Surface mesh", MakePathRow("ref_mesh", QString::fromUtf8("Choose a .stl file\xe2\x80\xa6"), false),
           "A 3D scan of the whole hemisphere taken before cutting.");
    QLineEdit* idx = MakeNumberEdit("mesh_reorient_with_indices", "999,1010,333333", true);
    idx->setValidator(nullptr);
    AddRow(ml, "mesh_reorient_with_indices", "Reorientation vertices", MakeUnitRow(idx, QString(), 192),
           QString::fromUtf8("Three vertex indices, separated by commas \xe2\x80\x94 see the PhotoTools wiki."));
  }
  fl->addWidget(m_pMeshBlock);

  m_pMriBlock = new QWidget;
  {
    QVBoxLayout* ml = new QVBoxLayout(m_pMriBlock);
    ml->setContentsMargins(0, 4, 0, 0);
    ml->setSpacing(8);
    AddRow(ml, "ref_mri", "MRI scan", MakePathRow("ref_mri", QString::fromUtf8("Choose a .nii or .mgz file\xe2\x80\xa6"), false),
           "An MRI of this brain. A SynthSeg segmentation is made for it if one is not given.");
    AddRow(ml, "ref_mri_synthseg", "SynthSeg", MakePathRow("ref_mri_synthseg", QString::fromUtf8("Choose a segmentation\xe2\x80\xa6"), false),
           "Written into the recon folder unless you choose an existing one. Computed if the file does not exist yet.");
    AddRow(ml, "ref_mri_synthsr", "SynthSR", MakePathRow("ref_mri_synthsr", QString::fromUtf8("Choose a SynthSR volume\xe2\x80\xa6"), false),
           "Written into the recon folder unless you choose an existing one. Computed if the file does not exist yet.");
    AddRow(ml, "low_field_synthsr", "", MakeCheck("low_field_synthsr", "Use low-field version of SynthSR"),
           "For low-field scans (for example Hyperfine).");
  }
  fl->addWidget(m_pMriBlock);

  AddSection(fl, "RUN SETTINGS");
  AddRow(fl, "threads", "Threads", MakeUnitRow(MakeNumberEdit("threads", "1", true), QString(), 192),
         "More threads can speed up processing. Set to -1 to use all cores available");
#ifndef Q_OS_MACOS
  // NVIDIA (CUDA) only, which Macs do not have.
  AddRow(fl, "gpu", "Graphics card",
         MakeUnitRow(MakeNumberEdit("gpu", QString::fromUtf8("none \xe2\x80\x94 use the processor"), true), QString(), 192),
         "The index of an NVIDIA card (0 for the first). Much faster if one is available.");
#endif

  fl->addSpacing(8);
  m_pAdvancedToggle = new QPushButton;
  m_pAdvancedToggle->setProperty("reconKind", "disclosure");
  m_pAdvancedToggle->setCursor(Qt::PointingHandCursor);
  fl->addWidget(m_pAdvancedToggle, 0, Qt::AlignLeft);
  connect(m_pAdvancedToggle, &QPushButton::clicked, this, [this]()
  {
    SetAdvancedOpen(!m_pAdvanced->isVisibleTo(m_pForm));
    if (!m_bLoading)
      m_pSave->start();
  });

  m_pAdvanced = new QWidget;
  {
    QVBoxLayout* al = new QVBoxLayout(m_pAdvanced);
    al->setContentsMargins(0, 0, 0, 0);
    al->setSpacing(8);
    AddSection(al, "IMAGES", true);
    AddRow(al, "equalize_images", "", MakeCheck("equalize_images", "Equalise photo contrast"),
           "Helps with low-contrast photos, but can reduce the quality of the machine learning imputation.");

    AddSection(al, "GEOMETRY");
    AddRow(al, "no_z_stretch", "", MakeCheck("no_z_stretch", "Lock slice thickness"),
           "Does not stretch the stack front to back, even when a mesh or MRI reference is used.");
    AddRow(al, "initial_stretch_factor_lr_photos", "Left-right stretch of photos",
           MakeUnitRow(MakeNumberEdit("initial_stretch_factor_lr_photos", "command default (1.0)", false), QString(), 192),
           "Useful for fresh tissue specimens. Default 1.0 does nothing.");
    AddRow(al, "stretch_factor_lr_mesh", "Left-right stretch of mesh",
           MakeUnitRow(MakeNumberEdit("stretch_factor_lr_mesh", "command default (1.0)", false), QString(), 192),
           "As above, applied to the mesh reference.");

    AddSection(al, "SLAB WEIGHTS");
    AddRow(al, "weights", "Slab weights file", MakePathRow("weights", QString::fromUtf8("Choose a .csv file\xe2\x80\xa6"), false),
           "One weight per slab (including anterior/posterior-most slab without a face), in photo order, separated by commas.");
    AddRow(al, "thickness_cap", "Maximum thickness", MakeUnitRow(MakeNumberEdit("thickness_cap", "no cap", false), "mm", 96),
           "Only used with a weights file. Caps the estimated thickness of any slab.");

    AddSection(al, "EXTRA VOLUMES");
    AddRow(al, "deform_recon_dir", "FreeSurfer recon-all folder of the template MRI", MakePathRow("deform_recon_dir", "recon-all subject folder", true),
           "Deforms that subject's surfaces and aparc+aseg, and the MRI's SynthSR, into photo space.");
    AddRow(al, "input_roi_dir", "Regions to deform", MakePathRow("input_roi_dir", "Folder of ROI masks", true),
           "Extra .nii or .mgz volumes to deform into the reconstruction.");

    AddSection(al, "OPTIMISATION");
    AddRow(al, "skip_bfgs", "", MakeCheck("skip_bfgs", "Skip BFGS"), "Use to skip BFGS finetuning.");

    al->addSpacing(4);
    m_pExpertToggle = new QPushButton;
    m_pExpertToggle->setProperty("reconKind", "disclosure");
    m_pExpertToggle->setCursor(Qt::PointingHandCursor);
    al->addWidget(m_pExpertToggle, 0, Qt::AlignLeft);
    m_pExpert = new QWidget;
    QVBoxLayout* xl = new QVBoxLayout(m_pExpert);
    xl->setContentsMargins(0, 0, 0, 0);
    xl->setSpacing(8);
    for (int i = 0; kExpertKeys[i]; ++i)
      AddRow(xl, kExpertKeys[i], kExpertLabels[i], MakeUnitRow(MakeExpertEdit(kExpertKeys[i]), QString(), 192),
             QString());
    al->addWidget(m_pExpert);
    connect(m_pExpertToggle, &QPushButton::clicked, this, [this]()
    {
      SetExpertOpen(!m_pExpert->isVisibleTo(m_pAdvanced));
      if (!m_bLoading)
        m_pSave->start();
    });
  }
  fl->addWidget(m_pAdvanced);
  SetAdvancedOpen(false);
  SetExpertOpen(false);

  sl->addStretch(1);
  return content;
}

// Pinned under the scrolling form, the way the other steps pin their primary
// action: the readiness caption, then the one primary button.
QWidget* ReconPage::BuildFooter()
{
  QWidget* footer = new QWidget;
  QVBoxLayout* v = new QVBoxLayout(footer);
  v->setContentsMargins(0, 0, 0, 16);
  v->setSpacing(8);
  QFrame* rule = new QFrame;
  rule->setObjectName("reconRule");
  rule->setFixedHeight(1);
  v->addWidget(rule);

  m_pReady = new QLabel(" ");
  m_pReady->setObjectName("reconReady");
  m_pReady->setWordWrap(true);
  v->addWidget(m_pReady);

  QHBoxLayout* h = new QHBoxLayout;
  h->setSpacing(8);
  h->addStretch(1);

  m_pRun = new QPushButton("Run Reconstruction");
  m_pRun->setObjectName("reconRun");
  m_pRun->setProperty("reconKind", "primary");
  m_pRun->setCursor(Qt::PointingHandCursor);
  m_pRun->setFixedHeight(36);
  connect(m_pRun, &QPushButton::clicked, this, &ReconPage::OnPrimary);
  h->addWidget(m_pRun);
  v->addLayout(h);
  return footer;
}

QWidget* ReconPage::BuildPane()
{
  // Parentless until MainWindow's layout adopts it (see ~ReconPage).
  QWidget* pane = new QWidget;
  pane->setObjectName("reconPane");
  QVBoxLayout* pl = new QVBoxLayout(pane);
  pl->setContentsMargins(16, 16, 16, 16);
  pl->setSpacing(16);

  QFrame* cmd = new QFrame;
  cmd->setObjectName("reconCard");
  QVBoxLayout* cl = new QVBoxLayout(cmd);
  cl->setContentsMargins(16, 14, 16, 12);
  cl->setSpacing(8);
  {
    QHBoxLayout* ch = new QHBoxLayout;
    ch->setSpacing(8);
    QLabel* ce = new QLabel("COMMAND");
    ce->setObjectName("reconSection");
    ch->addWidget(ce);
    ch->addStretch(1);
    m_pCopy = new QPushButton("Copy command");
    m_pCopy->setProperty("reconKind", "ghost");
    m_pCopy->setCursor(Qt::PointingHandCursor);
    m_pCopy->setFixedHeight(kControlHeight);
    {
      QFontMetrics fm(m_pCopy->font());
      m_pCopy->setMinimumWidth(qMax(fm.horizontalAdvance("Copy command"), fm.horizontalAdvance("Copied")) + 40);
    }
    connect(m_pCopy, &QPushButton::clicked, this, &ReconPage::OnCopyCommand);
    ch->addWidget(m_pCopy);
    cl->addLayout(ch);
  }
  m_pPreview = new QPlainTextEdit;
  m_pPreview->setObjectName("reconPreview");
  m_pPreview->setReadOnly(true);
  m_pPreview->setFrameShape(QFrame::NoFrame);
  m_pPreview->setLineWrapMode(QPlainTextEdit::WidgetWidth);
  m_pPreview->setWordWrapMode(QTextOption::WrapAnywhere);
  m_pPreview->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  m_pPreview->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  m_pPreview->setFont(MonoFont(font()));
  m_pPreview->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
  m_pPreview->viewport()->installEventFilter(this);
  cl->addWidget(m_pPreview);
  m_pPreviewNote = new QLabel(" ");
  m_pPreviewNote->setObjectName("reconCaption");
  m_pPreviewNote->setWordWrap(true);
  cl->addWidget(m_pPreviewNote);
  pl->addWidget(cmd, 0);

  QFrame* log = new QFrame;
  log->setObjectName("reconCard");
  m_pLogLayout = new QVBoxLayout(log);
  m_pLogLayout->setContentsMargins(16, 14, 16, 16);
  m_pLogLayout->setSpacing(8);
  {
    QHBoxLayout* lh = new QHBoxLayout;
    lh->setSpacing(8);
    QLabel* le = new QLabel("LOG");
    le->setObjectName("reconSection");
    lh->addWidget(le);
    lh->addStretch(1);
    m_pPhase = new QLabel(" ");
    m_pPhase->setObjectName("reconPhase");
    lh->addWidget(m_pPhase);
    m_pElapsed = new QLabel(" ");
    m_pElapsed->setObjectName("reconElapsed");
    m_pElapsed->setFont(MonoFont(font()));
    m_pElapsed->setMinimumWidth(48);
    m_pElapsed->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    lh->addWidget(m_pElapsed);
    m_pLogLayout->addLayout(lh);
  }
  m_pProgress = new QProgressBar;
  m_pProgress->setObjectName("reconProgress");
  m_pProgress->setRange(0, 0);
  m_pProgress->setTextVisible(false);
  m_pProgress->setFixedHeight(4);
  {
    QSizePolicy sp = m_pProgress->sizePolicy();
    sp.setRetainSizeWhenHidden(true);
    m_pProgress->setSizePolicy(sp);
  }
  m_pProgress->hide();
  m_pLogLayout->addWidget(m_pProgress);

  QWidget* status = new QWidget;
  status->setObjectName("reconStatusBlock");
  {
    QVBoxLayout* sv = new QVBoxLayout(status);
    m_pStatusLayout = sv;
    sv->setContentsMargins(0, 0, 0, 0);
    sv->setSpacing(6);
    // Grows only while the log is hidden (Idle), centring the idle line
    // vertically in the card; see SetState().
    sv->addStretch(0);
    m_pStatusLead = new QLabel(" ");
    m_pStatusLead->setObjectName("reconStatusLead");
    m_pStatusLead->setWordWrap(true);
    sv->addWidget(m_pStatusLead);
    m_pStatus = new QLabel(" ");
    m_pStatus->setObjectName("reconStatus");
    m_pStatus->setWordWrap(true);
    m_pStatus->setTextInteractionFlags(Qt::TextSelectableByMouse);
    sv->addWidget(m_pStatus);
    m_pStatusDetail = new QLabel(" ");
    m_pStatusDetail->setObjectName("reconStatusDetail");
    m_pStatusDetail->setFont(MonoFont(font()));
    m_pStatusDetail->setWordWrap(true);
    m_pStatusDetail->setTextInteractionFlags(Qt::TextSelectableByMouse);
    sv->addWidget(m_pStatusDetail);
    QHBoxLayout* bh = new QHBoxLayout;
    bh->setSpacing(8);
    m_pShowOutput = new QPushButton("Show output folder");
    m_pUnavailCopy = new QPushButton("Copy path");
    m_pCheckAgain = new QPushButton("Check again");
    foreach (QPushButton* b, QList<QPushButton*>() << m_pShowOutput << m_pUnavailCopy << m_pCheckAgain)
    {
      b->setProperty("reconKind", "ghost");
      b->setCursor(Qt::PointingHandCursor);
      b->setFixedHeight(kControlHeight);
      bh->addWidget(b);
      b->hide();
    }
    bh->addStretch(1);
    sv->addLayout(bh);
    connect(m_pShowOutput, &QPushButton::clicked, this, &ReconPage::OnShowOutput);
    connect(m_pUnavailCopy, &QPushButton::clicked, this, [this]() { QApplication::clipboard()->setText(m_strUnavailCopy); });
    connect(m_pCheckAgain, &QPushButton::clicked, this, &ReconPage::CheckAvailability);
    sv->addStretch(0);
  }
  m_pLogLayout->addWidget(status, 1);

  m_pLog = new QPlainTextEdit;
  m_pLog->setObjectName("reconLog");
  m_pLog->setReadOnly(true);
  m_pLog->setFrameShape(QFrame::NoFrame);
  m_pLog->setLineWrapMode(QPlainTextEdit::WidgetWidth);
  m_pLog->setMaximumBlockCount(kMaxLogBlocks);
  m_pLog->setFont(MonoFont(font()));
  m_pLog->setMinimumHeight(kLogMinHeight);
  m_pLog->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard);
  m_pLog->hide();
  m_pLogLayout->addWidget(m_pLog, 0);
  pl->addWidget(log, 1);
  return pane;
}

// ==== Form state ====

void ReconPage::ResetForm()
{
  const bool bWas = m_bLoading;
  m_bLoading = true;
  for (QMap<QString, QLineEdit*>::iterator it = m_edits.begin(); it != m_edits.end(); ++it)
    it.value()->clear();
  for (QMap<QString, QCheckBox*>::iterator it = m_checks.begin(); it != m_checks.end(); ++it)
    it.value()->setChecked(false);
  m_edits["output_name"]->setText(kOutputSubdir);
  m_strRefDefaultsDir.clear();
  m_pHemisphere->SetCurrentIndex(2, false);
  m_pFace->SetCurrentIndex(0, false);
  m_pOrder->SetCurrentIndex(0, false);
  m_pTissue->SetCurrentIndex(0, false);
  m_pRefNone->setChecked(true);
  m_bLoading = bWas;
}

void ReconPage::SetProject(const QString& strRoot)
{
  // The form is what the running command was started from; refilling it
  // underneath the run would make the preview and the log disagree.
  if (m_state == Running || m_state == Stopping)
    return;
  if (m_pSave->isActive())
  {
    m_pSave->stop();
    SaveSettings();
  }
  m_strRoot = strRoot.isEmpty() ? QString() : QDir(strRoot).absolutePath();
  m_bLoading = true;
  ResetForm();
  SetAdvancedOpen(false);
  SetExpertOpen(false);
  LoadSettings();
  FillDefaultReferencePaths();
  m_bLoading = false;
  if (m_bSettingsOutdated)
    SaveSettings();
  UpdateVisibility();
  UpdatePreview();
  if (m_state != Unavailable)
    SetState(Idle);
  CheckAvailability();
  OnValidate();
  // A run of this project started earlier in this session, perhaps from another
  // window: show it, still running or with its result.
  Attach(ReconRunner::ForRoot(m_strRoot));
}

// mri_3d_photo_recon runs SynthSR only when --ref_mri_synthsr names a file that
// doesn't exist yet (and only then does --low_field_synthsr apply). Default
// both reference outputs into the recon folder; only empty fields are filled.
void ReconPage::FillDefaultReferencePaths()
{
  const QString outDir = OutputDir();
  if (outDir.isEmpty())
    return;
  const bool bWas = m_bLoading;
  m_bLoading = true;
  // Defaults still pointing into the previous output folder follow a rename;
  // a path the user chose stays.
  const char* files[][2] = { {"ref_mri_synthsr", "synthsr.nii.gz"}, {"ref_mri_synthseg", "synthseg.nii.gz"} };
  for (const auto& f : files)
  {
    const QString cur = Text(f[0]);
    if (cur.isEmpty() || (!m_strRefDefaultsDir.isEmpty() && cur == QDir(m_strRefDefaultsDir).filePath(f[1])))
      m_edits[f[0]]->setText(QDir(outDir).filePath(f[1]));
  }
  m_strRefDefaultsDir = outDir;
  m_bLoading = bWas;
}

QVariantMap ReconPage::FormValues() const
{
  QVariantMap v;
  for (QMap<QString, QLineEdit*>::const_iterator it = m_edits.constBegin(); it != m_edits.constEnd(); ++it)
    v[it.key()] = it.value()->text().trimmed();
  for (QMap<QString, QCheckBox*>::const_iterator it = m_checks.constBegin(); it != m_checks.constEnd(); ++it)
    v[it.key()] = it.value()->isChecked();
  v["hemisphere"] = Hemisphere();
  v["photos_of_posterior_side"] = (m_pFace->CurrentIndex() == 1);
  v["order_posterior_to_anterior"] = (m_pOrder->CurrentIndex() == 1);
  v["fresh_tissue"] = (m_pTissue->CurrentIndex() == 1);
  v["ref_mode"] = ReferenceMode();
  return v;
}

void ReconPage::SetFormValues(const QVariantMap& v)
{
  const bool bWas = m_bLoading;
  m_bLoading = true;
  for (QVariantMap::const_iterator it = v.constBegin(); it != v.constEnd(); ++it)
  {
    const QString& k = it.key();
    if (m_edits.contains(k))
      m_edits[k]->setText(it.value().toString());
    else if (m_checks.contains(k))
      m_checks[k]->setChecked(it.value().toBool());
    else if (k == "hemisphere")
    {
      const QString h = it.value().toString();
      m_pHemisphere->SetCurrentIndex(h == "left" ? 0 : h == "right" ? 1 : 2, false);
    }
    else if (k == "photos_of_posterior_side")
      m_pFace->SetCurrentIndex(it.value().toBool() ? 1 : 0, false);
    else if (k == "order_posterior_to_anterior")
      m_pOrder->SetCurrentIndex(it.value().toBool() ? 1 : 0, false);
    else if (k == "fresh_tissue")
      m_pTissue->SetCurrentIndex(it.value().toBool() ? 1 : 0, false);
    else if (k == "ref_mode")
    {
      const QString m = it.value().toString();
      (m == "mesh" ? m_pRefMesh : m == "mri" ? m_pRefMri : m == "both" ? m_pRefBoth : m_pRefNone)->setChecked(true);
    }
  }
  m_bLoading = bWas;
  UpdateVisibility();
  UpdatePreview();
  OnValidate();
  if (!m_bLoading)
    m_pSave->start();
}

void ReconPage::SetAdvancedOpen(bool bOpen)
{
  m_pAdvanced->setVisible(bOpen);
  m_pAdvancedToggle->setText(bOpen ? QString::fromUtf8("\xe2\x96\xbe  Hide advanced settings")
                                   : QString::fromUtf8("\xe2\x96\xb8  Show advanced settings"));
}

void ReconPage::SetExpertOpen(bool bOpen)
{
  m_pExpert->setVisible(bOpen);
  m_pExpertToggle->setText(bOpen ? QString::fromUtf8("\xe2\x96\xbe  Hide expert weights")
                                 : QString::fromUtf8("\xe2\x96\xb8  Show expert weights"));
}

QString ReconPage::ReferenceMode() const
{
  if (m_pRefMesh->isChecked()) return "mesh";
  if (m_pRefMri->isChecked()) return "mri";
  if (m_pRefBoth->isChecked()) return "both";
  return "none";
}

QString ReconPage::Hemisphere() const
{
  switch (m_pHemisphere->CurrentIndex())
  {
  case 0: return "left";
  case 1: return "right";
  default: return "both";
  }
}

QString ReconPage::Text(const char* key) const
{
  QLineEdit* e = m_edits.value(key, nullptr);
  return e ? e->text().trimmed() : QString();
}

QString ReconPage::PhotoDir() const
{
  return m_strRoot.isEmpty() ? QString() : QDir(m_strRoot).filePath(kPhotoSubdir);
}

QString ReconPage::MaskDir() const
{
  return m_strRoot.isEmpty() ? QString() : QDir(m_strRoot).filePath(kMaskSubdir);
}

QString ReconPage::OutputDir() const
{
  const QString name = Text("output_name");
  return (m_strRoot.isEmpty() || name.isEmpty()) ? QString() : QDir(m_strRoot).filePath(name);
}

void ReconPage::OnFormEdited()
{
  if (m_bLoading)
    return;
  if (OutputDir() != m_strRefDefaultsDir)
    FillDefaultReferencePaths();
  UpdateVisibility();
  UpdatePreview();
  m_pValidate->start();
  m_pSave->start();
}

void ReconPage::UpdateVisibility()
{
  const QString mode = ReferenceMode();
  const bool bMesh = (mode == "mesh" || mode == "both");
  const bool bMri = (mode == "mri" || mode == "both");
  m_pMeshBlock->setVisible(bMesh);
  m_pMriBlock->setVisible(bMri);
  m_rows["stretch_factor_lr_mesh"]->setVisible(bMesh);
  m_rows["deform_recon_dir"]->setVisible(bMri);
  const bool bWeights = !Text("weights").isEmpty();
  m_edits["thickness_cap"]->setEnabled(bWeights);
  m_edits["thickness_cap"]->setToolTip(bWeights ? QString() : QString("Set a slab weights file first"));
  // The expert fields show what the command picks when they are left empty
  // (photo_reconstruction/utils.py, adjust_settings): it depends on the tissue
  // and the reference.
  const bool bFresh = m_pTissue->CurrentIndex() == 1;
  QMap<QString, QString> defaults;
  defaults["cp_spacing_2d"] = bFresh ? "10" : "25";
  defaults["cp_spacing_3d"] = bMri ? (bMesh ? "15" : "20") : (bMesh ? "not used with a mesh alone" : "10");
  defaults["k_lncc_mri"] = defaults["k_dice_mri"] = defaults["k_mesh_loss"] = "1";
  defaults["k_dif_slice_loss"] = "3";
  defaults["k_regularizer"] = "0.1";
  defaults["k_regularizer_nonlin"] = bFresh ? "1" : "3";
  defaults["k_regularizer_nonlin3d"] = bMri ? (bMesh ? "1.5" : "2") : (bMesh ? "0" : "1");
  defaults["k_regularizer_sz"] = "0.0005";
  for (int i = 0; kExpertKeys[i]; ++i)
    m_edits[kExpertKeys[i]]->setPlaceholderText(defaults[kExpertKeys[i]]);
  const bool bLowField = !Text("ref_mri").isEmpty() && !Text("ref_mri_synthsr").isEmpty();
  m_checks["low_field_synthsr"]->setEnabled(bLowField);
  m_checks["low_field_synthsr"]->setToolTip(bLowField ? QString() : QString("Needs an MRI scan and a SynthSR path"));
}

QStringList ReconPage::BuildArguments() const
{
  QStringList a;
  auto emitIf = [&a, this](const char* k)
  {
    const QString v = Text(k);
    if (!v.isEmpty())
      a << QString("--") + k << v;
  };
  // The four the parser requires with no default are always emitted; an
  // empty one is caught by Validate() before Run is enabled.
  a << "--input_photo_dir" << PhotoDir();
  a << "--input_segmentation_dir" << MaskDir();
  a << "--output_directory" << OutputDir();
  a << "--hemisphere" << Hemisphere();
  emitIf("slice_thickness");
  a << "--photo_resolution" << kDefaultPhotoResolution;
  if (m_pFace->CurrentIndex() == 1) a << "--photos_of_posterior_side";
  if (m_pOrder->CurrentIndex() == 1) a << "--order_posterior_to_anterior";
  if (m_pTissue->CurrentIndex() == 1) a << "--fresh_tissue";
  const QString mode = ReferenceMode();
  const bool bMesh = (mode == "mesh" || mode == "both");
  const bool bMri = (mode == "mri" || mode == "both");
  if (bMesh)
  {
    emitIf("ref_mesh");
    emitIf("mesh_reorient_with_indices");
    emitIf("stretch_factor_lr_mesh");
  }
  if (bMri)
  {
    emitIf("ref_mri");
    emitIf("ref_mri_synthseg");
    emitIf("ref_mri_synthsr");
    if (m_checks["low_field_synthsr"]->isChecked() && !Text("ref_mri").isEmpty() && !Text("ref_mri_synthsr").isEmpty())
      a << "--low_field_synthsr";
    emitIf("deform_recon_dir");
  }
  emitIf("input_roi_dir");
  emitIf("weights");
  // The cap only means something with weights (utils.py: "only used with
  // --weights"), so it is not sent on its own.
  if (!Text("weights").isEmpty())
    emitIf("thickness_cap");
  emitIf("initial_stretch_factor_lr_photos");
  if (m_checks["no_z_stretch"]->isChecked()) a << "--no_z_stretch";
  if (m_checks["equalize_images"]->isChecked()) a << "--equalize_images";
  if (m_checks["skip_bfgs"]->isChecked()) a << "--skip_bfgs";
  // Empty reads as 1 and is passed: FreeSurfer 8.1 actually defaults --threads
  // to -1 (all cores) despite its help text.
  a << "--threads" << (Text("threads").isEmpty() ? QString("1") : Text("threads"));
  emitIf("gpu");
  for (int i = 0; kExpertKeys[i]; ++i)
    emitIf(kExpertKeys[i]);
  return a;
}

QString ReconPage::CommandPreview(bool bMultiline) const
{
  QString prog = Program();
  prog = prog.isEmpty() ? QString("mri_3d_photo_recon") : QFileInfo(prog).fileName();
  const QStringList args = BuildArguments();
  if (!bMultiline)
  {
    QStringList q;
    q << ShellQuote(prog);
    foreach (const QString& s, args)
      q << ShellQuote(s);
    return q.join(" ");
  }
  QStringList lines;
  lines << ShellQuote(prog);
  for (int i = 0; i < args.size(); ++i)
  {
    QString line = "  " + args[i];
    if (i + 1 < args.size() && !args[i + 1].startsWith("--"))
      line += " " + ShellQuote(args[++i]);
    lines << line;
  }
  return lines.join(" \\\n");
}

void ReconPage::UpdatePreview()
{
  const QString str = CommandPreview(true);
  if (m_pPreview->toPlainText() != str)
    m_pPreview->setPlainText(str);
  // At least seven lines (so the Log card only moves when a flag is added),
  // at most what the text needs, counting wrapped lines.
  QTimer::singleShot(0, this, [this]()
  {
    int nLines = 0;
    for (QTextBlock b = m_pPreview->document()->begin(); b.isValid(); b = b.next())
      nLines += qMax(1, b.layout() ? b.layout()->lineCount() : 1);
    nLines = qMax(kPreviewMinLines, nLines);
    const QFontMetrics fm(m_pPreview->font());
    const int h = nLines * fm.lineSpacing() + int(m_pPreview->document()->documentMargin() * 2)
                  + 2 * m_pPreview->frameWidth() + 8;
    if (m_pPreview->height() != h)
      m_pPreview->setFixedHeight(h);
  });
}

QStringList ReconPage::ListPhotos(const QString& strDir, QString* pExt, int* pOther) const
{
  // The command's own glob order (utils.py get_photo_and_seg_lists): the
  // first pattern that matches anything wins, lower case before upper.
  static const char* exts[] = {"jpg", "jpeg", "tif", "tiff", "png", "JPG", "JPEG", "TIF", "TIFF", "PNG", nullptr};
  QStringList found;
  QString used;
  int other = 0;
  if (!strDir.isEmpty())
  {
    QDir d(strDir);
    for (int i = 0; exts[i]; ++i)
    {
      const QStringList l = d.entryList(QStringList() << QString("*.") + exts[i],
                                        QDir::Files | QDir::CaseSensitive, QDir::Name);
      if (l.isEmpty())
        continue;
      if (found.isEmpty())
      {
        found = l;
        used = exts[i];
      }
      else
        other += l.size();
    }
  }
  if (pExt) *pExt = used;
  if (pOther) *pOther = other;
  return found;
}

QStringList ReconPage::ListMasks(const QString& strDir) const
{
  static const char* exts[] = {"mat", "npy", "MAT", "NPY", nullptr};
  if (strDir.isEmpty())
    return QStringList();
  QDir d(strDir);
  for (int i = 0; exts[i]; ++i)
  {
    const QStringList l = d.entryList(QStringList() << QString("*.") + exts[i],
                                      QDir::Files | QDir::CaseSensitive, QDir::Name);
    if (!l.isEmpty())
      return l;
  }
  return QStringList();
}

QList<ReconPage::Issue> ReconPage::Validate() const
{
  QList<Issue> out;
  auto err = [&out](const QString& k, const QString& m) { Issue i; i.key = k; i.message = m; i.blocking = true; out << i; };
  auto warn = [&out](const QString& k, const QString& m) { Issue i; i.key = k; i.message = m; i.blocking = false; out << i; };

  // inputs: this project's own folders
  const QString photoDir = PhotoDir();
  const QString maskDir = MaskDir();
  QString ext;
  int other = 0;
  if (m_strRoot.isEmpty()) err("inputs", "No project is open");
  else
  {
    const QStringList photos = QDir(photoDir).exists() ? ListPhotos(photoDir, &ext, &other) : QStringList();
    const QStringList masks = QDir(maskDir).exists() ? ListMasks(maskDir) : QStringList();
    if (photos.isEmpty()) err("inputs", "This project has no corrected photos. Go back to Correction");
    else if (masks.isEmpty()) err("inputs", "This project has no saved slabs. Go back to Connected components");
    else if (masks.size() != photos.size())
      err("inputs", QString::fromUtf8("%1 photos but %2 masks \xe2\x80\x94 every photo needs its slabs saved in Connected components")
          .arg(photos.size()).arg(masks.size()));
    else if (other > 0) warn("inputs", QString("Only the %1 .%2 files will be used").arg(photos.size()).arg(ext));
  }

  // output: a plain folder name inside the project
  const QString name = Text("output_name");
  static const QRegularExpression badName("[/\\\\:]");
  if (name.isEmpty()) err("output_name", "Enter a folder name");
  else if (name.contains(badName) || name == "." || name == ".." || name.startsWith('.'))
    err("output_name", "Use a plain folder name, without slashes");
  else if (QStringList({kPhotoSubdir, kMaskSubdir, "mask", "photos"}).contains(name, Qt::CaseInsensitive))
    err("output_name", QString::fromUtf8("That folder holds the project\xe2\x80\x99s own files. Choose another name"));
  else if (!m_strRoot.isEmpty())
  {
    const QString outDir = OutputDir();
    QDir d(outDir);
    if (QFileInfo(outDir).exists() && !QFileInfo(outDir).isDir()) err("output_name", "A file has this name. Choose another");
    else if (!QFileInfo(d.exists() ? outDir : m_strRoot).isWritable()) err("output_name", "Cannot write here. Check the project folder's permissions");
    else if (d.exists() && !d.entryList(QStringList() << "photo_recon.*.mgz", QDir::Files).isEmpty())
      warn("output_name", "Existing reconstruction will be overwritten");
  }

  // slabs
  double v = 0;
  const QString th = Text("slice_thickness");
  if (th.isEmpty()) err("slice_thickness", "Enter the slab thickness in millimetres");
  else if (!ParseNumber(th, &v) || v <= 0) err("slice_thickness", "Enter a number above 0");
  else if (v < 2 || v > 20) warn("slice_thickness", "Thickness is usually between 2 and 20 mm");

  // reference
  const QString mode = ReferenceMode();
  const bool bMesh = (mode == "mesh" || mode == "both");
  const bool bMri = (mode == "mri" || mode == "both");
  if (bMesh)
  {
    const QString f = Text("ref_mesh");
    if (f.isEmpty()) err("ref_mesh", "Choose the mesh file");
    else if (!QFileInfo(f).isFile()) err("ref_mesh", "That file was not found");
    // The command refuses a mesh without them (utils.py: "You provided mesh
    // but no indices to reorient it").
    const QString idx = Text("mesh_reorient_with_indices");
    if (idx.isEmpty()) err("mesh_reorient_with_indices", "A mesh needs three reorientation vertex indices");
    else
    {
      static const QRegularExpression three("^\\s*\\d+\\s*,\\s*\\d+\\s*,\\s*\\d+\\s*$");
      if (!three.match(idx).hasMatch())
        err("mesh_reorient_with_indices", "Enter three whole numbers separated by commas");
    }
    const QString st = Text("stretch_factor_lr_mesh");
    if (!st.isEmpty() && (!ParseNumber(st, &v) || v <= 0)) err("stretch_factor_lr_mesh", "Enter a number above 0");
  }
  if (bMri)
  {
    const QString f = Text("ref_mri");
    if (f.isEmpty()) err("ref_mri", "Choose the MRI file");
    else if (!QFileInfo(f).isFile()) err("ref_mri", "That file was not found");
    // An existing file is reused as it is, even if it came from another MRI
    // or was made without the low-field model, so say so.
    const char* outputs[] = {"ref_mri_synthsr", "ref_mri_synthseg"};
    for (const char* k : outputs)
    {
      const QString o = Text(k);
      if (!o.isEmpty() && QFileInfo(o).isFile())
        warn(k, "This file exists and will be used as it is. Delete it to compute it again");
      else if (!o.isEmpty() && !QFileInfo(QFileInfo(o).absolutePath()).isDir() &&
               QFileInfo(QFileInfo(o).absolutePath()).absoluteFilePath() != QDir(OutputDir()).absolutePath())
        err(k, "Folder not found");
    }
    const QString dr = Text("deform_recon_dir");
    if (!dr.isEmpty() && !QDir(dr + "/surf").exists())
      err("deform_recon_dir", "Expected a FreeSurfer subject folder with a surf/ subfolder");
  }

  // advanced
  const QString roi = Text("input_roi_dir");
  if (!roi.isEmpty() && !QDir(roi).exists()) err("input_roi_dir", "Folder not found");
  const QString w = Text("weights");
  if (!w.isEmpty() && !QFileInfo(w).isFile()) err("weights", "That file was not found");
  const QString cap = Text("thickness_cap");
  if (!w.isEmpty() && !cap.isEmpty() && (!ParseNumber(cap, &v) || v <= 0)) err("thickness_cap", "Enter a number above 0");
  const QString sp = Text("initial_stretch_factor_lr_photos");
  if (!sp.isEmpty() && (!ParseNumber(sp, &v) || v <= 0)) err("initial_stretch_factor_lr_photos", "Enter a number above 0");
  int n = 0;
  const QString thr = Text("threads");
  if (!thr.isEmpty() && (!ParseInt(thr, &n) || (n < 1 && n != -1)))
    err("threads", QString::fromUtf8("Enter a whole number of threads, or \xe2\x88\x92" "1 for all cores"));
  const QString gpu = Text("gpu");
  if (!gpu.isEmpty() && (!ParseInt(gpu, &n) || n < 0)) err("gpu", "Enter the index of the card (0 for the first)");
  for (int i = 0; kExpertKeys[i]; ++i)
  {
    const QString x = Text(kExpertKeys[i]);
    if (!x.isEmpty() && !ParseNumber(x, &v)) err(kExpertKeys[i], "Enter a number");
  }

  // Last, so any other problem with the same field is the one shown.
  const QStringList emitted = BuildArguments();
  if (m_strRoot.contains(' '))
    warn("inputs", QString::fromUtf8("The project folder\xe2\x80\x99s path has a space in it. FreeSurfer\xe2\x80\x99s launcher splits it there, so the run may fail"));
  for (int i = 0; kPathKeys[i]; ++i)
  {
    const QString p = Text(kPathKeys[i]);
    if (!p.isEmpty() && p.contains(QRegularExpression("\\s")) && emitted.contains(p))
      warn(kPathKeys[i], QString::fromUtf8("This path has a space in it. FreeSurfer\xe2\x80\x99s launcher splits it there, so the run may fail"));
  }
  return out;
}

void ReconPage::SetCaptionText(QLabel* l, const QString& strText, const char* role)
{
  if (!l)
    return;
  if (l->text() != strText)
    l->setText(strText);
  if (l->isHidden() != strText.isEmpty())
    l->setVisible(!strText.isEmpty());
  SetKind(l, "level", role);
}

void ReconPage::ApplyIssues(const QList<Issue>& issues)
{
  QMap<QString, bool> seen;
  QString firstBlocking, firstBlockingKey;
  foreach (const Issue& i, issues)
  {
    if (seen.contains(i.key))
      continue;
    seen[i.key] = true;
    SetCaptionText(m_captions.value(i.key, nullptr), i.message, i.blocking ? "error" : "warning");
    if (i.blocking && firstBlocking.isEmpty())
    {
      firstBlocking = i.message;
      firstBlockingKey = i.key;
    }
  }
  for (QMap<QString, QString>::const_iterator it = m_help.constBegin(); it != m_help.constEnd(); ++it)
    if (!seen.contains(it.key()))
      SetCaptionText(m_captions.value(it.key(), nullptr), it.value(), "");

  const bool bValid = firstBlocking.isEmpty();
  if (m_state == Idle || m_state == Succeeded || m_state == Failed || m_state == Stopped)
  {
    m_pRun->setEnabled(bValid);
    const QString cap = (firstBlockingKey == "slice_thickness" && Text("slice_thickness").isEmpty())
                        ? QString("Enter the slab thickness to run.") : firstBlocking;
    m_pReady->setText(bValid ? QString(" ") : cap);
  }
  SetCaptionText(m_pPreviewNote, bValid ? QString(" ") : QString("Fix the fields marked on the left to run."),
                 bValid ? "" : "warning");
}

void ReconPage::OnValidate()
{
  RefreshInventory();
  ApplyIssues(Validate());
}

void ReconPage::RefreshInventory()
{
  const int nPhotos = ListPhotos(PhotoDir()).size();
  const int nMasks = ListMasks(MaskDir()).size();
  m_pInventory->setText(QString::fromUtf8("%1 photo%2 \xc2\xb7 %3 mask%4")
                        .arg(nPhotos).arg(nPhotos == 1 ? "" : "s")
                        .arg(nMasks).arg(nMasks == 1 ? "" : "s"));
}

void ReconPage::SetNarrow(bool bNarrow)
{
  if (bNarrow == m_bNarrow)
    return;
  m_bNarrow = bNarrow;
  foreach (ReconFormRow* r, m_listRows)
    r->SetNarrow(bNarrow);
}

bool ReconPage::eventFilter(QObject* o, QEvent* e)
{
  if (e->type() == QEvent::Resize)
  {
    // Decided from the viewport, not the form: the form cannot shrink below
    // its own minimum, so measuring it would never report "narrow".
    if (o == m_pScroll->viewport())
      SetNarrow(static_cast<QResizeEvent*>(e)->size().width() < kWideMin);
    else if (m_pPreview && o == m_pPreview->viewport())
      QTimer::singleShot(0, this, [this]() { UpdatePreview(); });
  }
  return QWidget::eventFilter(o, e);
}

void ReconPage::showEvent(QShowEvent* e)
{
  QWidget::showEvent(e);
  CheckAvailability();
}

// ==== Settings: per project root, paths inside the project stored relative ====

QString ReconPage::SettingsKey() const
{
  return SettingsKeyFor(m_strRoot);
}

QString ReconPage::SettingsKeyFor(const QString& strRoot)
{
  // SetProject() keeps the root as QDir::absolutePath().
  const QString strAbs = strRoot.isEmpty() ? QString() : QDir(strRoot).absolutePath();
  const QByteArray h = QCryptographicHash::hash(strAbs.toUtf8(), QCryptographicHash::Md5).toHex();
  return "Recon/" + QString::fromLatin1(h);
}

void ReconPage::SaveSettings()
{
  if (m_strRoot.isEmpty() || m_bLoading)
    return;
  QVariantMap v = FormValues();
  for (QVariantMap::iterator it = v.begin(); it != v.end(); ++it)
  {
    if (it.value().userType() == QMetaType::QString)
    {
      const QString s = it.value().toString();
      if (s == m_strRoot || s.startsWith(m_strRoot + "/"))
        it.value() = QString("$PROJECT") + s.mid(m_strRoot.length());
    }
  }
  QSettings s;
  s.setValue(SettingsKey() + "/settings", v);
  s.setValue(SettingsKey() + "/advanced_open", m_pAdvanced->isVisibleTo(m_pForm));
  s.setValue(SettingsKey() + "/expert_open", m_pExpert->isVisibleTo(m_pAdvanced));
}

void ReconPage::LoadSettings()
{
  if (m_strRoot.isEmpty())
    return;
  QSettings s;
  QVariantMap saved = s.value(SettingsKey() + "/settings").toMap();
  // Older versions saved the input and output paths, the photo resolution and
  // the GPU. An output folder inside the project becomes the folder name; the
  // rest is dropped, and the stored settings are rewritten with only the
  // current fields (SetProject() saves once loading is done).
  if (!saved.contains("output_name") && saved.value("output_directory").toString().startsWith("$PROJECT/"))
  {
    const QString name = saved.value("output_directory").toString().mid(9);
    if (!name.isEmpty() && !name.contains('/'))
      saved["output_name"] = name;
  }
  const QVariantMap current = FormValues();
  m_bSettingsOutdated = false;
  QVariantMap v;
  for (QVariantMap::const_iterator it = saved.constBegin(); it != saved.constEnd(); ++it)
  {
    if (!current.contains(it.key()))
      m_bSettingsOutdated = true;
    else if (it.value().userType() == QMetaType::Bool)
      v[it.key()] = it.value();
    else
    {
      QString val = it.value().toString();
      if (val.startsWith("$PROJECT"))
        val = m_strRoot + val.mid(8);
      v[it.key()] = val;
    }
  }
  if (!v.isEmpty())
    SetFormValues(v);
  SetAdvancedOpen(s.value(SettingsKey() + "/advanced_open", false).toBool());
  SetExpertOpen(s.value(SettingsKey() + "/expert_open", false).toBool());
}

// ==== Availability and state ====

QString ReconPage::Program() const
{
  const QString env = QProcessEnvironment::systemEnvironment().value("DISSECTION_PHOTO_RECON_CMD");
  if (!env.isEmpty())
    return env;
  const QString home = FreesurferHome();
  if (home.isEmpty())
    return QString();
  return QDir(home).filePath("bin/mri_3d_photo_recon");
}

bool ReconPage::IsRunning() const
{
  return m_state == Running || m_state == Stopping || (m_pRunner && m_pRunner->IsActive());
}

void ReconPage::CheckAvailability()
{
  if (m_state == Running || m_state == Stopping)
    return;
  const QString home = FreesurferHome();
  const QString prog = Program();
  const bool bOverridden = Overridden();
  QString lead, body, detail, copy;
  if (prog.isEmpty())
  {
    lead = "FreeSurfer is not set up";
    body = "Open a terminal, set up FreeSurfer, then start Dissection Photo from that terminal.";
    detail = "source $FREESURFER_HOME/SetUpFreeSurfer.sh";
  }
  else if (!QFileInfo(prog).isExecutable())
  {
    lead = QString::fromUtf8("FreeSurfer\xe2\x80\x99s reconstruction command was not found");
    body = bOverridden ? QString("The command this app was told to use does not exist.")
                       : QString("This FreeSurfer installation has no mri_3d_photo_recon. It ships with FreeSurfer 7.4 and later.");
    detail = bOverridden ? QString("Expected %1").arg(prog) : QString("Looked in %1").arg(home);
    copy = prog;
  }
  else if (!bOverridden && !QFileInfo(QDir(home).filePath("bin/fspython")).exists())
  {
    lead = QString::fromUtf8("FreeSurfer\xe2\x80\x99s Python runtime (fspython) was not found");
    body = "This FreeSurfer installation is incomplete. Reinstall FreeSurfer or ask whoever set it up.";
    detail = QString("Expected %1/bin/fspython").arg(home);
    copy = QDir(home).filePath("bin/fspython");
  }
  else if (!bOverridden && !QFileInfo(QDir(home).filePath("models/photo_imputation_unet.pth")).exists())
  {
    lead = "One model file is missing";
    body = QString::fromUtf8("Download photo_imputation_unet.pth and place it in the models folder, then check again "
                             "\xe2\x80\x94 no restart needed.");
    detail = QDir(home).filePath("models");
    copy = detail;
  }

  if (lead.isEmpty())
  {
    if (m_state == Unavailable)
      SetState(Idle);
    return;
  }
  m_strUnavailCopy = copy;
  m_pStatusLead->setText(lead);
  m_pStatus->setText(body);
  m_pStatusDetail->setText(detail);
  m_pUnavailCopy->setText(copy == QDir(home).filePath("models") ? "Copy destination path" : "Copy path");
  SetState(Unavailable);
}

void ReconPage::SetState(State s)
{
  m_state = s;
  const bool bRunning = (s == Running || s == Stopping);
  m_pForm->setEnabled(!bRunning);
  m_pBack->setEnabled(!bRunning);
  m_pRun->setEnabled(true);
  SetKind(m_pRun, "reconKind", s == Running ? "danger" : "primary");
  m_pProgress->setVisible(bRunning);
  m_pShowOutput->setVisible(s == Succeeded);
  m_pUnavailCopy->setVisible(s == Unavailable && !m_strUnavailCopy.isEmpty());
  m_pCheckAgain->setVisible(s == Unavailable);
  m_pStatusDetail->setVisible(s == Unavailable && !m_pStatusDetail->text().trimmed().isEmpty());
  const QString strElapsed = ElapsedText(m_nFinalMs);

  switch (s)
  {
  case Idle:
    m_pRun->setText("Run Reconstruction");
    m_pStatusLead->hide();
    m_pStatus->setText(QString::fromUtf8("Nothing has run yet. The command\xe2\x80\x99s output will appear here."));
    m_pStatus->setAlignment(Qt::AlignCenter);
    SetCaptionText(m_pPhase, " ", "");
    m_pElapsed->setText(" ");
    break;
  case Running:
    m_pRun->setText("Stop");
    m_pStatusLead->show();
    m_pStatusLead->setText(QString::fromUtf8("Reconstructing\xe2\x80\xa6"));
    m_pStatus->setText("Checking inputs");
    m_pStatus->setAlignment(Qt::AlignLeft);
    SetCaptionText(m_pPhase, QString::fromUtf8("Starting\xe2\x80\xa6"), "");
    m_pElapsed->setText("00:00");
    m_pReady->setText("You can leave this page open; the run continues.");
    break;
  case Stopping:
    m_pRun->setText(QString::fromUtf8("Stopping\xe2\x80\xa6"));
    m_pRun->setEnabled(false);
    SetCaptionText(m_pPhase, QString::fromUtf8("Stopping\xe2\x80\xa6"), "");
    break;
  case Succeeded:
    m_pRun->setText("Run Reconstruction");
    m_pStatusLead->show();
    m_pStatusLead->setText("Reconstruction finished");
    m_pStatus->setAlignment(Qt::AlignLeft);
    {
      QDir d(OutputDir());
      const int nFiles = d.entryList(QStringList() << "photo_recon*.mgz" << "mri.deformed*.mgz" << "*.deformed"
                                     << "registered.surf", QDir::Files).size();
      m_pStatus->setText(QString::fromUtf8("Saved to %1 \xc2\xb7 %2 file%3 \xc2\xb7 %4")
                         .arg(d.dirName()).arg(nFiles).arg(nFiles == 1 ? "" : "s").arg(strElapsed));
    }
    SetCaptionText(m_pPhase, QString::fromUtf8("Done \xc2\xb7 ") + strElapsed, "success");
    m_pElapsed->setText(" ");
    break;
  case Failed:
    m_pRun->setText("Try again");
    m_pStatusLead->show();
    m_pStatusLead->setText("Reconstruction did not finish");
    m_pStatus->setAlignment(Qt::AlignLeft);
    SetCaptionText(m_pPhase, QString::fromUtf8("Failed \xc2\xb7 ") + strElapsed, "error");
    m_pElapsed->setText(" ");
    break;
  case Stopped:
    m_pRun->setText("Run Reconstruction");
    m_pStatusLead->show();
    m_pStatusLead->setText("Reconstruction stopped");
    m_pStatus->setText("Nothing was saved. You can run it again later.");
    m_pStatus->setAlignment(Qt::AlignLeft);
    SetCaptionText(m_pPhase, QString::fromUtf8("Stopped \xc2\xb7 ") + strElapsed, "");
    m_pElapsed->setText(" ");
    break;
  case Unavailable:
    m_pRun->setText("Run Reconstruction");
    m_pRun->setEnabled(false);
    m_pStatusLead->show();
    m_pStatus->setAlignment(Qt::AlignLeft);
    SetCaptionText(m_pPhase, " ", "");
    m_pElapsed->setText(" ");
    m_pReady->setText(m_pStatusLead->text());
    break;
  }
  SetKind(m_pStatusLead, "level", s == Failed ? "error" : "");

  // Idle and Unavailable have no output: the status block fills the card.
  // Other states show the log below a compact status block.
  const bool bShowLog = !(s == Idle || s == Unavailable);
  m_pLog->setVisible(bShowLog);
  m_pLogLayout->setStretchFactor(m_pLog, bShowLog ? 1 : 0);
  m_pLogLayout->setStretchFactor(m_pStatus->parentWidget(), bShowLog ? 0 : 1);
  // Idle: equal stretch above and below centres the line. Otherwise the
  // text hugs the top of its block.
  m_pStatusLayout->setStretch(0, s == Idle ? 1 : 0);
  m_pStatusLayout->setStretch(m_pStatusLayout->count() - 1, bShowLog ? 0 : 1);

  if (!bRunning && s != Unavailable)
    ApplyIssues(Validate());
}

// ==== Running ====

void ReconPage::OnPrimary()
{
  if (m_state == Running)
  {
    QMessageBox box(window());
    box.setWindowTitle("Stop the reconstruction?");
    box.setText("Stop the reconstruction?");
    box.setInformativeText("Progress so far is lost. You can run it again later.");
    QPushButton* stop = box.addButton("Stop", QMessageBox::DestructiveRole);
    QPushButton* keep = box.addButton("Keep running", QMessageBox::RejectRole);
    box.setDefaultButton(keep);
    box.exec();
    if (box.clickedButton() == stop)
      StopWithoutConfirm();
    return;
  }
  if (m_state == Stopping || m_state == Unavailable)
    return;
  const QString prog = Program();
  if (prog.isEmpty() || !QFileInfo(prog).isExecutable())
  {
    CheckAvailability();
    return;
  }
  // Checked here as well as by the disabled button: a click can arrive in
  // the 150 ms before a debounced validation catches up with an edit.
  const QList<Issue> issues = Validate();
  foreach (const Issue& i, issues)
    if (i.blocking)
    {
      ApplyIssues(issues);
      return;
    }

  // The command writes here; creating it first turns "the folder does not
  // exist" from a failure several minutes in into nothing at all.
  QDir().mkpath(OutputDir());

  m_pSave->stop();
  SaveSettings();
  m_nFinalMs = 0;
  Attach(ReconRunner::Start(m_strRoot, prog, BuildArguments(), CommandPreview(false)));
}

void ReconPage::StopWithoutConfirm()
{
  if (m_state != Running || !m_pRunner || !m_pRunner->IsActive())
    return;
  SetState(Stopping);
  m_pRunner->Stop();
}

void ReconPage::Attach(ReconRunner* r)
{
  if (m_pRunner)
    disconnect(m_pRunner, nullptr, this, nullptr);
  m_pRunner = r;
  m_pTick->stop();
  if (!r)
    return;
  connect(r, &ReconRunner::LineAdded, this, [this](const QString& strLine) { m_pLog->appendPlainText(strLine); });
  connect(r, &ReconRunner::PhaseChanged, this, [this](const QString& strPhase)
  {
    if (m_state == Running)
      SetCaptionText(m_pPhase, strPhase, "");
  });
  connect(r, &ReconRunner::Finished, this, &ReconPage::OnRunFinished);
  // Everything so far in one go: a long log appended line by line is slow.
  m_pLog->setPlainText(r->Lines().join("\n"));
  m_pLog->moveCursor(QTextCursor::End);
  if (!r->IsActive())
  {
    OnRunFinished();
    return;
  }
  SetState(Running);
  if (r->IsStopping())
    SetState(Stopping);
  else if (!r->Phase().isEmpty())
    SetCaptionText(m_pPhase, r->Phase(), "");
  m_pTick->start();
  OnTick();
}

void ReconPage::OnRunFinished()
{
  m_pTick->stop();
  if (!m_pRunner)
    return;
  m_nFinalMs = m_pRunner->ElapsedMs();
  switch (m_pRunner->GetOutcome())
  {
  case ReconRunner::Stopped:
    SetState(Stopped);
    break;
  case ReconRunner::Succeeded:
    SetState(Succeeded);
    break;
  case ReconRunner::Failed:
    if (m_pRunner->FailedToStart())
      m_pStatus->setText("The reconstruction tool could not be started.");
    else if (m_pRunner->Crashed())
      m_pStatus->setText("The reconstruction tool crashed. No volume was written.");
    else
      m_pStatus->setText(QString("No volume was written (exit code %1). The most common causes are a wrong hemisphere or photo "
                                 "order, or too little memory. The log below says where it stopped.").arg(m_pRunner->ExitCode()));
    SetState(Failed);
    break;
  case ReconRunner::Active:
    break;
  }
}

void ReconPage::OnCopyCommand()
{
  QApplication::clipboard()->setText(CommandPreview(false));
  m_pCopy->setText("Copied");
  m_pCopied->start();
}

void ReconPage::OnShowOutput()
{
  const QString str = OutputDir();
  if (!str.isEmpty() && QDir(str).exists())
    QDesktopServices::openUrl(QUrl::fromLocalFile(str));
}

QString ReconPage::ElapsedText(qint64 ms) const
{
  const qint64 s = ms / 1000;
  const qint64 h = s / 3600, m = (s / 60) % 60, sec = s % 60;
  if (h > 0)
    return QString("%1:%2:%3").arg(h).arg(m, 2, 10, QChar('0')).arg(sec, 2, 10, QChar('0'));
  return QString("%1:%2").arg(m, 2, 10, QChar('0')).arg(sec, 2, 10, QChar('0'));
}

void ReconPage::OnTick()
{
  const qint64 ms = m_pRunner ? m_pRunner->ElapsedMs() : 0;
  m_pElapsed->setText(ElapsedText(ms));
  if (m_state == Running || m_state == Stopping)
  {
    const qint64 s = ms / 1000;
    const QString forText = s < 60 ? QString("%1 s").arg(s) : QString("%1 min").arg(s / 60);
    m_pStatus->setText(QString::fromUtf8("Running for %1 \xc2\xb7 this usually takes 10 to 40 minutes").arg(forText));
  }
}

// ==== Theme ====

void ReconPage::SetDarkMode(bool bDark)
{
  const QString strCaption = bDark ? "#a8a8a8" : "#6a6a6a";
  const QString strBody    = bDark ? "#f0f0f0" : "#262626";
  const QString strRule    = bDark ? "#3f3f3f" : "#d5d5d5";
  const QString strFieldBg = bDark ? "#333333" : "#ffffff";
  const QString strFieldBd = bDark ? "#5f5f5f" : "#c0c0c0";
  const QString strFocus   = bDark ? kPurpleInkDark : kPurple;
  const QString strCardBg  = bDark ? "#1e1e1e" : "#f7f7f7";
  const QString strCardBd  = bDark ? "#333333" : "#b3b3b3";
  const QString strDanger  = bDark ? "#e08a8d" : "#a2293c";
  const QString strWarning = bDark ? "#e0b252" : "#8a6100";
  const QString strSuccess = bDark ? "#8bd48b" : "#3f7f3f";
  const QString strGhostBd = bDark ? "#777777" : "#999999";
  const QString strGhostTx = bDark ? "#dddddd" : "#555555";
  const QString strGhostHv = bDark ? "#4a4a4a" : "#e4e4e4";
  const QString strDisTx   = bDark ? "#6a6a6a" : "#bbbbbb";
  const QString strDisBd   = bDark ? "#3f3f3f" : "#dddddd";
  const QString strDisBg   = bDark ? "#3a3a3a" : "#999999";
  const QString strSel     = kPurple;
  const QString strDisRun  = bDark ? "#7a7a7a" : "#ffffff";
  const QString strTrack   = bDark ? "#333333" : "#d8d8d8";

  // Heights are set in the sheet, not with setFixedHeight(): a re-polish (e.g.
  // on theme switch) restores pre-polish size constraints, undoing fixed
  // heights and breaking pill radii.
  //
  // The widget's own sheet beats the app sheet, so the bare QWidget rule stops
  // style.qss's app-wide background from banding every child. Surfaces are
  // declared by name below.
  QString css;
  css += "QWidget {background: transparent;}";
  css += QString("QLabel {color: %1; font-size: 13px;}").arg(strBody);
  css += QString("QCheckBox, QRadioButton {color: %1; font-size: 13px; spacing: 6px;}"
                 "QCheckBox:disabled, QRadioButton:disabled {color: %2;}").arg(strBody, strCaption);
  // Colours only: style.qss keeps every state's indicator the same size and
  // border width, which is what stops the label shifting on toggle.
  css += QString("QCheckBox::indicator:checked, QRadioButton::indicator:checked"
                 " {background-color: %1; border: 1px solid %1;}").arg(strFocus);
  css += QString("#reconSubtitle {font-size: 14px; color: %1;}").arg(strBody);
  css += QString("#reconSection {font-size: 12px; font-weight: bold; color: %1;}").arg(strCaption);
  css += QString("#reconLabel {font-size: 13px; color: %1;}").arg(strBody);
  css += QString("#reconLabel:disabled {color: %1;}").arg(strCaption);
  css += QString("#reconUnit, #reconInventory, #reconReady, #reconPhase, #reconElapsed, #reconStatusDetail"
                 " {font-size: 12px; color: %1;}").arg(strCaption);
  css += QString("#reconCaption {font-size: 12px; color: %1;}").arg(strCaption);
  css += QString("#reconCaption[level=\"error\"], #reconPhase[level=\"error\"], #reconStatusLead[level=\"error\"]"
                 " {color: %1;}").arg(strDanger);
  css += QString("#reconCaption[level=\"warning\"] {color: %1;}").arg(strWarning);
  css += QString("#reconPhase[level=\"success\"] {color: %1;}").arg(strSuccess);
  css += QString("#reconStatusLead {font-size: 14px; font-weight: bold; color: %1;}").arg(strBody);
  css += QString("#reconStatus {font-size: 13px; color: %1;}").arg(strBody);
  css += QString("#reconRule {background: %1;}").arg(strRule);
  css += QString("#reconCard {background: %1; border: 1px solid %2; border-radius: 8px;}").arg(strCardBg, strCardBd);
  css += QString("QLineEdit {background: %1; border: 1px solid %2; border-radius: 4px; padding: 2px 6px;"
                 " color: %3; font-size: 13px; min-height: 22px; max-height: 22px;"
                 " selection-background-color: %7; selection-color: white;}"
                 "QLineEdit:focus {border: 1px solid %4;}"
                 "QLineEdit:disabled {color: %5; border: 1px solid %6;}")
         .arg(strFieldBg, strFieldBd, strBody, strFocus, strCaption, strRule, strSel);
  // A stylesheet font rule replaces the setFont() font, so the mono family is
  // named here as well as the size (style.qss sets 18px app-wide).
  const QString strMono = "font-family: \"Menlo\", \"DejaVu Sans Mono\", \"Noto Sans Mono\", \"Liberation Mono\", monospace;";
  css += "#reconPreview, #reconLog, #reconInventory, #reconElapsed, #reconStatusDetail {font-size: 12px; " + strMono + "}";
  css += "#reconMonoField {font-size: 13px; " + strMono + "}";
  css += "SegmentedControl {font-size: 13px;}";
  css += QString("#reconPreview, #reconLog {color: %1; selection-background-color: %2; selection-color: white;}")
         .arg(strBody, strSel);
  css += QString("#reconProgress {background: %1; border: none; border-radius: 2px;}"
                 "#reconProgress::chunk {background: %2; border-radius: 2px;}").arg(strTrack, QString(kPurple));
  // The standard primary pill (pushButtonCC's geometry) in purple; Stop uses
  // the Correction red.
  css += QString("QPushButton[reconKind=\"primary\"] {background: %3; padding: 0 18px; font-size: 14px; font-weight: bold;"
                 " border: none; border-radius: 18px; color: white; min-height: 36px; max-height: 36px;}"
                 "QPushButton[reconKind=\"primary\"]:hover {background: %4;}"
                 "QPushButton[reconKind=\"primary\"]:pressed {background: %5;}"
                 "QPushButton[reconKind=\"primary\"]:disabled, QPushButton[reconKind=\"danger\"]:disabled"
                 " {background: %1; color: %2;}"
                 "QPushButton[reconKind=\"danger\"] {background: #a2293c; padding: 0 18px; font-size: 14px;"
                 " border: none; border-radius: 18px; color: white; min-height: 36px; max-height: 36px;}"
                 "QPushButton[reconKind=\"danger\"]:hover {background: #8c2334;}"
                 "QPushButton[reconKind=\"danger\"]:pressed {background: #761d2b;}")
         .arg(strDisBg, strDisRun, QString(kPurple), QString(kPurpleHover), QString(kPurplePressed));
  // The ghost pill of the other pages' Back / Undo buttons
  // (MainWindow::RethemeGhostPillButtons), at control height.
  css += QString("QPushButton[reconKind=\"ghost\"] {background: transparent; padding: 0 12px; font-size: 13px;"
                 " border: 1px solid %1; border-radius: 14px; color: %2; min-height: 26px; max-height: 26px;}"
                 "QPushButton[reconKind=\"ghost\"]:hover {background: %3; color: %4;}"
                 "QPushButton[reconKind=\"ghost\"]:disabled {color: %5; border: 1px solid %6;}")
         .arg(strGhostBd, strGhostTx, strGhostHv, strBody, strDisTx, strDisBd);
  css += QString("QPushButton[reconKind=\"disclosure\"] {background: transparent; border: none; padding: 4px 0;"
                 " font-size: 13px; color: %1; text-align: left;}"
                 "QPushButton[reconKind=\"disclosure\"]:hover {color: %2;}").arg(strCaption, strBody);

  // Wide enough for its longest label, so Run -> Stop -> Try again never
  // resizes it (or shifts the footer).
  {
    QFont f = m_pRun->font();
    f.setPixelSize(14);
    const QFontMetrics fm(f);
    int w = 0;
    foreach (const QString& s, QStringList() << "Run Reconstruction" << "Stop" << QString::fromUtf8("Stopping\xe2\x80\xa6") << "Try again")
      w = qMax(w, fm.horizontalAdvance(s));
    m_pRun->setMinimumWidth(w + 40);
  }

  setStyleSheet(css);
  // The pane is no longer our child once MainWindow has inserted it, so it
  // does not inherit the sheet above; it gets its own copy.
  if (m_pPane)
    m_pPane->setStyleSheet(css);
}

// ==== ReconRunner ====

namespace
{
// One entry per project; parented to the app, so they outlive every window.
QMap<QString, ReconRunner*>& Runners()
{
  static QMap<QString, ReconRunner*> map;
  return map;
}

QString RootKey(const QString& strRoot)
{
  if (strRoot.isEmpty())
    return QString();
  const QString strCanonical = QFileInfo(strRoot).canonicalFilePath();
  return strCanonical.isEmpty() ? QDir::cleanPath(QFileInfo(strRoot).absoluteFilePath()) : strCanonical;
}
}

ReconRunner::ReconRunner(const QString& strRoot, QObject* parent) : QObject(parent), m_strRoot(strRoot)
{
  // A tool that ignores SIGTERM must not stay on "Stopping".
  m_pKill = new QTimer(this);
  m_pKill->setSingleShot(true);
  m_pKill->setInterval(5000);
  connect(m_pKill, &QTimer::timeout, this, [this]()
  {
    if (m_proc && m_proc->state() != QProcess::NotRunning)
      m_proc->kill();
  });
}

ReconRunner* ReconRunner::ForRoot(const QString& strRoot)
{
  return Runners().value(RootKey(strRoot), nullptr);
}

ReconRunner* ReconRunner::Start(const QString& strRoot, const QString& strProgram, const QStringList& args,
                                const QString& strCommandLine)
{
  const QString strKey = RootKey(strRoot);
  ReconRunner* old = Runners().take(strKey);
  if (old)
  {
    // Only a finished run is replaced (the page starts none while one is active).
    old->disconnect();
    if (old->m_proc)
      old->m_proc->disconnect(old);
    old->deleteLater();
  }
  ReconRunner* r = new ReconRunner(strRoot, qApp);
  Runners().insert(strKey, r);
  r->m_strProgram = strProgram;
  r->m_listLines << "$ " + strCommandLine;
  r->m_proc = new QProcess(r);
  connect(r->m_proc, &QProcess::readyReadStandardOutput, r, [r]()
  {
    r->AppendChunk(QString::fromLocal8Bit(r->m_proc->readAllStandardOutput()), r->m_strBufOut);
  });
  connect(r->m_proc, &QProcess::readyReadStandardError, r, [r]()
  {
    r->AppendChunk(QString::fromLocal8Bit(r->m_proc->readAllStandardError()), r->m_strBufErr);
  });
  connect(r->m_proc, static_cast<void (QProcess::*)(int, QProcess::ExitStatus)>(&QProcess::finished),
          r, &ReconRunner::OnProcFinished);
  connect(r->m_proc, &QProcess::errorOccurred, r, &ReconRunner::OnProcError);
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  env.insert("PYTHONUNBUFFERED", "1");   // line-by-line streaming; changes buffering only
  r->m_proc->setProcessEnvironment(env);
  r->m_timer.start();
  // start(program, args): the command-line overload no longer splits on Qt 6,
  // and argv needs no quoting.
  r->m_proc->start(strProgram, args);
  return r;
}

QStringList ReconRunner::RunningRoots()
{
  QStringList l;
  foreach (ReconRunner* r, Runners())
    if (r->IsActive())
      l << r->m_strRoot;
  return l;
}

bool ReconRunner::ConfirmQuit(QWidget* parent)
{
  const QStringList roots = RunningRoots();
  if (roots.isEmpty())
    return true;
  QStringList names;
  foreach (const QString& strRoot, roots)
  {
    ProjectInfo info;
    names << (ProjectManager::LoadProject(strRoot, info) && !info.name.isEmpty() ? info.name : QFileInfo(strRoot).fileName());
  }
  QMessageBox box(parent);
  box.setWindowTitle("Stop the reconstruction?");
  box.setText(roots.size() == 1 ? QString("A reconstruction is still running for %1.").arg(names.first())
                                : QString("Reconstructions are still running for %1.").arg(names.join(", ")));
  box.setInformativeText("Quitting stops it, and nothing is saved. To keep it running, leave the app open; "
                         "you can work on other projects meanwhile (File > Switch Project).");
  QPushButton* stop = box.addButton("Stop and quit", QMessageBox::DestructiveRole);
  QPushButton* keep = box.addButton("Keep running", QMessageBox::RejectRole);
  box.setDefaultButton(keep);
  box.exec();
  return box.clickedButton() == stop;
}

void ReconRunner::StopAll()
{
  // Terminate, then kill, then give up. Disconnected first: nothing is shown any more.
  foreach (ReconRunner* r, Runners())
  {
    if (!r->m_proc || r->m_proc->state() == QProcess::NotRunning)
      continue;
    r->m_proc->disconnect(r);
    r->m_proc->terminate();
    if (!r->m_proc->waitForFinished(3000))
    {
      r->m_proc->kill();
      r->m_proc->waitForFinished(1000);
    }
    // No longer running: nothing treats the project as reconstructing.
    r->m_bStopRequested = true;
    r->Finish(Stopped);
  }
}

void ReconRunner::Stop()
{
  if (!IsActive() || !m_proc || m_proc->state() == QProcess::NotRunning)
    return;
  m_bStopRequested = true;
  m_proc->terminate();
  m_pKill->start();
}

// Whole lines only: chunks can end mid-line. A carriage return (progress
// counters) ends a line too.
void ReconRunner::AppendChunk(const QString& strChunk, QString& strBuffer)
{
  strBuffer += strChunk;
  int start = 0;
  for (int i = 0; i < strBuffer.size(); ++i)
  {
    const QChar c = strBuffer.at(i);
    if (c != '\n' && c != '\r')
      continue;
    const QString line = strBuffer.mid(start, i - start);
    if (!(c == '\r' && line.isEmpty()))
    {
      m_listLines << line;
      emit LineAdded(line);
      NoteLine(line);
    }
    if (c == '\r' && i + 1 < strBuffer.size() && strBuffer.at(i + 1) == '\n')
      ++i;
    start = i + 1;
  }
  strBuffer = strBuffer.mid(start);
}

// The command's own progress lines, turned into the phase caption.
void ReconRunner::NoteLine(const QString& strLine)
{
  if (!IsActive() || m_bStopRequested)
    return;
  QString phase;
  if (strLine.startsWith("Adjusting settings")) phase = "Adjusting settings";
  else if (strLine.startsWith("Summary of settings") || strLine.startsWith("Optimization")) phase = "Reading photos";
  else if (strLine.contains("Working on resolution")) phase = "Aligning slabs";
  else if (strLine.startsWith("Machine learning interpolation")) phase = "Filling in slices";
  else if (strLine.startsWith("Writing deformed references") || strLine.startsWith("Writing uncorrected photo recon"))
    phase = "Writing outputs";
  if (!phase.isEmpty())
  {
    m_strPhase = phase;
    emit PhaseChanged(phase);
  }
}

void ReconRunner::OnProcFinished(int nExitCode, QProcess::ExitStatus status)
{
  if (!m_strBufOut.isEmpty()) AppendChunk("\n", m_strBufOut);
  if (!m_strBufErr.isEmpty()) AppendChunk("\n", m_strBufErr);
  m_nExitCode = nExitCode;
  m_bCrashed = (status == QProcess::CrashExit);
  if (m_bStopRequested)
    Finish(Stopped);
  else
    Finish(status == QProcess::NormalExit && nExitCode == 0 ? Succeeded : Failed);
}

void ReconRunner::OnProcError(QProcess::ProcessError err)
{
  if (err != QProcess::FailedToStart)
    return;   // every other error also arrives as finished()
  m_bFailedToStart = true;
  const QString line = "Could not start " + m_strProgram;
  m_listLines << line;
  emit LineAdded(line);
  Finish(Failed);
}

void ReconRunner::Finish(Outcome o)
{
  if (!IsActive())
    return;
  m_pKill->stop();
  m_nFinalMs = m_timer.elapsed();
  m_outcome = o;
  emit Finished();
}
