#include <QtWidgets>
#include "window.h"

// Construtora
Window::Window(QWidget *parent) : QWidget(parent)
{
    // Instancia logo e título da app
    QLabel *logo = new QLabel;
    logo->setAlignment(Qt::AlignHCenter);
    logo->setPixmap(QPixmap(":/images/nice_logo_mono.bmp"));

    QLabel *title = new QLabel;
    title->setText("Bio App v2.0");
    QFont f = title->font();
    f.setPointSize(12);
    f.setBold(true);
    title->setFont(f);
    title->setAlignment(Qt::AlignHCenter);

    // Instancia objetos dp, finger_img etc
    dp = new DataProvider;
    finger_img = new FingerImg;
    cb_mode   = new QCheckBox("Mode");
    bt_enroll = new QPushButton("Enroll");
    bt_search = new QPushButton("Search");
    timer     = new QElapsedTimer;
    lb_enroll_cnt = new QLabel;

    // Conecta signal de imagem capturada ao slot de tratamento
    QObject::connect(dp, &DataProvider::img_captured, this, &Window::handleImgCaptured);
    // Demais connects
    QObject::connect(cb_mode, &QCheckBox::clicked, this, &Window::cb_mode_clicked);
    QObject::connect(bt_enroll, &QPushButton::clicked, this, &Window::bt_enroll_ButtonClicked);
    QObject::connect(bt_search, &QPushButton::clicked, this, &Window::bt_search_ButtonClicked);

    // Define os layouts
    QVBoxLayout *layout = new QVBoxLayout;
    QHBoxLayout *buttons = new QHBoxLayout;
    QHBoxLayout *h_lay   = new QHBoxLayout;

    buttons->addWidget(bt_enroll);
    buttons->addWidget(bt_search);

    h_lay->addWidget(cb_mode);
    h_lay->addWidget(lb_enroll_cnt);

    layout->addWidget(logo);
    layout->addWidget(title);
    layout->addLayout(h_lay);
    layout->addWidget(finger_img);
    layout->addLayout(buttons);

    setLayout(layout);

    setWindowTitle(tr("Fingerprint capture"));

    // Inicia checkbox de modo como checked (enroll)
    cb_mode->setChecked(true);
    bt_enroll->setEnabled(true);
    bt_search->setEnabled(false);
}

// Slots
void Window::handleImgCaptured(uint8_t *buf)
{
    timer->start();
    finger_img->handleImgCaptured(buf);
    qDebug() << "HandleImgCapture() levou " << timer->elapsed() << "milliseconds";
    // Atualiza o label de templates aprendidos
    lb_enroll_cnt->setText(QString::number(finger_img->enroll_cnt));
}

void Window::cb_mode_clicked()
{
    if (cb_mode->checkState()) {
        finger_img->m_mode = ENROLL_MODE;
        bt_enroll->setEnabled(true);
        bt_search->setEnabled(false);
    }
    else {
        finger_img->m_mode = SEARCH_MODE;
        bt_enroll->setEnabled(false);
        bt_search->setEnabled(true);
    }
}

void Window::bt_enroll_ButtonClicked()
{
    dp->capture_img();

}

void Window::bt_search_ButtonClicked()
{
    dp->capture_img();
}
