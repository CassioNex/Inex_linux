#include <QtWidgets>
#include "window.h"
#include "finger_img.h"

Window::Window(QWidget *parent) : QWidget(parent)
{
    QLabel *logo = new QLabel;
    logo->setAlignment(Qt::AlignHCenter);
    logo->setPixmap(QPixmap(":/images/nice_logo_mono.bmp"));

    QLabel *title = new QLabel;
    QLabel *subtitle = new QLabel;
    title->setText("Biometrics App");
    subtitle->setText("v1.0");
    QFont f = title->font();
    f.setPointSize(12);
    f.setBold(true);
    title->setFont(f);
    title->setAlignment(Qt::AlignHCenter);
    subtitle->setFont(f);
    subtitle->setAlignment(Qt::AlignHCenter);

    finger_img = new FingerImg;
    QVBoxLayout *layout = new QVBoxLayout;
    QHBoxLayout *buttons = new QHBoxLayout;

    QPushButton *bt1 = new QPushButton("Enroll");
    QPushButton *bt2 = new QPushButton("Match");

    QObject::connect(bt1, &QPushButton::clicked, this, &Window::bt1ButtonClicked);
    QObject::connect(bt2, &QPushButton::clicked, this, &Window::bt2ButtonClicked);

    buttons->addWidget(bt1);
    buttons->addWidget(bt2);

    layout->addWidget(logo);
    layout->addWidget(title);
    layout->addWidget(subtitle);
    layout->addWidget(finger_img);
    layout->addLayout(buttons);

    setLayout(layout);

    setWindowTitle(tr("Fingerprint capture"));

    // Abre descriptors para os leds de enrol e match
    m_fd_led_enrol = open(LED_B_BRIGHTNESS, O_WRONLY);
    if (m_fd_led_enrol < 0) {
        exit(1);
    }
    m_fd_led_match = open(LED_Y_BRIGHTNESS, O_WRONLY);
    if (m_fd_led_match < 0) {
        exit(1);
    }
    // Default - enrol mode
    write(m_fd_led_enrol, "1", 2);
    write(m_fd_led_match, "0", 2);
}

void Window::handleImgCaptured(uint8_t *buf)
{
    finger_img->handleImgCaptured(buf);
}

void Window::bt1ButtonClicked()
{
    //finger_img->hide();
    finger_img->m_enroll_match = ENROL_MODE;
    write(m_fd_led_enrol, "1", 2);
    write(m_fd_led_match, "0", 2);
}

void Window::bt2ButtonClicked()
{
    //finger_img->show();
    finger_img->m_enroll_match = MATCH_MODE;
    write(m_fd_led_enrol, "0", 2);
    write(m_fd_led_match, "1", 2);
}
