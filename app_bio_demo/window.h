#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <unistd.h>
#include "data_provider.h"

// Leds
#define LED_B_BRIGHTNESS    "/sys/class/leds/led_b/brightness"
#define LED_Y_BRIGHTNESS    "/sys/class/leds/led_y/brightness"

// Forward declaration
class FingerImg;

class Window : public QWidget
{
    Q_OBJECT

public:
    explicit Window(QWidget *parent = nullptr);

signals:

public slots:
    void handleImgCaptured(uint8_t *buf);

private slots:
    void bt1ButtonClicked();
    void bt2ButtonClicked();

private:
     FingerImg *finger_img;
     int m_fd_led_enrol;
     int m_fd_led_match;
};

#endif // WINDOW_H
