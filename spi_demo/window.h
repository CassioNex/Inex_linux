#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include "data_provider.h"

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
};

#endif // WINDOW_H
