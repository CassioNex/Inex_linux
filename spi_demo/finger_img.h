#ifndef FINGER_IMG_H
#define FINGER_IMG_H

#include <QWidget>
#include "data_provider.h"

class QLabel;

class FingerImg : public QWidget
{
    Q_OBJECT

public:
    explicit FingerImg(QWidget *parent = nullptr);
    /* Debug - Buffer (x,y) final - imagem
    uint8_t img_array[N_LIN][N_COL];
    */

signals:

public slots:
    void handleImgCaptured(uint8_t *buf);

private:
    QLabel *finger;
};

#endif // FINGER_IMG_H
