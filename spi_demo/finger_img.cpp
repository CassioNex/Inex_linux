#include <QtWidgets>
#include "finger_img.h"

FingerImg::FingerImg(QWidget *parent) : QWidget(parent)
{
    /* Debug
    int i, j;

    // Inicializa buffer de imagem
    for (j = 0; j < N_LIN; j++) {
        for (i = 0; i < N_COL; i++) {
            img_array[j][i] = 0;
        }
    }
    */

    QVBoxLayout *layout = new QVBoxLayout;

    QLabel *title = new QLabel(tr("Finger Image"));
    title->setAlignment(Qt::AlignHCenter);

    finger = new QLabel();
    finger->setAlignment(Qt::AlignHCenter);
    finger->setPixmap(QPixmap(N_COL, N_LIN));

    layout->addWidget(title);
    layout->addWidget(finger);

    setLayout(layout);
}

void FingerImg::handleImgCaptured(uint8_t *buf)
{
    uint8_t img_buf[IMG_SIZE];

    /* Debug
    volatile int i, j, k;
    // Transfere imagem para array (x,y)
    for (j = 0; j < N_LIN; j++) {
        for (i = 0; i < N_COL; i++) {
            img_array[j][i] = *buf++;
        }
    }
    */

    // Transfere imagem para img_buf
    for (int i = 0; i < IMG_SIZE; i++) {
        img_buf[i] = *buf++;
    }

    // Carrega buffer de imagem em objeto QImage,
    // e plota via QPixmap no label finger
    QImage img(img_buf, N_COL, N_LIN, QImage::Format_Grayscale8);
    finger->setPixmap(QPixmap::fromImage(img));
}
